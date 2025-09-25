library(sf)
library(ggplot2)
library(data.table)
library(lubridate)
library(RColorBrewer)

# Load spatial data
COMID_df <- read.csv("NatalStream_COMIDs.csv")
cids <- unique(COMID_df$COMID)
salmon_snake.streams <- read_sf("salmon_snake.streams.shp")
study_reaches <- subset(salmon_snake.streams, COMID %in% cids)
streams <- unique(study_reaches$GNIS_NAME)

ggplot()+
  geom_sf(data = salmon_snake.streams, col = "gray")+
  geom_sf(data = study_reaches, col = "red")+
  theme_light()+
  theme(legend.position = "right")

# Load predicted stream temperature data
pred_dir <- "G:/My Drive/daily-st-pnw/predictions"
files <- dir(pred_dir)[grep(170602, dir(pred_dir))]
stdata <- NULL
for(h in files){
  td <- fread(paste0(pred_dir, "/", h))[, c("COMID", "tim.date", "prd.stream_temp")]
  stdata <- rbind.data.frame(stdata, td)
}
stdata$years <- lubridate::year(stdata$tim.date)
stdata$doy <- lubridate::yday(stdata$tim.date)
fwrite(stdata, "salmon_river_stream_temps.csv")

# Cumulative exposure (for temperature, this is degrees-days)
fnc_Cum.Exposure <- function(frame, site, start.date, end.date, st.col, site.col = "COMID", date.col = "tim.date", show.na = T){
  frame <- as.data.frame(frame)
  frame <- frame[frame[, site.col] == site,]
  frame <- frame[frame[, date.col] >= start.date & frame[, date.col] <= end.date,]
  frame <- frame[order(frame[, date.col]),]
  datelist <- sort(unique(frame[, date.col]))
  na.list <- which(is.na(frame[, st.col]))
  na.datelist <- unique(frame[, date.col][na.list])
  
  if(show.na){cat("There are", length(na.datelist),"missing days  \n", as.character(na.datelist), "\n")}
  
  if(any(!is.na(frame[, st.col]))){ #proceed only if there are data
    cumexp <- aggregate(frame[,st.col], list(frame$COMID), sum, na.rm = T)
    colnames(cumexp) <- c("COMID", st.col)
    
    return(cumexp[, st.col])
  } else{ return (NA)}
}

# Proportion of days with temperatures within temperature range XX to YY during facet period
fnc_pDays.in.Range = function(frame, site, start.date, end.date, XX, YY, st.col, site.col = "COMID", date.col = "tim.date", show.na = T){
  frame <- as.data.frame(frame)
  frame <- frame[frame[, site.col] == site,]
  frame <- frame[frame[, date.col] >= start.date & frame[, date.col] <= end.date,]
  frame <- frame[order(frame[, date.col]),]
  datelist <- sort(unique(frame[, date.col]))
  na.list <- which(is.na(frame[, st.col]))
  na.datelist <- unique(frame[, date.col][na.list])
  
  if(show.na){cat("There are", length(na.datelist),"missing days  \n", as.character(na.datelist), "\n")}
  
  period <- as.numeric(end.date - start.date) + 1
  
  if(any(!is.na(frame[, st.col]))){ #proceed only if there are data
    return(nrow(frame[frame[,st.col] >= XX & frame[,st.col] <= YY,]) / period)
  }
}

# Proportion of days with values above or below a threshold 
fnc_pDays.Unsuitable <- function(frame, site, start.date, end.date, XX, st.col, site.col = "COMID", date.col = "tim.date", sign = "GT", show.na = T){
  
  frame <- as.data.frame(frame)
  frame <- frame[frame[, site.col] == site,]
  frame <- frame[frame[, date.col] >= start.date & frame[, date.col] <= end.date,]
  frame <- frame[order(frame[, date.col]),]
  datelist <- sort(unique(frame[, date.col]))
  na.list <- which(is.na(frame[, st.col]))
  na.datelist <- unique(frame[, date.col][na.list])
  
  if(show.na){cat("There are", length(na.datelist),"missing days  \n", as.character(na.datelist), "\n")}
  
  period <- as.numeric(end.date - start.date) + 1
  
  if(any(!is.na(frame[, st.col]))){ #proceed only if there are data
    if(sign == "GT"){
      dayslist <- frame[!is.na(frame[, st.col]) & frame[, st.col] >= XX,]
    } else if(sign == "LT"){
      dayslist <- frame[!is.na(frame[, st.col]) & frame[, st.col] <= XX,]
    }
    
    return(length(unique(dayslist[, date.col]))/period)
  } else{ return(NA)}
}


# Calculate thermal metrics
years <- unique(stdata$years)
temp.lethal <- 23
temp.opt.min <- 10
temp.opt.max <- 17

metrics <- NULL
for(s in streams){
  for(y in years){
    # parameters
    stream_cids <- study_reaches$COMID[study_reaches$GNIS_NAME %in% s]
    td <- stdata[stdata$COMID %in% stream_cids,]
    stdate <- as.Date(paste0(y, "-01-01"))
    endate <- as.Date(paste0(y, "-12-31"))

    # compute
    ce <- lapply(X = unique(td$COMID), FUN = fnc_Cum.Exposure, frame = td, start.date = stdate, end.date = endate, st.col = "prd.stream_temp", show.na = F)
    ce[sapply(ce, is.null)] <- NA; ce <- unlist(ce)
    ps <- lapply(X = unique(td$COMID), FUN = fnc_pDays.in.Range, frame = td, start.date = stdate, end.date = endate, XX = temp.opt.min, YY = temp.opt.max, st.col = "prd.stream_temp", show.na = F)
    ps[sapply(ps, is.null)] <- NA; ps <- unlist(ps)
    pu <- lapply(X = unique(td$COMID), FUN = fnc_pDays.Unsuitable, frame = td, start.date = stdate, end.date = endate, XX = temp.lethal, st.col = "prd.stream_temp", sign = "GT", show.na = F)
    pu[sapply(pu, is.null)] <- NA; pu <- unlist(pu)
    
    # combine
    mets <- cbind.data.frame("stream" = s, "year" = y, "COMID" = unique(td$COMID), "cumulative_exposure" = ce, "prop_days_suitable" = ps, "prop_days_unsuitable" = pu)
    metrics <- rbind.data.frame(metrics, mets)
  }  
}
fwrite(metrics, "salmon_river_thermal_metrics.csv")

# Plot

y <- 2010
results <- dplyr::left_join(study_reaches, metrics[metrics$year %in% y,], by = "COMID")
colnames(results)

fncPlotMap <- function(var){
  plot_prd <- ggplot()+
    geom_sf(data = salmon_snake.streams, col = "gray")+
    geom_sf(data = results, aes(color = {{var}}), lwd = 2)+
    scale_color_gradientn(
      colors = rev(brewer.pal(n = 11, name = "Spectral")), 
      na.value = "gray80")+
    
    theme_light()+
    theme(legend.position = "right")
  return(plot_prd)
}
(p1 <- fncPlotMap(cumulative_exposure))
(p2 <- fncPlotMap(prop_days_suitable))
(p3 <- fncPlotMap(prop_days_unsuitable)) # none
