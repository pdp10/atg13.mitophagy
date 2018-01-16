# The MIT License
# 
# Copyright (c) 2017 Piero Dalle Pezze
# 
# Permission is hereby granted, free of charge, 
# to any person obtaining a copy of this software and 
# associated documentation files (the "Software"), to 
# deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, 
# merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom 
# the Software is furnished to do so, 
# subject to the following conditions:
# 
# The above copyright notice and this permission notice 
# shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(Hmisc)
library(scales)
library(plyr)


# A simple theme without grid
theme_basic <- function(base_size = 28){
  theme_bw(base_size) %+replace%
    theme(
      aspect.ratio=0.6,
      panel.border = element_rect(colour = "black", fill=NA, size=1.0),
      panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank()
    )
}


# https://uk.mathworks.com/help/stats/lognrnd.html?requestedDomain=true 

meanlog <- function(mu, v) {
  log((mu^2)/sqrt(v+mu^2))
}


sdlog <- function(mu, v) {
  sqrt(log(v/(mu^2)+1))
}


skewness <- function(x, na.rm = FALSE, ...) {
  if (na.rm) x = x[!is.na(x)]
  return(sum((x-mean(x))^3/length(x))/sqrt(var(x))^3)
}


kurtosis <- function(x, na.rm = FALSE, ...) {
  if (na.rm) x = x[!is.na(x)]
  return(sum((x-mean(x))^4/length(x))/sqrt(var(x))^4)
}


# return the linear model equation
equation <- function(x, digits=4, show.r2=FALSE) {
  lm_coef <- list(a = round(coef(x)[1], digits=digits),
                  b = round(coef(x)[2], digits=digits),
                  r2 = round(summary(x)$r.squared, digits=digits));
  if(show.r2) {
    if(lm_coef$b >= 0) {
      lm_eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)    
    } else {
      lm_coef$b <- abs(lm_coef$b)
      lm_eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)
    }
  } else {
    if(lm_coef$b >= 0) {
      lm_eq <- substitute(italic(y) == a + b %.% italic(x),lm_coef)    
    } else {
      lm_coef$b <- abs(lm_coef$b)
      lm_eq <- substitute(italic(y) == a - b %.% italic(x),lm_coef)
    }
  }
  return(as.character(as.expression(lm_eq)))
}


plot_correlation <- function(df, xlab, ylab) {
  #print(df)
  pcorr <- cor(df[,c(xlab)], df[,c(ylab)])  
  #print(pcorr)
  g <- ggplot(df, aes_string(xlab, ylab)) +
    geom_point() +     # Use hollow circles
    geom_smooth(method=lm) +  # Add linear regression line (by default includes 95% confidence region)
    ggtitle(paste0("R=", round(pcorr, digits = 3))) + theme_basic()
}


plot_correlation_main <- function(dfnt, points, xlab, ylab, filenameout, location.results) {
  # Correlations for profiles at their initial time point
  time <- data.frame(x=points*10, y=t(dfnt[1,]))
  #print(time)
  colnames(time) <- c('time', ylab)
  #print(time)
  
  g <- plot_correlation(time, 'time', ylab) + 
    labs(x=xlab, y=ylab)
  ggsave(file.path(location.results, paste0(filenameout, "_corr.png")),
         dpi=300,  width=6, height=4)
}





qqplot <- function(vec, distribution=stats::qnorm) {
  d <- data.frame(resids=vec)
  g <- ggplot(d, aes(sample=resids)) + 
    stat_qq(distribution=distribution) + 
    stat_qq_line(distribution=distribution) +
    theme_basic()
  g
}




# overlay histogram, normal and lognormal densities
hist_w_dist <- function(vec) {
  df <- data.frame(x=vec)
  
  g <- ggplot(df, aes(x)) +
    geom_histogram(aes(y = ..density..), binwidth = 1, color = "black", fill = 'white') +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(df$x), sd = sd(df$x)), 
                  lwd = 1,
                  aes(colour='norm')) +
    stat_function(fun = dlnorm, 
                  args = list(meanlog = meanlog(mean(df$x), var(df$x)), sdlog = sdlog(mean(df$x), var(df$x))), 
                  lwd = 1, 
                  aes(colour='lnorm')) +
    scale_colour_manual('distrib', values=c("blue", "red")) +
    theme_basic() + 
    theme(legend.position="top", legend.title=element_blank(), legend.text=element_text(size=22))
  g
}


normality_analysis <- function(data, title, filenameout, location.results) {
  
  q1 <- qqplot(data, stats::qnorm) + ggtitle(paste0('norm'))
  q2 <- qqplot(data, stats::qlnorm) + ggtitle(paste0('lnorm'))
  q3 <- hist_w_dist(data) + labs(x='sample')
  
  plots <- list(plots=c(list(q1), list(q2), list(q3)))
  plot.arrange <- do.call(grid.arrange, c(plots$plots, ncol=3, top=title))
  ggsave(file.path(location.results, paste0(filenameout, '_distrib_analysis.png')), plot=plot.arrange, width=16, height=4, dpi=300)
}


compute_time_vec <- function(time0_pos, length, step=1) {
  pre <- seq(time0_pos-1, 1)*(-step)
  post <- seq(1, length-time0_pos)*step
  c(pre, 0, post)
}



plot_curves <- function(df) {
  mdf <- melt(df,id.vars="time",variable.name="rep",value.name="intensity")
  g <- ggplot() + geom_line(data=mdf,aes(x=time,y=intensity,color=rep), size=1.0) + 
    theme_basic() + theme(legend.position="none")
  return(g)
}


plot_curves2 <-function(df, xlab) {
  names <- colnames(df)
  # remove xlab from names
  names <- names[! names %in% c(xlab)]
  g <- ggplot(df, aes_string(x=xlab)) + 
    theme_basic() + theme(legend.position="none")
  for(i in names) {
    g <- g + geom_line(aes_string(y=i), size=0.2)
  }
  g
}


synchronise_timecourse <- function(dfnt, init_offsets_times, latest_peak_time, ylab, location.data, filenameout, location.results) {
  # Calculate the maximum offset (this line will be the synchronised point)
  max_offset <- max(init_offsets_times)
  #print(max_offset)
  
  # Calculate the number of offset (NAs) to add at the bottom 
  bottom_offset <- max_offset - init_offsets_times
  #print(bottom_offset)
  
  # Create the new data frame adding offset at the top and bottom
  sync_df <- data.frame(na=c(rep(NA, max_offset + nrow(dfnt))))
  for(i in 1:ncol(dfnt)) {
    sync_col_i <- c(rep(NA,init_offsets_times[i]), dfnt[,i], rep(NA, bottom_offset[i]))
    sync_df[ ,i] <- sync_col_i
  }
  # copy column names 
  colnames(sync_df) <- colnames(dfnt)
  
  # remove rows with all NAs
  ind <- apply(sync_df, 1, function(x) all(is.na(x)))
  sync_df <- sync_df[ !ind, ]
  
  #print(sync_df)
  
  # Compute time vector. 0 is the maximum time point.
  #time_vec <- compute_time_vec(latest_peak_time, nrow(sync_df), step=10)
  # Compute time vector. 0 is the first row.
  time_vec <- seq(0, nrow(sync_df)-1) *10
  #print(time_vec)
  sync_df <- cbind(time_vec, sync_df)
  colnames(sync_df)[colnames(sync_df) == 'time_vec'] <- 'time'
  
  #print(sync_df)
  
  # plot synchronised profiles
  g <- plot_curves(sync_df)
  g <- g + 
    geom_vline(xintercept=(latest_peak_time*10)-10, size=1.5, color="black", linetype="dashed") +    
    #geom_vline(xintercept=0, size=1, color="magenta", linetype="dashed") +
    labs(title=paste0("n=", ncol(sync_df)), x="Time [s]", y=ylab) + 
    ggsave(file.path(location.results, paste0(filenameout, ".png")), 
           dpi=300,  width=6, height=4)
  
  write.table(sync_df, file=file.path(location.data, paste0(filenameout,"_sync.csv")), sep=",", row.names=FALSE)
  return(sync_df)
}


comp_mean_error <- function(df, ylab, filenameout, location.results) {
  dfnt <- df[,-1]
  time_vec <- subset(df, select=c('time'))
  
  means <- apply(dfnt, 1, mean, na.rm=TRUE)
  stdevs <- sqrt(apply(dfnt, 1, var,na.rm=TRUE))
  lengths <- apply(dfnt, 1, function(x) length(which(!is.na(x))))  
  stderrs <- stdevs / sqrt(lengths)
  ci95 <- qnorm(.975)*(stderrs)
  
  statdf <- data.frame(a=time_vec, b=means, c=stdevs, d=ci95)
  colnames(statdf)[1] <- 'a'
  #print(statdf)
  
  # plot mean+sd
  g <- ggplot(statdf, aes(x=a, y=b)) + 
    geom_errorbar(aes(ymin=b-c, ymax=b+c)) + 
    geom_line(aes(x=a, y=b), color="black", size=1.0) +
    labs(x="Time [s]", y=ylab) + 
    theme_basic()
  ggsave(file.path(location.results, paste0(filenameout, "_sd.png")), dpi=300,  width=6, height=4)
  
  g <- g +     
    geom_errorbar(aes(ymin=b-d, ymax=b+d), colour="magenta") +
    geom_line(aes(x=a, y=b), color="black", size=1.0) + 
    theme_basic()
  ggsave(file.path(location.results, paste0(filenameout, "_ci95.png")), dpi=300,  width=6, height=4)
  
  colnames(statdf) <- c("time", "mean", "sd", "ci95")
  write.table(statdf, file=file.path(location.results, paste0(filenameout, "_stats.csv")), sep=",", row.names=FALSE)  
}




sync_tc_fun <- function(df, location.data, filenameout, ylab, location.results) {
  
  # discard the first column (time)
  dfnt <- df[,-1]
  #print(dfnt)
  
  # for each column (par:2) get the row containing the max value
  peak_times <- apply(dfnt, 2, which.max)
  #print(peak_times)
  # get the peak values
  peak_intensities <- apply(dfnt, 2, function(x) max(x, na.rm=TRUE))
  #print(peak_intensities)
  
  # Get the maximum latest maximum peak
  latest_peak_time <- max(peak_times)
  #print(latest_peak_time)
  # Create a simple data frame containing the intensities
  df_peaks <- data.frame(t(peak_intensities))
  #print(df_peaks)
  
  # Calculate the number of offset (NAs) to add at the top 
  init_offsets_times <- latest_peak_time - peak_times
  #print(init_offsets_times)
  init_offsets_intensities <- unlist(dfnt[1,])
  #print(init_offsets_intensities)
  
  # calculate correlations and plot
  plot_correlation_main(dfnt, init_offsets_times, "Init. sign. [s]", ylab, paste0(filenameout, "_offsets_times_vs_intensities"), location.results)
  plot_correlation_main(df_peaks, peak_times, "peak times [s]", ylab, paste0(filenameout, "_peak_times_vs_intensities"), location.results)  
  
  # save correlation data on file
  corr_df <- data.frame(init_offsets_times*10, peak_times*10, t(dfnt[1,]))
  corr_df[,4] <- rownames(corr_df)
  colnames(corr_df) <- c('init_offsets_times', 'peaks_times', paste0(ylab, "_at_t0"), 'repeat_id')
  write.table(corr_df, file=file.path(location.results, paste0(filenameout, "_corr.csv")), sep=",", row.names=FALSE)
  # save correlation stats on file
  corr_stats_df <- t(data.frame(mean(init_offsets_times*10), sqrt(var(init_offsets_times*10)), 
                                skewness(init_offsets_times*10), kurtosis(init_offsets_times*10), 
                                mean(peak_times*10), sqrt(var(peak_times*10)),
                                skewness(peak_times*10), kurtosis(peak_times*10),
                                mean(init_offsets_intensities), sqrt(var(init_offsets_intensities)),
                                skewness(init_offsets_intensities), kurtosis(init_offsets_intensities),                                
                                mean(peak_intensities), sqrt(var(peak_intensities)),
                                skewness(peak_intensities), kurtosis(peak_intensities),
                                meanlog(mean(init_offsets_times*10),var(init_offsets_times*10)), sdlog(mean(init_offsets_times*10),var(init_offsets_times*10)), 
                                meanlog(mean(peak_times*10), var(peak_times*10)), sdlog(mean(peak_times*10),var(peak_times*10)),
                                meanlog(mean(init_offsets_intensities), var(init_offsets_intensities)), sdlog(mean(init_offsets_intensities*10),var(init_offsets_intensities)),
                                meanlog(mean(peak_intensities), var(peak_intensities)), sdlog(mean(peak_intensities),var(peak_intensities))))
  rownames(corr_stats_df) <- c('init_offsets_times_mean', 'init_offsets_times_sd', 
                               'init_offsets_times_skew', 'init_offsets_times_kurt', 
                               'peak_times_mean', 'peak_times_sd',
                               'peak_times_skew', 'peak_times_kurt',
                               'init_offsets_intensities_mean', 'init_offsets_intensities_sd',
                               'init_offsets_intensities_skew', 'init_offsets_intensities_kurt',
                               'peak_intensities_mean', 'peak_intensities_sd',
                               'peak_intensities_skew', 'peak_intensities_kurt',
                               'init_offsets_times_meanlog', 'init_offsets_times_sdlog',
                               'peak_times_meanlog', 'peak_times_sdlog',
                               'init_offsets_intensities_meanlog', 'init_offsets_intensities_sdlog',
                               'peak_intensities_meanlog', 'peak_intensities_sdlog')
  write.table(corr_stats_df, file=file.path(location.results, paste0(filenameout, "_corr_stats.csv")), sep=",", row.names=TRUE, col.names=FALSE)
  
  
  # check normality and log-normality
  normality_analysis(init_offsets_times, 'offsets t=0', paste0(filenameout, "_init_offsets_times"), location.results)
  normality_analysis(peak_times, 'peak times', paste0(filenameout, "_peak_times"), location.results)
  normality_analysis(init_offsets_intensities, 'offsets int', paste0(filenameout, "_init_offsets_intensities"), location.results)  
  normality_analysis(peak_intensities, 'peak int', paste0(filenameout, "_peak_intensities"), location.results)
  
  
  # synchronise the time courses by maximum peak
  sync_df <- synchronise_timecourse(dfnt, init_offsets_times, latest_peak_time, ylab, location.data, filenameout, location.results)
  
  comp_mean_error(sync_df, ylab, filenameout, location.results)
}



plot_original_tc <- function(df, filenameout, ylab, location.results) {
  names(df)[1] <- "time"
  g <- plot_curves(df) + 
    labs(title=paste0("n=", ncol(df)), x="Time [s]", y=ylab)
  ggsave(file.path(location.results, paste0(filenameout, "_orig.png")), 
         dpi=300,  width=6, height=4)
}


sync_tc_main <- function(location.data, location.results, csv.file, readout) {
  df <- read.table( file.path(location.data, paste0(csv.file, '.csv')), header=TRUE, na.strings="NA", dec=".", sep="," )
  # remove rows and columns containing only NA
  df <- df[colSums(!is.na(df)) > 0]
  filenameout <- paste0(csv.file)
  # plot original time courses
  plot_original_tc(df, filenameout, readout, location.results)    
  # synchronise time courses
  sync_tc_fun(df, location.data, filenameout, readout, location.results)
}



### -------------------- ###


# return the upper values of the oscillations (the peaks)
hv <- function(vec, thres=0.5, ini.tp=0) {
  vec <- vec[!is.na(vec) & as.numeric(names(vec))>=ini.tp & vec>thres]
  return(vec)
}

# return the lower values of the oscillations
lv <- function(vec, thres=0.5, ini.tp=0) {
  vec <- vec[!is.na(vec) & as.numeric(names(vec))>=ini.tp & vec<thres]
  return(vec)  
}


# extract lower and upper values
extract_min_max <- function(location.data, location.results, filename, thres.hv, thres.lv) {
  
  suffix <- '.csv'
  # import the data
  data <- read.table( file.path(location.data, paste0(filename, suffix)), header=TRUE, na.strings="NA", dec=".", sep=",", row.names=1)
  
  # EXTRACT THE UPPER VALUES
  
  # extract the delays from the time courses
  data.plot.hv <- data.frame(time=numeric(0), val=numeric(0))
  for(j in 1:ncol(data)) {
    col.j <- data[,j]
    names(col.j) <- row.names(data)
    hv.j <- hv(col.j, thres=thres.hv)
    hv.tc.j <- as.numeric(names(hv.j))
    data.plot.hv <- rbind(data.plot.hv, data.frame(time=hv.tc.j, val=hv.j))
  }
  data.plot.hv <- cbind(data.plot.hv, pos=rep('top', nrow(data.plot.hv)))
  
  # plot
  g <- ggplot() + 
    geom_point(data=data.plot.hv, aes(x=time, y=val)) +  
    labs(title='Upper Values', x='Time [s]', y='Intensity [a.u.]') +
    theme_basic(base_size = 16)
  ggsave(file.path(location.results, paste0(filename, '_upper_values.png')), width=4, height=3, dpi=300)
  write.table(data.plot.hv, file=file.path(location.results, paste0(filename, '_upper_values', suffix)), row.names=FALSE, quote=FALSE, sep=',')
  
  # density plot
  g <- ggplot(data.plot.hv, aes(x=val)) + 
    geom_density(colour = "black", fill = "#56B4E9", alpha=0.5) +
    labs(title='Density of upper values', x='signal intensity [a.u.]') +
    theme_basic(base_size = 16)
  ggsave(file.path(location.results, paste0(filename, '_upper_values_density.png')), width=4, height=3, dpi=300)
  
  
  # EXTRACT THE LOWER VALUES
  
  # extract the delays from the time courses
  data.plot.lv <- data.frame(time=numeric(0), val=numeric(0))
  for(j in 1:ncol(data)) {
    col.j <- data[,j]
    names(col.j) <- row.names(data)
    lv.j <- lv(col.j, thres=thres.lv)
    lv.tc.j <- as.numeric(names(lv.j))
    data.plot.lv <- rbind(data.plot.lv, data.frame(time=lv.tc.j, val=lv.j))
  }
  data.plot.lv <- cbind(data.plot.lv, pos=rep('bottom', nrow(data.plot.lv)))
  
  # plot
  g <- ggplot() + 
    geom_point(data=data.plot.lv, aes(x=time, y=val)) +  
    labs(title='Lower values', x='Time [s]', y='Intensity [a.u.]') +
    theme_basic(base_size = 16)
  ggsave(file.path(location.results, paste0(filename, '_lower_values.png')), width=4, height=3, dpi=300)
  write.table(data.plot.lv, file=file.path(location.results, paste0(filename, '_lower_values', suffix)), row.names=FALSE, quote=FALSE, sep=',')
  
  # density plot
  g <- ggplot(data.plot.lv, aes(x=val)) + 
    geom_density(colour = "black", fill = "#56B4E9", alpha=0.5) +
    labs(title='Density of lower values', x='signal intensity [a.u.]') +
    theme_basic(base_size = 16)
  ggsave(file.path(location.results, paste0(filename, '_lower_values_density.png')), width=4, height=3, dpi=300)
  
  
  # write the peaks stats
  
  peaks.stats <- data.frame(type=c('top', 'bottom'), mean=c(mean(data.plot.hv$val), mean(data.plot.lv$val)), 
                            sd=c(sd(data.plot.hv$val), sd(data.plot.lv$val)))
  write.table(peaks.stats, file=file.path(location.results, paste0(filename, '_peaks_stats', suffix)), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=',')
  
}



### -------------------- ###


# Generate data sets for Copasi
generate_copasi_datasets <- function(location.data, filename, location.results, filename.out, observable, treatment) {
  suffix <- '.csv'
  # DATA READING
  data <- read.table(file.path(location.data, paste0(filename, suffix)), header=TRUE, na.strings="NA", dec=".", sep=",")
  
  # MELT
  data.copasi <- melt(data, id.vars=c('time'), value.name = observable)
  data.copasi <- subset(data.copasi, select=-c(variable))
  if(treatment == FALSE) {
    data.copasi <- cbind(data.copasi, trtm=rep(0, nrow(data.copasi)))
  } else {
    data.copasi <- cbind(data.copasi, trtm=rep(1, nrow(data.copasi)))
  }
  
  # DATA WRITING
  write.table(data.copasi, file=file.path(location.results, paste0(filename.out, suffix)), sep="\t", na="", dec=".", row.names=FALSE, quote=FALSE)
  print(paste0('Copasi data set saved in: ', file.path(location.results, paste0(filename.out, suffix))))
}