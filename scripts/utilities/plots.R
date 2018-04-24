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

library(plotly)
library(tools)
library(gplots)
# palette
library(colorRamps) # matlab.like
library(RColorBrewer)




# A simple theme without grid
theme_basic <- function(base_size = 12){
  theme_bw(base_size) %+replace%
    theme(
      panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank()
    )
}





### CORRELATIONS ###



plot_correlation_main <- function(dfnt, points, xlab, ylab, filenameout, location.results) {
  # Correlations for profiles at their initial time point
  time <- data.frame(x=points*10, y=t(dfnt[1,]))
  #print(time)
  colnames(time) <- c('time', ylab)
  #print(time)
  
  g <- plot_correlation(time, 'time', ylab) + 
    theme_basic(base_size=28) + 
    labs(x=xlab, y=ylab)
  ggsave(file.path(location.results, paste0(filenameout, "_corr.png")),
         dpi=300,  width=6, height=4)
}


normality_analysis <- function(data, title, filenameout, location.results) {
  
  q1 <- qqplot(data, stats::qnorm) + ggtitle(paste0('norm')) + theme_basic(base_size=28)
  q2 <- qqplot(data, stats::qlnorm) + ggtitle(paste0('lnorm')) + theme_basic(base_size=28)
  q3 <- hist_w_dist(data) + labs(x='sample') + theme_basic(base_size=28) +
    theme(legend.position="top", legend.title=element_blank(), legend.text=element_text(size=22))
  
  plots <- list(plots=c(list(q1), list(q2), list(q3)))
  plot.arrange <- do.call(grid.arrange, c(plots$plots, ncol=3, top=title))
  ggsave(file.path(location.results, paste0(filenameout, '_distrib_analysis.png')), plot=plot.arrange, width=14, height=4, dpi=300)
}


compute_time_vec <- function(time0_pos, length, step=1) {
  pre <- seq(time0_pos-1, 1)*(-step)
  post <- seq(1, length-time0_pos)*step
  c(pre, 0, post)
}




###-------------###



### TIME COURSES ###



plot_curves <- function(df) {
  mdf <- melt(df,id.vars="time",variable.name="rep",value.name="intensity")
  g <- ggplot() + geom_line(data=mdf,aes(x=time,y=intensity,color=rep), size=1.0) + 
    theme_basic(base_size=28) + theme(legend.position="none")
  return(g)
}


plot_curves2 <-function(df, xlab) {
  names <- colnames(df)
  # remove xlab from names
  names <- names[! names %in% c(xlab)]
  g <- ggplot(df, aes_string(x=xlab)) + 
    theme_basic(base_size=28) + theme(legend.position="none")
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
  g <- plot_curves(sync_df) +
    geom_vline(xintercept=(latest_peak_time*10)-10, size=1.5, color="black", linetype="dashed") +    
    labs(x="Time [s]", y=ylab, title='') +  # title=paste0("n=", ncol(sync_df))) +
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
  
  # plot mean+sd+ci95
  g <- ggplot(statdf, aes(x=a, y=b)) + 
    geom_errorbar(aes(ymin=b-c, ymax=b+c)) + 
    geom_errorbar(aes(ymin=b-d, ymax=b+d), colour="magenta") +
    geom_line(aes(x=a, y=b), color="black", size=1)
  
  colnames(statdf) <- c("time", "mean", "sd", "ci95")
  write.table(statdf, file=file.path(location.results, paste0(filenameout, "_stats.csv")), sep=",", row.names=FALSE)  
  return(g)
}


# plot the mean with error bars
comp_mean_error_w_linear_model <- function(df, filename, ylab, show.linear.model=FALSE) {
  
  dfnt <- df[,-1]
  time_vec <- subset(df, select=c('Time'))
  
  means <- apply(dfnt, 1, mean, na.rm=TRUE)
  stdevs <- sqrt(apply(dfnt, 1, var,na.rm=TRUE))
  lengths <- apply(dfnt, 1, function(x) length(which(!is.na(x))))  
  stderrs <- stdevs / sqrt(lengths)
  ci95 <- qnorm(.975)*(stderrs)
  
  statdf <- data.frame(Time=time_vec, means=means, stdevs=stdevs, ci95=ci95)
  colnames(statdf)[1] <- 'Time'
  #print(statdf)
  
  if(show.linear.model) {
    fit <- lm(means ~ Time, data = statdf)
  }
  
  # plot mean+sd+ci95
  g <- ggplot(statdf, aes(x=Time, y=means)) + 
    # SD
    geom_errorbar(aes(ymin=means-stdevs, ymax=means+stdevs)) + 
    # CI 95%
    geom_errorbar(aes(ymin=means-ci95, ymax=means+ci95), colour="magenta") +
    # mean
    geom_line(aes(x=Time, y=means), color="black", size=1.0) +
    labs(x="Time [s]", y=ylab
         , title=paste0('')
         #, title=paste0('n=', ncol(df)-1)
    ) +
    theme_basic(base_size=40)
  
  if(show.linear.model) {
    x0 <- 550 #500
    y0 <- 500 #-0.1
    g <- g +
      stat_smooth(method = "lm", se=TRUE, color="blue", aes(group=1)) + 
      annotate("text", x=x0, y=y0, label = equation(fit), size=8, parse = TRUE)
    
    # export regression data
    data.regr <- data.frame(regression='meansVStime', slope=coef(fit)[2], intercept=coef(fit)[1], check.names = FALSE)
    write.table(data.regr, file=file.path(paste0(filename, '_linear_regression_data.csv')), row.names=FALSE, quote=FALSE, sep=',')
  }
  
  colnames(statdf) <- c("Time", "mean", "sd", "ci95")
  write.table(statdf, file=file.path(paste0(filename, "_stats.csv")), sep=",", row.names=FALSE)
  return (g)
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
                               'init_offsets_times_skew', 'init_offsets_times_kurt (excess)', 
                               'peak_times_mean', 'peak_times_sd',
                               'peak_times_skew', 'peak_times_kurt (excess)',
                               'init_offsets_intensities_mean', 'init_offsets_intensities_sd',
                               'init_offsets_intensities_skew', 'init_offsets_intensities_kurt (excess)',
                               'peak_intensities_mean', 'peak_intensities_sd',
                               'peak_intensities_skew', 'peak_intensities_kurt (excess)',
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
  

  #################
  # Normality tests
  #################
  swt <- shapiro.test(init_offsets_times)
  kt <- kurtosis.test(init_offsets_times)
  st <- skew.test(init_offsets_times)
  df.stats <- data.frame(test=c('shapiro.test', 'kurtosis.test', 'skew.test'), 
                         pvalue=c(swt$p.value, kt, st))
  write.table(df.stats, file=file.path(location.results, paste0(filenameout, '_init_offsets_times_normality_tests.csv')), row.names=FALSE, quote=FALSE, sep=',')
  
  swt <- shapiro.test(peak_times)
  kt <- kurtosis.test(peak_times)
  st <- skew.test(peak_times)
  df.stats <- data.frame(test=c('shapiro.test', 'kurtosis.test', 'skew.test'), 
                         pvalue=c(swt$p.value, kt, st))
  write.table(df.stats, file=file.path(location.results, paste0(filenameout, '_peak_times_normality_tests.csv')), row.names=FALSE, quote=FALSE, sep=',')
  
  swt <- shapiro.test(init_offsets_intensities)
  kt <- kurtosis.test(init_offsets_intensities)
  st <- skew.test(init_offsets_intensities)
  df.stats <- data.frame(test=c('shapiro.test', 'kurtosis.test', 'skew.test'), 
                         pvalue=c(swt$p.value, kt, st))
  write.table(df.stats, file=file.path(location.results, paste0(filenameout, '_init_offsets_intensities_normality_tests.csv')), row.names=FALSE, quote=FALSE, sep=',')
  
  swt <- shapiro.test(peak_intensities)
  kt <- kurtosis.test(peak_intensities)
  st <- skew.test(peak_intensities)
  df.stats <- data.frame(test=c('shapiro.test', 'kurtosis.test', 'skew.test'), 
                         pvalue=c(swt$p.value, kt, st))
  write.table(df.stats, file=file.path(location.results, paste0(filenameout, '_peak_intensities_normality_tests.csv')), row.names=FALSE, quote=FALSE, sep=',')
  
  
  # synchronise the time courses by maximum peak
  sync_df <- synchronise_timecourse(dfnt, init_offsets_times, latest_peak_time, ylab=paste(ylab, '[a.u.]'), location.data, filenameout, location.results)
  
  g <- comp_mean_error(sync_df, ylab, filenameout, location.results) + 
    labs(x="Time [s]", y=paste(ylab, '[a.u.]'), title='') +     
    theme_basic(base_size=28)
  ggsave(file.path(location.results, paste0(filenameout, "_ci95.png")), dpi=300,  width=6, height=4)
}



plot_original_tc <- function(df, filenameout, ylab, location.results) {
  names(df)[1] <- "time"
  g <- plot_curves(df) + 
    labs(x="Time [s]", y=ylab, title="")  #, title=paste0("n=", ncol(df))) 
  ggsave(file.path(location.results, paste0(filenameout, "_orig.png")), 
         dpi=300,  width=6, height=4)
}


sync_tc_main <- function(location.data, location.results, csv.file, readout) {
  df <- read.table( file.path(location.data, paste0(csv.file, '.csv')), header=TRUE, na.strings="NA", dec=".", sep="," )
  # remove rows and columns containing only NA
  df <- df[colSums(!is.na(df)) > 0]
  filenameout <- paste0(csv.file)
  # plot original time courses
  plot_original_tc(df, filenameout, ylab=paste(readout, '[a.u.]'), location.results)    
  # synchronise time courses
  sync_tc_fun(df, location.data, filenameout, readout, location.results)
}






# Apply a smooth.spline to a data frame
spline.data.frame <- function(data, spar) {
  # create the spline for time
  data.spline <- data.frame(Time=smooth.spline(data[,1], spar=spar)$y)
  # create the splines for each event
  for(i in 2:ncol(data)) {
    sp <- smooth.spline(na.omit(data[,i]), spar=spar)$y
    # extract $y and add NAs
    sp <- c(sp, rep(NA, nrow(data.spline)-length(sp)))
    data.spline <- cbind(data.spline, sp)
  }
  colnames(data.spline) <- colnames(data)
  return(data.spline)
} 


# plot the time course repeats
plot_tc_repeats <- function(df, ylab) {
  df.melt <- melt(df,id.vars=c("Time"), variable.name = 'repeats')
  g <- ggplot() + geom_line(data=df.melt,aes(x=Time,y=value,color=repeats), size=1) +
    labs(x="Time [s]", y=ylab) + #, title=paste0('n=', ncol(df)-1)) +
    theme_basic(base_size=40)
  return(g)
}


# Plot the time course
plot_tc <- function(df, title='title', xlab='Time [s]', ylab='Sign. Int. [a.u.]') {
  g <- ggplot(data=df, aes(x=x, y=y)) + 
    geom_line() + geom_point() +
    labs(title=title, x=xlab, y=ylab) + 
    theme_basic(15)
    #theme_basic(12)
  return(g)
}


plot_tc_with_spline <- function(data, file, spline.order, title, xlab='Time [s]', ylab='Green (Ch2) Sign. Int. [a.u.]') {
  data.spline <- spline(data$X, data$Y, n=spline.order)
  data.spline <- data.frame(X=data.spline$x, Y=data.spline$y)
  
  g <- ggplot() + 
    geom_point(data=data, aes(x=X, y=Y), shape=1) +
    geom_line(data=data.spline, aes(x=X, y=Y), color="red", size=1.0) +
    theme_basic() + 
    labs(title=title, x=xlab, y=ylab)
  
  data.spline$X <- data.spline$X
  write.csv(data.spline, file=file, row.names=FALSE)
  
  return(g)
}


plot_tc_with_regr_line <- function(data, file, title, xlab='Time [s]', ylab='Green (Ch2) Sign. Int. [a.u.]') {
  # linear model for my data
  fit <- lm(Y ~ X, data = data)
  digits <- 2
  
  g <- ggplot(data=data, aes(x=X, y=Y)) + 
    geom_point(shape=1) +
    stat_smooth(method = "lm", se=TRUE, color="red") + 
    annotate("text", x=75, y=640, label = equation(fit, digits=digits), parse = TRUE, size=2.6) +
    theme_basic() + 
    labs(title=title, x=xlab, y=ylab)
  
  # export regression data
  data.regr <- data.frame(file=title, slope=coef(fit)[2], intercept=coef(fit)[1], check.names = FALSE)
  if(!file.exists(file)){
    write.table(data.regr, file=file, row.names=FALSE, quote=FALSE, sep=',')
  } else {
    write.table(data.regr, file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',') 
  }
  return(g)
}



plot_combined_tc <- function(data, expand.xaxis=TRUE) {
  samples <- colnames(data)
  PLOTS <- list() 
  for(i in 2:ncol(data)) {
    if(expand.xaxis) {
      non.na <- sum(!is.na(data[,i]))
      #print(non.na)
      df <- data.frame(x=data[1:non.na,1], y=data[1:non.na,i])
    } else {
      df <- data.frame(x=data[,1], y=data[,i])      
    }
    p <- plot_tc(df, title=samples[i], xlab='Time [s]', ylab='Intensity Mean [a.u.]') +
      theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "in")) # , axis.text.x=element_text(angle=-45, hjust=0, vjust=0.8))
    PLOTS <- c(PLOTS, list(p))
  }
  return( list(plots=PLOTS) )
}



# plot the time courses
plot_synchronised_tc <- function(df, filename, ylab) {
  # plot time courses
  plot.tc <- plot_tc_repeats(df, ylab=paste(ylab, '[a.u.]'))
  # plot mean, sd, ci95, and save stats for the tc
  plot.tc.err <- comp_mean_error_w_linear_model(df, filename, ylab, show.linear.model=FALSE)
  # plot mean, sd, ci95, and save stats for the tc. Also add linear model information
  plot.tc.err.abline <- comp_mean_error_w_linear_model(df, filename, ylab, show.linear.model=TRUE)
  
  
  # Comment if legend is inserted
  plot.tc <- plot.tc + theme(legend.position="none")
  # Uncomment if legend is inserted. This adds extra margins
  #plot.tc.err <- plot.tc.err + theme(plot.margin=unit(c(0.3,3.8,0.3,0.8),"cm"))
  #plot.tc.err.abline <- plot.tc.err.abline + theme(plot.margin=unit(c(0.3,3.8,0.3,0.8),"cm"))
  
  PLOTS <-list(plots=c(list(g1=plot.tc, g2=plot.tc.err, g3=plot.tc.err.abline)))
  
  plot.arrange <- do.call(grid.arrange, c(PLOTS$plots, nrow=3, bottom=paste0(filename, suffix)))
  # Comment if legend is inserted 
  ggsave(file.path(paste0(filename, '_all.png')), plot=plot.arrange, width=8, height=18, dpi=300)    
  # Uncomment if legend is inserted 
  #ggsave(file.path(paste0(filename, '_all.png')), plot=plot.arrange, width=8, height=18, dpi=300)  
  return(plot.arrange)
}


# Filter the data set
data_filtering <- function(df, remove.cols, remove.row.head, remove.row.tail) {
  # remove the unwanted columns
  df.filt <- df[, !(names(df) %in% c('Time'))]
  df.filt <- df.filt[, !(names(df.filt) %in% remove.cols)]
  
  # remove first and last entries
  df.filt <- df.filt[1:remove.row.tail, ]
  df.filt <- df.filt[remove.row.head:nrow(df.filt), ]  
  
  # log my data.
  # FluorescenceIntensity/ProteinConcentration is a sigmoid curve, but the central 
  # part, which represents our data, is nearly-linear. Therefore, we do not need to 
  # log our data. 
  # df.filt <- log10(df.filt)
  
  # apply min max (DISCARD min-max rescaling)
  #df.filt <- data.frame(Time=10 * 0:(nrow(df.filt)-1),
  #                      apply(df.filt, 2, normalise), 
  #                      check.names=FALSE)
  
  df.filt <- data.frame(Time=10 * 0:(nrow(df.filt)-1),
                        df.filt, 
                        check.names=FALSE)  
  
  return (df.filt)
}





# plot the table as heatmap.Rows are sorted by maximum increasing
# :param filename: the filename containing the table to plot
tc_heatmap <- function(folder, filename, labCol, title="Time courses", df.thres=17) {
  df <- read.table(file.path(folder, paste0(filename)), header=TRUE, 
                   na.strings=c("nan", "NA", ""), 
                   dec = ".", sep=',')
  #print(df)
  
  # rename repeat names
  colnames(df) <- c('time', paste0("x",1:(ncol(df)-1)))
  
  # traspose the data
  df <- t(df)
  
  # the first row becomes the header
  colnames(df) = df[1, ]
  # remove the first row.
  df = df[-1, ]          
  repeats <- nrow(df)
  
  # reduce the data
  if(ncol(df) > df.thres) {
    df <- df[1:df.thres,]
  }
  
  # plot the heatmap
  
  # creates a 5 x 5 inch image
  png(paste(file_path_sans_ext(filename), ".png", sep=""),    # create PNG for the heat map
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  # normalise each row within [0,1]
  df <- t(apply(df, 1, normalise))
  
  # replace NA with 0
  df[is.na(df)] <- 0
  
  par(cex.main=1.5, cex.axis=0.9)
  g <- heatmap.2(as.matrix(df),
                 main = title, # heat map title
                 col=matlab.like(256),
                 scale="none",
                 dendrogram="none", Rowv = FALSE, Colv = FALSE,
                 density.info="none",  # turn off density plot inside color legend
                 key=FALSE, symkey=FALSE, # turn off key                 
                 trace="none",        # turn off trace lines inside the heat map
                 labRow=FALSE,        # turn off row labels
                 margins = c(7, 4),   # reduce the right margin
                 lwid=c(4,25), lhei=c(4,25), # adjust the margin when key=FALSE
                 
                 labCol=labCol,
                 srtCol=45,
                 cexCol=2.5,
                 cexRow=1.6
                 # don't use this now as I haven't found a way to increase this fonts..
                 #xlab="Time (s)", ylab="Repeats by peak time"
  )
  # This works but is not elegant..
  mtext("Time [s]", side=1, line=4, cex=2.8)
  mtext(paste0("normalised repeats (n=",repeats, ")"), side=4, line=1, cex=2.8)
  
  # close the PNG device  
  dev.off()
  return(g)
}





### -------------------- ###



### STATISTICS ###


# Basic scatter plot
scatter_plot <- function(X, Y, se=TRUE, annot.size=5.5, annot.x=6.5, annot.y=1) {
  df <- data.frame(X, Y)
  # linear model for my data
  fit <- lm(Y ~ X, data = df)
  digits <- 2
  g <- ggplot(df, aes(x=X, y=Y, group=X)) +
    geom_point() +
    geom_smooth(method=lm, se=se, aes(group=1)) +
    annotate("text", x=annot.x, y=annot.y, label = equation(fit, digits=digits), parse=TRUE, size=annot.size) +
    annotate("text", x=annot.x, y=annot.y-0.04, label = corr.coef(X, Y, digits=digits), parse=TRUE, size=annot.size) +
    theme_basic(24)
  
  return(g)
}


# Basic violin plot
violin_plot <- function(X, Y) {
  df <- data.frame(X, Y)
  g <- ggplot(df, aes(x=factor(X), y=Y, group=X)) +
    geom_violin(trim = FALSE, size=0.3) + 
    geom_jitter(height = 0, width = 0.1, size=0.5) +
    theme_basic(24)
}



### -------------------- ###




### UPPER AND LOWER VALUES ###


# return the upper values of the peakillations (the peaks)
hv <- function(vec, thres=0.5, ini.tp=0) {
  vec <- vec[!is.na(vec) & as.numeric(names(vec))>=ini.tp & vec>thres]
  return(vec)
}

# return the lower values of the peakillations
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
    labs(x='Time [s]', y='Intensity [a.u.]') +
    theme_basic(base_size=28)
  ggsave(file.path(location.results, paste0(filename, '_upper_values.png')), width=6, height=4, dpi=300)
  write.table(data.plot.hv, file=file.path(location.results, paste0(filename, '_upper_values', suffix)), row.names=FALSE, quote=FALSE, sep=',')
  
  # density plot
  g <- ggplot(data.plot.hv, aes(x=val)) + 
    geom_density(colour = "black", fill = "#56B4E9", alpha=0.5) +
    labs(x='Intensity [a.u.]') +
    theme_basic(base_size=28)
  ggsave(file.path(location.results, paste0(filename, '_upper_values_density.png')), width=6, height=4, dpi=300)
  
  
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
    labs(x='Time [s]', y='Intensity [a.u.]') +
    theme_basic(base_size=28)
  ggsave(file.path(location.results, paste0(filename, '_lower_values.png')), width=6, height=4, dpi=300)
  write.table(data.plot.lv, file=file.path(location.results, paste0(filename, '_lower_values', suffix)), row.names=FALSE, quote=FALSE, sep=',')
  
  # density plot
  g <- ggplot(data.plot.lv, aes(x=val)) + 
    geom_density(colour = "black", fill = "#56B4E9", alpha=0.5) +
    labs(x='Intensity [a.u.]') +
    theme_basic(base_size=28)
  ggsave(file.path(location.results, paste0(filename, '_lower_values_density.png')), width=6, height=4, dpi=300)
  
  
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