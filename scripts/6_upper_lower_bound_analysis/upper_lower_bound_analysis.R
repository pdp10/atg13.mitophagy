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



# analyse the delay


library(ggplot2)

source('../utilities/plots.R')


# return the upper values of the oscillations (the peaks)
osc.hv <- function(vec, thres=0.5, ini.tp=0) {
  vec <- vec[as.numeric(names(vec))>=ini.tp]
  hv <- c()
  prev <- 0
  prev.prev <- 0
  for(i in 1:length(vec)) {
    # prev must be above a certain threshold
    # prev must be gt than prev.prev and vec[i]
    if(!is.na(prev.prev) & !is.na(prev) & !is.na(vec[i])) {
      if(prev > thres & (vec[i] < prev & prev > prev.prev)) {
        hv <- c(hv, prev)
      }
    }
    prev.prev <- prev
    prev <- vec[i]
  }
  return(hv)
}

# return the lower values of the oscillations
osc.lv <- function(vec, thres=0.5, ini.tp=0) {
  vec <- vec[as.numeric(names(vec))>=ini.tp]
  lv <- c()
  prev <- 0
  prev.prev <- 0
  for(i in 1:length(vec)) {
    # prev must be below a certain threshold
    # prev must be lt than prev.prev and vec[i]
    if(!is.na(prev.prev) & !is.na(prev) & !is.na(vec[i])) {
      if(prev < thres & vec[i] > prev & prev < prev.prev) {
        lv <- c(lv, prev)
      }
    }
    prev.prev <- prev
    prev <- vec[i]
  }
  return(lv)
}



##########
# DATA SET
##########

suffix <- '.csv'
location <- file.path('..', '..', 'data')
filename <- 'mitophagy_summary_intensity_mean_ch2__synchronised_filtered_regularised'

data <- read.table( file.path(location, paste0(filename, suffix)), header=TRUE, na.strings="NA", dec=".", sep=",", row.names=1)



# ##########################
# # PLOT SPLINE TIME COURSES
# ##########################
# 
# # we also plot the splines of these time courses
# # spar: smoothing parameter, typically (but not necessarily) in (0,1]
# # spars <- c(0.3, 0.4, 0.5, 0.6)
# # 0.4 is the best trade-off
# spar <- 0.4
# 
# data.spline <- spline.data.frame(data, spar)
# write.table(data.spline, file=paste0(location, file, '_spline', suffix), sep=",", row.names=FALSE)
#   
# # Plot all time course splines. x axis is expanded. Plots have different timings
# plots <- plot_combined_tc(data.spline, expand.xaxis=TRUE)
# plot.arrange <- do.call(grid.arrange, c(plots$plots, ncol=5, bottom=paste0(file, '_spline', suffix)))
# ggsave(paste0(filename, '_spline_all.png'), plot=plot.arrange, width=13, height=13, dpi=300)
# 
# data <- data.spline


####################
# Extract the delays
####################

# we manually extracted thresholds for each time course. These are used to determine the highest and lowest values of 
# the oscillations.
print(colnames(data))
#thres.hv <- c(1500,1500,1500,1500,1500,1400,1500,1500,2500,1700,1500,1600,1600,1700,1400,1300,1200)
#thres.lv <- c(1300,1300,1300,1400,1400,1100,1200,1200,1500,1500,1200,1500,1500,1400,1200,1200,1100)

thres.hv <- c(1500,1500,1500,1500,1500,1400,1500,1500,2500,1700,1500,1600,1600,1700,1400,1300,1200)
thres.lv <- c(1300,1300,1300,1300,1300,1100,1200,1200,1300,1300,1200,1300,1300,1300,1200,1200,1100)

# EXTRACT THE UPPER VALUES (the peaks)

# extract the delays from the time courses
data.plot.hv <- data.frame(time=numeric(0), frame=character(0), val=numeric(0))
for(i in 1:ncol(data)) {
 col.i <- data[,i]
 names(col.i) <- row.names(data)
 ## with min-max scaling
 #hv.i <- osc.hv(col.i, thres=0.50)
 ## without min-max scaling
 hv.i <- osc.hv(col.i, thres=thres.hv[i])
 hv.tc.i <- as.numeric(names(hv.i))
 data.plot.hv <- rbind(data.plot.hv, data.frame(time=hv.tc.i, frame=colnames(data)[i], val=hv.i))
}
data.plot.hv <- cbind(data.plot.hv, pos=rep('top', nrow(data.plot.hv)))

# plot
g <- ggplot() + 
  geom_point(data=data.plot.hv, aes(x=time, y=val)) +  
  labs(title='Oscillation upper bounds', x='Time [s]', y='Intensity [a.u.]') +
  theme_basic()
ggsave(paste0(filename, '_upper_bounds.png'), width=4, height=3, dpi=300)
write.table(data.plot.hv, file=paste0(filename, '_upper_bounds', suffix), row.names=FALSE, quote=FALSE, sep=',')

# density plot
g <- ggplot(data.plot.hv, aes(x=val)) + 
  geom_density(colour = "black", fill = "#56B4E9", alpha=0.5) +
  labs(title='Density of upper bounds', x='signal intensity [a.u.]') +
  theme_basic()
ggsave(paste0(filename, '_upper_bounds_density.png'), width=4, height=3, dpi=300)


# EXTRACT THE LOWER VALUES (the lower peaks)

# extract the delays from the time courses
data.plot.lv <- data.frame(time=numeric(0), val=numeric(0))
for(i in 1:ncol(data)) {
  col.i <- data[,i]
  names(col.i) <- row.names(data)
  ## with min-max scaling
  #lv.i <- osc.lv(col.i, thres=0.30)
  ## without min-max scaling
  lv.i <- osc.lv(col.i, thres=thres.lv[i])
  lv.tc.i <- as.numeric(names(lv.i))
  data.plot.lv <- rbind(data.plot.lv, data.frame(time=lv.tc.i, val=lv.i))
}
data.plot.lv <- cbind(data.plot.lv, pos=rep('bottom', nrow(data.plot.lv)))

# plot
g <- ggplot() + 
  geom_point(data=data.plot.lv, aes(x=time, y=val)) +  
  labs(title='Oscillation lower bounds', x='Time [s]', y='Intensity [a.u.]') +
  theme_basic()
ggsave(paste0(filename, '_lower_bounds.png'), width=4, height=3, dpi=300)
write.table(data.plot.lv, file=paste0(filename, '_lower_bounds', suffix), row.names=FALSE, quote=FALSE, sep=',')

# density plot
g <- ggplot(data.plot.lv, aes(x=val)) + 
  geom_density(colour = "black", fill = "#56B4E9", alpha=0.5) +
  labs(title='Density of lower bounds', x='signal intensity [a.u.]') +
  theme_basic()
ggsave(paste0(filename, '_lower_bounds_density.png'), width=4, height=3, dpi=300)


# write the peaks stats

peaks.stats <- data.frame(type=c('top', 'bottom'), mean=c(mean(data.plot.hv$val), mean(data.plot.lv$val)), 
                          sd=c(sd(data.plot.hv$val), sd(data.plot.lv$val)))
write.table(peaks.stats, file=paste0(filename, '_peaks_stats', suffix), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=',')
