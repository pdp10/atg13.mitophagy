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

# calculate the delays between time points. This is zero between the first two time points
# delays are assumed to increase with time.
calculate.delays <- function(vec) {
  vec.delay <- c(0)
  delay1 <- 0
  if(length(vec) > 2) {
    delay1 <- vec[2]-vec[1]
  }
  for(i in 2:length(vec)) {
    vec.delay <- c(vec.delay, vec[i]-vec[i-1]-delay1)
  }
  return(vec.delay)
}


# merge the values having column name difference (delay) less than time.diff. 
# e.g. X1=560, X2=580 => X1=570
merge.values <- function(vec, time.diff=20) {
  vec.new <- c()
  vec.times <- as.numeric(names(vec))
  i <- 2
  while(i <= length(vec.times)) {
    curr.diff <- vec.times[i]-vec.times[i-1]
    if(curr.diff <= time.diff) {
      # the new element will be the average
      vec.new <- c(vec.new, (vec[i]+vec[i-1])/2)
      names(vec.new)[i-1] <- (vec.times[i]+vec.times[i-1])/2
      # we need to skip one element
      i <- i+1
    } else {
      vec.new <- c(vec.new, vec[i-1])
    }
    i <- i + 1
  }
  # add the last value
  if(vec.times[length(vec.times)]-as.numeric(names(vec.new))[length(vec.new)] > time.diff) {
    vec.new <- c(vec.new, vec[length(vec)])
  }
  return(vec.new)
}


##########
# DATA SET
##########

suffix <- '.csv'
location <- file.path('..', '..', 'data')
filename <- 'mitophagy_summary_intensity_mean_ch2__synchronised_filtered_regularised'

data <- read.table( file.path(location, paste0(filename, suffix)), header=TRUE, na.strings="NA", dec=".", sep=",", row.names=1)


####################
# Extract the delays
####################

# extract the delays from the mean of the time courses
# Although this method collects a little amount of data, these are accurate
# because the mean is smoother and the peaks are well defined. To add more data points, 
# we can calculate the delays from the lowest points of the mean oscillation. 
data.means <- apply(data, 1, mean, na.rm=TRUE)
# highest values time points and intensities of the mean oscillations
## with min-max scaling
#hv <- osc.hv(data.means, thres=0.4, ini.tp=0)
## without min-max scaling
hv <- osc.hv(data.means, thres=1500, ini.tp=0)
# merge values if the time difference is smaller than time.diff 
hv <- merge.values(hv, time.diff=20)
hv.tc <- as.numeric(names(hv))
delay <- calculate.delays(hv.tc)
data.plot.hv <- data.frame(time=hv.tc, val=hv, delay, type=rep('high', length(hv)))

# lowest values time points and intensities of the mean oscillations
## with min-max scaling
#lv <- osc.lv(data.means, thres=0.25, ini.tp=0)
## with min-max scaling
lv <- osc.lv(data.means, thres=1300, ini.tp=0)
# merge values if the time difference is smaller than time.diff 
lv <- merge.values(lv, time.diff=20)
lv.tc <- as.numeric(names(lv))
delay <- calculate.delays(lv.tc)
data.plot.lv <- data.frame(time=lv.tc, val=lv, delay, type=rep('low', length(lv)))



# ---------------------------------

# Save the mean of the highest values. THE ACTUAL HIGHEST VALUES, not from the mean.
# This corresponds to ATG13_max in the model
data.max <- apply(data, 1, max, na.rm=TRUE)
## with min-max scaling
#hv.max <- osc.hv(data.max, thres=0.7)
## without min-max scaling
hv.max <- osc.hv(data.max, thres=1400)
# we remove 0s
#hv.max[hv.max == 0] <- NA
hv.max.mean <- mean(hv.max, na.rm=TRUE)
names(hv.max.mean) <- 'mean_highest_vals'
write.table(hv.max.mean, file=paste0(filename, '_hv_mean', suffix), row.names=TRUE, col.names=FALSE, quote=FALSE, sep=',')

# Save the mean of the lowest values. THE ACTUAL LOWEST VALUES, not from the mean.
# This corresponds to ATG13_min in the model
data.min <- apply(data, 1, min, na.rm=TRUE)
## with min-max scaling
#lv.min <- osc.lv(data.min, thres=0.05)
## without min-max scaling
lv.min <- osc.lv(data.min, thres=1300)
# we remove 0s
#lv.min[lv.min == 0] <- NA
lv.min.mean <- mean(lv.min, na.rm=TRUE)
names(lv.min.mean) <- 'mean_lowest_vals'
write.table(lv.min.mean, file=paste0(filename, '_lv_mean', suffix), row.names=TRUE, col.names=FALSE, quote=FALSE, sep=',')


# ---------------------------------



## LET'S NOW CALCULATE the time from low to peak. We know that ATG13 peak times follow a normal distribution from ATG13 in generic autophagy data.
# Therefore, we compute these time differences, calculate the mean+sd. Then we will sample from the event delay for switching ATG13 kinetic rate constants 
# from N(mean, sd) .

data.plot.sorted.time <- data.plot[with(data.plot, order(time)), ]

# we cut off the first point (which is a TOP and we do not care), and points after 600s because they are too noisy.
data.plot.sorted.time <- data.plot.sorted.time[data.plot.sorted.time$time>40 & data.plot.sorted.time$time<600,]

# extract the high and the low
data.plot.sorted.time.high <- data.plot.sorted.time[data.plot.sorted.time$type=='high',]
data.plot.sorted.time.low <- data.plot.sorted.time[data.plot.sorted.time$type=='low',]
low.high.time.diff <- c()
for(i in 1:nrow(data.plot.sorted.time.high)) {
  low.high.time.diff <- c(low.high.time.diff, data.plot.sorted.time.high$time[i] - data.plot.sorted.time.low$time[i])
}

low.high.time.diff.mean <- mean(low.high.time.diff)
low.high.time.diff.sd <- sd(low.high.time.diff)

df.low.high <- data.frame(mean=low.high.time.diff.mean, sd=low.high.time.diff.sd)
write.table(df.low.high, file=paste0(filename, '__atg13_accum_time_stats', suffix), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=',')



# plot
data.plot <- rbind(data.plot.hv, data.plot.lv)
g <- ggplot() + 
  geom_point(data=data.plot, aes(x=time, y=delay, col=type), size=2.5) +  
  labs(title='Peak delays', x='Time [s]', y='Delay [s]', color='Peak') +
  annotate("text", x=320, y=155, label=paste0("high-low time diff"), parse=FALSE, size=5.5) +  
  annotate("text", x=185, y=130, label=paste0("n=", length(low.high.time.diff)), parse=FALSE, size=5.5) +  
  annotate("text", x=250, y=110, label=paste0("mu=", round(df.low.high$mean, digits=3)), parse=FALSE, size=5.5) +
  annotate("text", x=230, y=90, label=paste0("SD=", round(df.low.high$sd, digits=3)), parse=FALSE, size=5.5) +  
  theme_basic(base_size=24)
ggsave(paste0(filename, '_delay.png'), width=5, height=4, dpi=300)
write.table(data.plot, file=paste0(filename, '_delay', suffix), row.names=FALSE, quote=FALSE, sep=',')
