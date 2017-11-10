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



##########
# DATA SET
##########

suffix <- '.csv'
location <- '../data/'
filename <- 'mitophagy_summary_intensity_mean_ch2__synchronised_filtered_regularised'

data <- read.table( paste0(location, filename, suffix), header=TRUE, na.strings="NA", dec=".", sep=",", row.names=1)


####################
# Extract the delays
####################


# EXTRACT THE UPPER VALUES (the peaks)

# extract the delays from the time courses
# NOTE: this is not accurate! 
# - the single oscillations are not smooth, so setting the threshold (thres) is tricky.
# - smoothing the dataset is possible, but it is tricky to define a good smooth parameter which 
# consistently works for every time course
#data.plot.hv <- data.frame(time=numeric(0), val=numeric(0), delay=numeric(0))
#for(i in 1:ncol(data)) {
#  col.i <- data[,i]
#  names(col.i) <- row.names(data)
#  hv.i <- osc.hv(col.i, thres=0.6)
#  hv.tc.i <- as.numeric(names(hv.i))
#  delay.i <- calculate.delays(hv.tc.i)
#  data.plot.hv <- rbind(data.plot.hv, data.frame(time=hv.tc.i, val=hv.i, delay=delay.i))
#}


# extract the delays from the mean of the time courses
# Although this method collects a little amount of data, these are accurate
# because the mean is smoother and the peaks are well defined. To add more data points, 
# we can calculate the delays from the lowest points of the mean oscillation. 
data.means <- apply(data, 1, mean, na.rm=TRUE)
# highest values time points and intensities of the mean oscillations
hv <- osc.hv(data.means, thres=0.4, ini.tp=0)
hv.tc <- as.numeric(names(hv))
delay <- calculate.delays(hv.tc)
data.plot.hv <- data.frame(time=hv.tc, val=hv, delay, pos=rep('top', length(hv)))

# lowest values time points and intensities of the mean oscillations
lv <- osc.lv(data.means, thres=0.25, ini.tp=0)
lv.tc <- as.numeric(names(lv))
delay <- calculate.delays(lv.tc)
data.plot.lv <- data.frame(time=lv.tc, val=lv, delay, pos=rep('bottom', length(lv)))

# plot
data.plot <- rbind(data.plot.hv, data.plot.lv)
g <- ggplot() + 
  geom_point(data=data.plot, aes(x=time, y=delay, col=pos)) +  
  labs(title='Delays of Mean Oscillation', x='Time [s]', y='Delay [s]') +
  theme_basic()
ggsave(paste0(filename, '_delay.png'), width=4, height=3, dpi=300)
write.table(data.plot, file=paste0(filename, '_delay', suffix), row.names=FALSE, quote=FALSE, sep=',')



# Save the mean of the lowest values. THE ACTUAL LOWEST VALUES, not from the mean.
# This corresponds to ATG13LT in the model
data.min <- apply(data, 1, min, na.rm=TRUE)
lv.min <- osc.lv(data.min, thres=0.05)
# we remove 0s
#lv.min[lv.min == 0] <- NA
lv.min.mean <- mean(lv.min, na.rm=TRUE)
names(lv.min.mean) <- 'mean_lowest_vals'
write.table(lv.min.mean, file=paste0(filename, '_lv_mean', suffix), row.names=TRUE, col.names=FALSE, quote=FALSE, sep=',')

