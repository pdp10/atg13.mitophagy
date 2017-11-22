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





## LET'S NOW CALCULATE the time from bottom to peak. We know that ATG13 peak times follow a normal distribution from ATG13 in generic autophagy data.
# Therefore, we compute these time differences, calculate the mean+sd. Then we will sample from the event delay for switching ATG13 kinetic rate constants 
# from N(mean, sd) .

data.plot.sorted.time <- data.plot[with(data.plot, order(time)), ]

# we cut off the first point (which is a TOP and we do not care), and points after 600s because they are too noisy.
data.plot.sorted.time <- data.plot.sorted.time[data.plot.sorted.time$time>40 & data.plot.sorted.time$time<600,]

# extract the top and the bottom
data.plot.sorted.time.top <- data.plot.sorted.time[data.plot.sorted.time$pos=='top',]
data.plot.sorted.time.bottom <- data.plot.sorted.time[data.plot.sorted.time$pos=='bottom',]
bottom.top.time.diff <- c()
for(i in 1:nrow(data.plot.sorted.time.top)) {
  bottom.top.time.diff <- c(bottom.top.time.diff, data.plot.sorted.time.top$time[i] - data.plot.sorted.time.bottom$time[i])
}

bottom.top.time.diff.mean <- mean(bottom.top.time.diff)
bottom.top.time.diff.sd <- sd(bottom.top.time.diff)

df.bottom.top <- data.frame(mean=bottom.top.time.diff.mean, sd=bottom.top.time.diff.sd)
write.table(df.bottom.top, file=paste0(filename, '__atg13_accum_time_stats', suffix), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=',')

