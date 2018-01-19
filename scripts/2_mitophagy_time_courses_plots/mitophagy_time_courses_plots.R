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



# plot raw time course data



source('../utilities/plots.R')
source('../utilities/statistics.R')




##########
# DATA SET
##########

location <- file.path('..', '..', 'data')
file <- 'mitophagy_summary_intensity_mean_ch2'
suffix <- '.csv'

# Load my data
data <- read.csv(file.path(location, paste0(file, suffix)))



###################
# PLOT TIME COURSES
###################

# Plot all time courses. x axis is expanded. Plots have different timings
plots <- plot_combined_tc(data, expand.xaxis=TRUE)
plot.arrange <- do.call(grid.arrange, c(plots$plots, ncol=5, bottom=paste0(file, suffix)))
ggsave(paste0(file, '_all.png'), plot=plot.arrange, width=13, height=13, dpi=300)



##########################
# PLOT SPLINE TIME COURSES
##########################

# we also plot the splines of these time courses
# spar: smoothing parameter, typically (but not necessarily) in (0,1]
# spars <- c(0.3, 0.4, 0.5, 0.6)
# 0.4 is the best trade-off
spars <- c(0.4)

for(spar in spars) {
  data.spline <- spline.data.frame(data, spar)
  write.table(data.spline, file=file.path(location, paste0(file, '_spline', suffix)), sep=",", row.names=FALSE)
  
  # Plot all time course splines. x axis is expanded. Plots have different timings
  plots <- plot_combined_tc(data.spline, expand.xaxis=TRUE)
  plot.arrange <- do.call(grid.arrange, c(plots$plots, ncol=5, bottom=paste0(file, '_spline', suffix)))
  ggsave(paste0(file, '_spline_all.png'), plot=plot.arrange, width=13, height=13, dpi=300)
}

