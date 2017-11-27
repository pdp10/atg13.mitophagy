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


library(ggplot2)
library(grid)
library(gridExtra)

source('../utilities/plots.R')





##########
# DATA SET
##########

suffix <- '.csv'
location <- '../data/'
files <- c('mitophagy_summary_intensity_mean_ch2__synchronised', 
           'mitophagy_summary_intensity_mean_ch2__synchronised_filtered') 

location.regression <- '../6_synchronised_time_courses/'
files.regression <- c('mitophagy_summary_intensity_mean_ch2__synchronised_linear_regression_data', 
                      'mitophagy_summary_intensity_mean_ch2__synchronised_filtered_linear_regression_data')



for(i in 1:length(files)) {
  
  ##############
  # DATA LOADING
  ##############
  data <- read.table( paste0(location, files[i], suffix), header=TRUE, na.strings="NA", dec=".", sep=",")
  data.regr <- read.table( paste0(location.regression, files.regression[i], suffix), header=TRUE, na.strings="NA", dec=".", sep=",") 
  
  print(files[i])
  print(data.regr)

    
  #########################################################################
  # REGULARISE TIME COURSES USING REGRESSION LINE COMPUTED FROM THE AVERAGE
  #########################################################################
  print('plotting combined tc')
  for(j in 2:ncol(data)) {
    data[,j] <- data[,j] - data.regr$slope * data[,1]
    #plot(x=data[,1], y=data[,j], type='l')
  }
  
  # re-apply min max 
  #data[, 2:ncol(data)] <- data.frame(apply(data[, 2:ncol(data)], 2, normalise), check.names=FALSE)
  
  write.table(data, file=paste0(location, files[i], '_regularised', suffix), row.names=FALSE, quote=FALSE, sep=',')
  
  
  ###################
  # PLOT TIME COURSES
  ###################
  # Plot all time courses. x axis is expanded. Plots have different timings
  # NOTE: we remove the synchronisation for this plots as it is not necessary
  filename <- paste0(gsub("_synchronised", "", files[i]), "_regularised")
  
  
  data.unsynch <- data.frame(Time=data[,1], apply(data[,2:ncol(data)], 2, function(x){c(x[which.min(is.na(x)):length(x)], rep(NA, which.min(is.na(x))-1))}))
  write.table(data.unsynch, file=paste0(location, filename, suffix), row.names=FALSE, quote=FALSE, sep=',')
  
  plots <- plot_combined_tc(data.unsynch, expand.xaxis=TRUE)
  plot.arrange <- do.call(grid.arrange, c(plots$plots, ncol=5, bottom=paste0(filename, suffix)))
  ggsave(paste0(filename, '_all.png'), plot=plot.arrange, width=13, height=13, dpi=300)
  
  
  
  ##################################
  # PLOT TIME COURSES (SYNCHRONISED)
  ##################################
  
  # rename the columns
  colnames(data) <- c("Time", gsub("MAX_Cell", "", tail(colnames(data), ncol(data)-1)))
  # plot the synchronised filtered time courses
  plot.arrange <- plot_synchronised_tc(data, paste0(files[i], '_regularised'), ylab='Normalised Intensity Mean [a.u.]')
  colnames(data)[2:ncol(data)] <- paste0("MAX_Cell", colnames(data)[2:ncol(data)])
}














