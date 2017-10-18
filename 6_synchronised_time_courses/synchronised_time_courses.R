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



# plot sinchronised original and splined time course data


library(ggplot2)
library(grid)
library(gridExtra)


source('../utilities/plots.R')






#########
# MY DATA
#########
location <- '../data/'
files <- c('mitophagy_summary_intensity_mean_ch2__synchronised', 'mitophagy_summary_intensity_mean_ch2_spline__synchronised')
suffix <- '.csv'

remove.cols <- c('3_a', '3_b', '5_a', '6_d', '6_f_2', '16s_a')
remove.row.head <- 54
remove.row.tail <- 164



#################
# PROCESS MY DATA
#################

for(f in files) {
  print(f)
  df <- read.table( paste0(location, f, suffix), header=TRUE, na.strings="NA", dec=".", sep=",")
  # rename the columns
  colnames(df) <- c("Time", gsub("MAX_Cell", "", tail(colnames(df), ncol(df)-1)))

  # plot the synchronised time courses
  plot.arrange <- plot_synchronised_tc(df, f, ylab='Intensity Mean [a.u.]')
  
  
  # data filtering
  f <- paste0(f, "_filtered")
  df.filt <- data_filtering(df, remove.cols, remove.row.head, remove.row.tail)

  # plot the synchronised filtered time courses
  plot.arrange <- plot_synchronised_tc(df.filt, f, ylab='Normalised Intensity Mean [a.u.]')

  colnames(df.filt)[2:ncol(df.filt)] <- paste0("MAX_Cell", colnames(df.filt)[2:ncol(df.filt)])
  write.csv(df.filt, file=paste0(location, f, suffix), quote=FALSE, row.names=FALSE)
}


