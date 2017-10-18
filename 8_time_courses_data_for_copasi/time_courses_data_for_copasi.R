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



# generate time course data for Copasi


library(reshape2)


##########
# DATA SET
##########

suffix <- '.csv'
location <- '../data/'
filenames <- c('mitophagy_summary_intensity_mean_ch2__synchronised_filtered_regularised', 
               'mitophagy_summary_intensity_mean_ch2__synchronised_regularised')

filenames.out <- c('mitophagy_atg13_tc_filtered__copasi', 
                   'mitophagy_atg13_tc__copasi')

location.copasi <- '../../Models/'



for(i in 1:length(filenames)) {
  # DATA LOADING
  data <- read.table( paste0(location, filenames[i], suffix), header=TRUE, na.strings="NA", dec=".", sep=",")
  
  # MELT
  data.copasi <- melt(data, id.vars=c('Time'), value.name = 'ATG13_obs')
  data.copasi <- subset(data.copasi, select=-c(variable))
  
  # DATA WRITING
  write.table(data.copasi, file=paste0(location.copasi, filenames.out[i], suffix), sep="\t", na="", dec=".", row.names=FALSE, quote=FALSE)
  print(paste0('Copasi data set saved in: ', location.copasi, filenames.out[i], suffix))
  
}

