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


######################
# TIME COURSE DATA SET
######################

suffix <- '.csv'

location.tc <- '../data/'
filenames.tc <- c('mitophagy_summary_intensity_mean_ch2__synchronised_filtered_regularised', 
               'mitophagy_summary_intensity_mean_ch2__synchronised_regularised')
filenames.tc.out <- c('mitophagy_atg13_tc_filtered__copasi', 
                   'mitophagy_atg13_tc__copasi')

location.delay <- '../8_delay_analysis/'
filename.delay <- 'mitophagy_summary_intensity_mean_ch2__synchronised_filtered_regularised_delay'
filename.delay.out <- 'mitophagy_atg13_tc_filtered_delay__copasi'

location.copasi <- '../../Models/'



######################
# Time course data set
######################

for(i in 1:length(filenames.tc)) {
  # DATA LOADING
  data <- read.table( paste0(location.tc, filenames.tc[i], suffix), header=TRUE, na.strings="NA", dec=".", sep=",")
  
  # MELT
  data.copasi <- melt(data, id.vars=c('Time'), value.name = 'ATG13_obs')
  data.copasi <- subset(data.copasi, select=-c(variable))
  
  # DATA WRITING
  write.table(data.copasi, file=paste0(location.copasi, filenames.tc.out[i], suffix), sep="\t", na="", dec=".", row.names=FALSE, quote=FALSE)
  print(paste0('Copasi data set saved in: ', location.copasi, filenames.tc.out[i], suffix))
  
}


################
# delay data set
################

data.delay <- read.table( paste0(location.delay, filename.delay, suffix), header=TRUE, na.strings="NA", dec=".", sep=",")

data.copasi.delay <- data.delay[,c(1,2,3)]
colnames(data.copasi.delay) <- c('Time', 'ATG13_mean_obs', 'delay_obs')

# set negative values to 0
#data.copasi.delay <- apply(data.copasi.delay, c(1,2), function(x) {if(x<0) x=0 else x})

# DATA WRITING
write.table(data.copasi.delay, file=paste0(location.copasi, filename.delay.out, suffix), sep="\t", na="", dec=".", row.names=FALSE, quote=FALSE)
print(paste0('Copasi data set saved in: ', location.copasi, filename.delay.out, suffix))
