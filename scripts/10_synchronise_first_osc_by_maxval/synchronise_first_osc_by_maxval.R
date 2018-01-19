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


# synchronise time courses by maximum value


# retrieve this script path 
args <- commandArgs(trailingOnly=FALSE)
SCRIPT_PATH <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
if(length(SCRIPT_PATH) > 0) {
  # we use this when Rscript is used from a different directory
  SCRIPT_PATH <- normalizePath(SCRIPT_PATH)
  source(file.path(SCRIPT_PATH, '../utilities/plots.R'))
  source(file.path(SCRIPT_PATH, '../utilities/statistics.R'))
} else {
  # we use this when Rscript is used from this directory
  source('../utilities/plots.R')
  source('../utilities/statistics.R')
}



###########
# LOAD DATA 
###########

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) {
  location.data <- file.path('..', '..', 'data')
  location.results <- '.'
  filenames <- c('mitophagy_summary_intensity_mean_ch2__synchronised_filtered_regularised_first_oscillation')
  yaxis.label <- c('IntensityMean')
  # Synchronise tc  
  for(i in 1:length(filenames)) {
    sync_tc_main(location.data, location.results, filenames[i], yaxis.label[i])
  }
} else {
  location.data <- dirname(normalizePath(args[1]))
  location.results <- SCRIPT_PATH  
  filename <- sub("^([^.]*).*", "\\1", basename(args[1]))
  yaxis.label <- args[2]
  sync_tc_main(location.data, location.results, filename, yaxis.label)
}


