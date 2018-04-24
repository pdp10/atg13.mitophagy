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


# Read raw data for Green Channel signal reduction, compute a spline, and export the data for this spline. 


library(ggplot2)

source('../utilities/plots.R')
source('../utilities/statistics.R')



plotting.spline <- function(files, spline.order=20) {
  PLOTS <- list()
  for(f in files) {
    print(f)
    data <- read.csv(file.path(location, f))
    f.noext <- strsplit(f, "\\.")[[1]][1]
    p <- plot_tc_with_spline( data, file=paste0(f.noext, '_spline_data', '.csv'), spline.order, title=gsub('_signal_reduction', '', f.noext), xlab='Time [s]', ylab='Green (Ch2) Sign. Int. [a.u.]')
    PLOTS <- c(PLOTS, list(p))
  }  
  return( list(plots=PLOTS))
}


plotting.regr <- function(files, file='linear_regression_data.csv') {
  PLOTS <- list()
  for(f in files) {
    print(f)
    data <- read.csv(file.path(location, f))
    f.noext <- strsplit(f, "\\.")[[1]][1]
    p <- plot_tc_with_regr_line( data, file=file, title=gsub('_signal_reduction', '', f.noext), xlab='Time [s]', ylab='Green (Ch2) Sign. Int. [a.u.]')
    PLOTS <- c(PLOTS, list(p))
  }  
  return( list(plots=PLOTS))
}



###########
# LOAD DATA
###########

location <- file.path('..', '..', 'data')
files <- list.files(path = location, pattern = "^[MAX_Cell]")
suffix <- '.csv'
filename <- 'mitophagy_summary_intensity_mean_ch2'
lr.filename <- 'linear_regression_data'

if(file.exists(paste0(lr.filename, suffix))){
  print(paste0('refreshing ', lr.filename, suffix))
  file.remove(paste0(lr.filename, suffix))
}


##########
# PLOTTING
##########

# Plot and arrange everything (spline)
print('plot signal reduction with spline curve')
spline.order <- 20
plots <- plotting.spline(files, spline.order)
plot.arrange <- do.call(grid.arrange, c(plots$plots, ncol=3))
ggsave(paste0('signal_reduction_spline_all.png'), plot=plot.arrange, width=8, height=8, dpi=300)


# Plot and arrange everything (regression)
print('plot signal reduction with regression line')
plots <- plotting.regr(files)
plot.arrange <- do.call(grid.arrange, c(plots$plots, ncol=3))
ggsave(paste0('signal_reduction_regr_all.png'), plot=plot.arrange, width=8, height=8, dpi=300)





##### *** THIS DOES NOT WORK WELL. The BACKGROUND WAS ALREADY SUBTRACTED. 
##### WE WILL SUBTRACT THE MEAN REGRESSION LINE INSTEAD.

############################################################################
# REGULARISE TIME COURSES USING REGRESSION LINES for each experimental image
############################################################################

# df.regr <- read.table( paste0(lr.filename, suffix), header=TRUE, na.strings="NA", dec=".", sep=",") 
# df <- read.table( paste0(location, filename, suffix), header=TRUE, na.strings="NA", dec=".", sep=",")
# 
# for(i in 2:ncol(df)) {
#   # retrieve the regression parameters for each column in df
#   regr.param <- df.regr[grepl( gsub("(.*)_.*", "\\1", colnames(df)[i]), df.regr$file ), ]
#   if(nrow(regr.param) == 0) {
#     regr.param <- df.regr[grepl( gsub("(.*)_.*_.*", "\\1", colnames(df)[i]), df.regr$file ), ]
#   }
#   #print(regr.param)
#   df[,i] <- df[,i] - regr.param$slope * df[,1]
#   plot(x=df[,1], y=df[,i], type='l', main=regr.param$file)
# }

###write.table(df, file=paste0(location, filename, '_regularised', suffix), row.names=FALSE, quote=FALSE, sep=',')

