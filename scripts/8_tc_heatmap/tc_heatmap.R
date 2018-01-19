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


# plot the table as heatmap.Rows are sorted by maximum increasing



source('../utilities/plots.R')
source('../utilities/statistics.R')


###########
# LOAD DATA
###########

location <- file.path('..', '..', 'data')
filename <- "mitophagy_summary_intensity_mean_ch2__synchronised_filtered_regularised"
suffix <- '.csv'



###############
# Plot Datasets
###############

# Generate x axis labels
expLabCol <- seq(0, 1000, 10)
expLabCol[expLabCol %% 100 != 0] <- NA

# Plot
tc_heatmap(location, paste0(filename, suffix), expLabCol, title="Exp ATG13", df.thres=17)

