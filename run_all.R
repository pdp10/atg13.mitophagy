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



#4
print('4_signal_decrease_data')
setwd('4_signal_decrease_data')
source('signal_decrease_data.R')
rm(list = ls())
setwd('../')

#5
print('5_mitophagy_time_courses_plots')
setwd('5_mitophagy_time_courses_plots')
source('mitophagy_time_courses_plots.R')
rm(list = ls())
setwd('../')

#6
print('6_synchronised_time_courses')
setwd('6_synchronised_time_courses')
source('synchronised_time_courses.R')
rm(list = ls())
setwd('../')

#7
print('7_time_courses_w_adjust_green_ch')
setwd('7_time_courses_w_adjust_green_ch')
source('time_courses_w_adjust_green_ch.R')
rm(list = ls())
setwd('../')

