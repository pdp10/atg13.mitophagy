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



###########
# LOAD DATA
###########

suffix <- '.csv'
location.data <- file.path('..', '..', 'data')
filename.frame <- 'selected_movies__17'
location.results <- file.path('.')

selected.frames <- read.table( file.path(location.data, paste0(filename.frame, suffix)), header=TRUE, na.strings="NA", dec=".", sep=",")



#####################################
# Extract statistics for MT diameters
#####################################

# Extract mean MT diameters:
# 1) MT diameters were measured for each frame
# 2) the mean of MT diameters was computed for each frame
# 3) the means of MT diameters were plotted
df.diam.means <- data.frame(File=character(), Oscillation=integer(), Mean_Diameter=double(), Repeat=integer())

# Extract all MT diameters (all measurements):
# 1) MT diameters were measured for each frame
# 2) The MT diameters were plotted
df.diam <- data.frame(File=character(), Oscillation=integer(), Diameter=double())

for(i in 1:nrow(selected.frames)) {
  frame <- selected.frames[i,1]
  oscillations <- selected.frames[i,2]
  diameters <- read.table( file.path(location.results, paste0(frame, suffix)), header=TRUE, na.strings="NA", dec=".", sep=",")$Length
  # add means
  df.diam.means <- rbind(df.diam.means, data.frame(File=frame, Oscillation=oscillations, Mean_Diameter=mean(diameters), Repeat=length(diameters)))
  # add raw data
  n = length(diameters)
  df.diam <- rbind(df.diam, data.frame(File=rep(frame,n), Oscillation=rep(oscillations,n), Diameter=diameters))
}

# save table
write.table(df.diam.means, file=file.path(location.data, paste0(filename.frame, suffix)), row.names=FALSE, quote=FALSE, sep=',')
# save raw data table
write.table(df.diam, file=file.path(location.data, paste0(filename.frame, '_raw_data', suffix)), row.names=FALSE, quote=FALSE, sep=',')

# Name shortening
df.diam$File <- gsub("MAX_Cell", "", df.diam$File)
#print(df.diam)


######
# Plot
######

# Plot mean diameters vs oscillations
g <- scatter_plot(df.diam.means$Oscillation, df.diam.means$Mean_Diameter, annot.x=3.0, annot.y=1.0, annot.size=5) + 
  labs(title='', x='oscillation #', y='MT diameter [um]') + 
  scale_x_discrete(limits=c("1","2","3","4","5","6"))
ggsave(file.path(location.results, paste0('mean_mt_diam_vs_peaks_scatterplot','.png')), width=4, height=4, dpi=300)

# Q-Q plot of mean diameters
q <- qq_plot(df.diam.means$Mean_Diameter)
ggsave(file.path(location.results, paste0('mean_mt_diam_qqplot','.png')), width=4, height=4, dpi=300)

# Density plot of mean diameters
q <- density_plot(df.diam.means$Mean_Diameter) + 
  labs(title='', x='MT diameter [um]') + 
  annotate("text", x=0.85, y=4.5, label = paste0('n=', length(df.diam.means$Mean_Diameter)), size=5.5) + 
  annotate("text", x=0.85, y=4, label = paste0('mu=', round(mean(df.diam.means$Mean_Diameter), digits=3)), size=5.5) + 
  annotate("text", x=0.85, y=3.5, label = paste0('SD=', round(sd(df.diam.means$Mean_Diameter), digits=3)), size=5.5) + 
ggsave(file.path(location.results, paste0('mean_mt_diam_density','.png')), width=4, height=4, dpi=300)

df.diam.means.stats <- data.frame(stat=c('n', 'mu', 'sd'), 
                                  value=c(length(df.diam.means$Mean_Diameter), mean(df.diam.means$Mean_Diameter), sd(df.diam.means$Mean_Diameter)))
# save table
write.table(df.diam.means.stats, file=file.path(location.results, paste0('mean_mt_diameter_stats', suffix)), row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')



###############
# Plot raw data
###############

# Plot raw diameter measurements vs oscillation number
g <- box_plot(df.diam$File, df.diam$Diameter) + 
  labs(title='MT diameter samples', x='Frames', y='MT diameter [um]') +
  #theme(axis.text.x=element_text(angle=90,hjust=1, vjust=0.3)) + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(file.path(location.results, paste0('mt_diam_density_per_frame__boxplot','.png')), width=6, height=4, dpi=300)
g <- violin_plot(df.diam$File, df.diam$Diameter) + 
  labs(title='MT diameter samples', x='Frames', y='MT diameter [um]') +
  #theme(axis.text.x=element_text(angle=90,hjust=1, vjust=0.3)) + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(file.path(location.results, paste0('mt_diam_density_per_frame__violinplot','.png')), width=6, height=4, dpi=300)


# diameter densities by frame
q <- density_plot_wcolour(df.diam$Diameter, df.diam$File) + 
  labs(x='MT diameter [um]')
# split the densities in subplots
q <- q + facet_wrap(~variable, ncol=4) + 
  ggtitle('Densities of MT diameters per frame') + 
  theme(legend.position="none")
q <- ggplotly(q)
ggsave(file.path(location.results, paste0('mt_diam_density_per_frame','.png')), width=9, height=9, dpi=300)


# # Plot raw diameter measurements vs oscillation number
# g <- violin_plot(df.diam$Oscillation, df.diam$Diameter) + 
#   labs(title='MT diameters', x='oscillation #', y='MT diameter [um]')
# ggsave(file.path(location.results, paste0('mt_diam_density_per_osc_num__violinplot','.png')), width=4, height=4, dpi=300)
# 
# # diameter densities by oscillation number
# q <- density_plot_wcolour(df.diam$Diameter, df.diam$Oscillation) + 
#   labs(x='MT diameter [um]')
# # split the densities in subplots
# q <- q + facet_wrap(~variable) + 
#   ggtitle('Densities of MT diameters per oscillation #') + 
#   theme(legend.position="none")
# ggplotly()
# ggsave(file.path(location.results, paste0('mt_diam_density_per_osc_num','.png')), width=7, height=4, dpi=300)




#################
# Normality tests
#################

# Shapiro-Wilk normality test 
swt <- shapiro.test(df.diam.means$Mean_Diameter)

# kurtosis test
kt <- kurtosis.test(df.diam.means$Mean_Diameter)

# skewness test 
st <- skew.test(df.diam.means$Mean_Diameter)

df.stats <- data.frame(test=c('shapiro.test', 'kurtosis.test', 'skew.test'), 
                       pvalue=c(swt$p.value, kt, st))
# save table
write.table(df.stats, file=file.path(location.results, paste0('mean_mt_diameter_normality_tests', suffix)), row.names=FALSE, quote=FALSE, sep=',')

