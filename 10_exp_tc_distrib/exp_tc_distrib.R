
library(tools)
library(gplots)

# palette
library(colorRamps) # matlab.like
library(RColorBrewer)



# plot the table as heatmap.Rows are sorted by maximum increasing
#
# :param filename: the filename containing the table to plot
tc_heatmap <- function(folder, filename, labCol) {
  df <- read.table(paste0(folder, filename), header=TRUE, 
                   na.strings=c("nan", "NA", ""), 
                   dec = ".", sep=',')
  #print(df)
  
  # rename repeat names
  colnames(df) <- c('time', paste0("x",1:(ncol(df)-1)))
  
  # traspose the data
  df <- t(df)
  
  # the first row becomes the header
  colnames(df) = df[1, ]
  # remove the first row.
  df = df[-1, ]          
  
  
  # reduce the data
  if(ncol(df) > 17) {
    df <- df[1:17,]
  }

  
  # plot the heatmap


  # creates a 5 x 5 inch image
  png(paste(file_path_sans_ext(filename), ".png", sep=""),    # create PNG for the heat map
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size

  
  # normalise each row within [0,1]
  #df <- t(apply(df, 1, function(x)(x-min(x))/(max(x)-min(x))))
    
  par(cex.main=1.5, cex.axis=0.9)
  g <- heatmap.2(as.matrix(df),
            main = "Time courses", # heat map title
            density.info="none",  # turns off density plot inside color legend
            key=TRUE, symkey=FALSE,
            trace="none",        # turns off trace lines inside the heat map
            margins =c(5,5),     # widens margins around plot
            col=matlab.like(256),
            scale="none",
            dendrogram="none", Rowv = FALSE, Colv = FALSE,
            labCol=labCol, srtCol=45,
            cexCol=1.6,
            cexRow=1.1
            # don't use this now as I haven't found a way to increase this fonts..
            #xlab="Time (s)", ylab="Repeats by peak time"
            )
  # This works but is not elegant..
  mtext("       Time (s)", side=1, line=4, cex=2 )
  mtext("Repeats by peak time              ", side=4, line=1, cex=2)
  
  # close the PNG device  
  dev.off()
  return(g)
}




# Generate x axis labels
expLabCol <- seq(0, 1000, 10)
expLabCol[expLabCol %% 50 != 0] <- NA


### EXPERIMENTAL DATA ###
# prepare datasets 
# THIS IS NOT NEEDED
#chunk_size <- 45000
#df_col <- 3
#split_df("ds_atg13_ctrl.csv", chunk_size, df_col)
#split_df("ds_atg13_wrtm.csv", chunk_size, df_col)

# plot datasets 
tc_heatmap("../data/", "mitophagy_summary_intensity_mean_ch2__synchronised_filtered_regularised.csv", expLabCol)
