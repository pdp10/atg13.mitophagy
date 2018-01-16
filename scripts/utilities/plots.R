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





library(reshape2)
library(ggplot2)
library(plotly)
library(grid)
library(gridExtra)

library(tools)
library(gplots)
# palette
library(colorRamps) # matlab.like
library(RColorBrewer)


# A simple theme without grid
theme_basic <- function(base_size = 12){
  theme_bw(base_size) %+replace%
    theme(
      panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank()
    )
}


# min max normalisation
# call: apply(df, 2, normalise) 
# to normalise each column within [0,1]
normalise <- function(x, na.rm = TRUE) {
  ranx <- range(x, na.rm = na.rm)
  (x - ranx[1]) / diff(ranx)
}


# return the linear model equation
equation <- function(x, digits=4, show.r2=FALSE) {
  lm_coef <- list(a = round(coef(x)[1], digits=digits),
                  b = round(coef(x)[2], digits=digits),
                  r2 = round(summary(x)$r.squared, digits=digits));
  if(show.r2) {
    if(lm_coef$b >= 0) {
      lm_eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)    
    } else {
      lm_coef$b <- abs(lm_coef$b)
      lm_eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)
    }
  } else {
    if(lm_coef$b >= 0) {
      lm_eq <- substitute(italic(y) == a + b %.% italic(x),lm_coef)    
    } else {
      lm_coef$b <- abs(lm_coef$b)
      lm_eq <- substitute(italic(y) == a - b %.% italic(x),lm_coef)
    }
  }
  return(as.character(as.expression(lm_eq)))
}

# return the linear model equation
corr.coef <- function(x, y, digits=4, method='pearson') {
  my.corr <- cor.test(x, y, method=method)
  coef <- list(estimate = round(as.double(my.corr$estimate[1]), digits=digits),
               p.value = round(as.double(my.corr$p.value), digits=7));
  coef.str <- substitute(italic(r) == estimate*","~~italic(p)*"-"*val == p.value,coef)    
  return(as.character(as.expression(coef.str)));
}


# plot the mean with error bars
compute_mean_error <- function(df, filename, ylab, show.linear.model=FALSE) {
  
  dfnt <- df[,-1]
  time_vec <- subset(df, select=c('Time'))
  
  means <- apply(dfnt, 1, mean, na.rm=TRUE)
  stdevs <- sqrt(apply(dfnt, 1, var,na.rm=TRUE))
  lengths <- apply(dfnt, 1, function(x) length(which(!is.na(x))))  
  stderrs <- stdevs / sqrt(lengths)
  ci95 <- qnorm(.975)*(stderrs)
  
  statdf <- data.frame(Time=time_vec, means=means, stdevs=stdevs, ci95=ci95)
  colnames(statdf)[1] <- 'Time'
  #print(statdf)
  
  if(show.linear.model) {
    fit <- lm(means ~ Time, data = statdf)
  }
  
  # plot mean+sd+ci95
  g <- ggplot(statdf, aes(x=Time, y=means)) + 
    # SD
    geom_errorbar(aes(ymin=means-stdevs, ymax=means+stdevs)) + 
    geom_line(aes(x=Time, y=means), color="black", size=0.5) +
    # CI 95%
    geom_errorbar(aes(ymin=means-ci95, ymax=means+ci95), colour="magenta") +
    geom_line(aes(x=Time, y=means), color="black", size=0.5) +
    labs(x="Time [s]", y=ylab
         , title=paste0('')
         #, title=paste0('n=', ncol(df)-1)
         ) +
    theme_basic(base_size=34)
  
  if(show.linear.model) {
    x0 <- 550 #500
    y0 <- 500 #-0.1
    g <- g +
      stat_smooth(method = "lm", se=TRUE, color="blue", aes(group=1)) + 
      annotate("text", x=x0, y=y0, label = equation(fit), size=6, parse = TRUE)
    
    # export regression data
    data.regr <- data.frame(regression='meansVStime', slope=coef(fit)[2], intercept=coef(fit)[1], check.names = FALSE)
    write.table(data.regr, file=file.path(paste0(filename, '_linear_regression_data.csv')), row.names=FALSE, quote=FALSE, sep=',')
  }
  
  colnames(statdf) <- c("Time", "mean", "sd", "ci95")
  write.table(statdf, file=file.path(paste0(filename, "_stats.csv")), sep=",", row.names=FALSE)
  return (g)
}


# Apply a smooth.spline to a data frame
spline.data.frame <- function(data, spar) {
  # create the spline for time
  data.spline <- data.frame(Time=smooth.spline(data[,1], spar=spar)$y)
  # create the splines for each event
  for(i in 2:ncol(data)) {
    sp <- smooth.spline(na.omit(data[,i]), spar=spar)$y
    # extract $y and add NAs
    sp <- c(sp, rep(NA, nrow(data.spline)-length(sp)))
    data.spline <- cbind(data.spline, sp)
  }
  colnames(data.spline) <- colnames(data)
  return(data.spline)
} 


# plot the time course repeats
plot_tc_repeats <- function(df, ylab) {
  df.melt <- melt(df,id.vars=c("Time"), variable.name = 'repeats')
  g <- ggplot() + geom_line(data=df.melt,aes(x=Time,y=value,color=repeats), size=0.5) +
    labs(x="Time [s]", y=ylab, title=paste0('n=', ncol(df)-1)) +
    theme_basic(base_size=34)
  return(g)
}


# Plot the time course
plot_tc <- function(df, title='title', xlab='Time [s]', ylab='Sign. Int. [a.u.]') {
  g <- ggplot(data=df, aes(x=x, y=y)) + 
    geom_line() + geom_point() +
    labs(title=title, x=xlab, y=ylab) + 
    theme_basic()
  return(g)
}


plot_tc_with_spline <- function(data, file, spline.order, title, xlab='Time [s]', ylab='Green (Ch2) Sign. Int. [a.u.]') {
  data.spline <- spline(data$X, data$Y, n=spline.order)
  data.spline <- data.frame(X=data.spline$x, Y=data.spline$y)
  
  g <- ggplot() + 
    geom_point(data=data, aes(x=X, y=Y), shape=1) +
    geom_line(data=data.spline, aes(x=X, y=Y), color="red", size=1.0) +
    theme_basic() + 
    labs(title=title, x=xlab, y=ylab)
  
  data.spline$X <- data.spline$X
  write.csv(data.spline, file=file, row.names=FALSE)
  
  return(g)
}


plot_tc_with_regr_line <- function(data, file, title, xlab='Time [s]', ylab='Green (Ch2) Sign. Int. [a.u.]') {
  # linear model for my data
  fit <- lm(Y ~ X, data = data)
  digits <- 2
  
  g <- ggplot(data=data, aes(x=X, y=Y)) + 
    geom_point(shape=1) +
    stat_smooth(method = "lm", se=TRUE, color="red") + 
    annotate("text", x=75, y=640, label = equation(fit, digits=digits), parse = TRUE, size=2.6) +
    theme_basic() + 
    labs(title=title, x=xlab, y=ylab)

  # export regression data
  data.regr <- data.frame(file=title, slope=coef(fit)[2], intercept=coef(fit)[1], check.names = FALSE)
  if(!file.exists(file)){
    write.table(data.regr, file=file, row.names=FALSE, quote=FALSE, sep=',')
  } else {
    write.table(data.regr, file=file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',') 
  }
  return(g)
}



plot_combined_tc <- function(data, expand.xaxis=TRUE) {
  samples <- colnames(data)
  PLOTS <- list() 
  for(i in 2:ncol(data)) {
    if(expand.xaxis) {
      non.na <- sum(!is.na(data[,i]))
      #print(non.na)
      df <- data.frame(x=data[1:non.na,1], y=data[1:non.na,i])
    } else {
      df <- data.frame(x=data[,1], y=data[,i])      
    }
    p <- plot_tc(df, title=samples[i], xlab='Time [s]', ylab='Intensity Mean [a.u.]') +
      theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "in"))
    PLOTS <- c(PLOTS, list(p))
  }
  return( list(plots=PLOTS) )
}



# plot the time courses
plot_synchronised_tc <- function(df, filename, ylab) {
  # plot time courses
  plot.tc <- plot_tc_repeats(df, ylab)
  # plot mean, sd, ci95, and save stats for the tc
  plot.tc.err <- compute_mean_error(df, filename, ylab, show.linear.model=FALSE)
  # plot mean, sd, ci95, and save stats for the tc. Also add linear model information
  plot.tc.err.abline <- compute_mean_error(df, filename, ylab, show.linear.model=TRUE)

  
  # Comment if legend is inserted
  plot.tc <- plot.tc + theme(legend.position="none")
  # Uncomment if legend is inserted. This adds extra margins
  #plot.tc.err <- plot.tc.err + theme(plot.margin=unit(c(0.3,3.8,0.3,0.8),"cm"))
  #plot.tc.err.abline <- plot.tc.err.abline + theme(plot.margin=unit(c(0.3,3.8,0.3,0.8),"cm"))
  
  PLOTS <-list(plots=c(list(g1=plot.tc, g2=plot.tc.err, g3=plot.tc.err.abline)))
  
  plot.arrange <- do.call(grid.arrange, c(PLOTS$plots, nrow=3, bottom=paste0(filename, suffix)))
  # Comment if legend is inserted 
  ggsave(file.path(paste0(filename, '_all.png')), plot=plot.arrange, width=6, height=18, dpi=300)    
  # Uncomment if legend is inserted 
  #ggsave(file.path(paste0(filename, '_all.png')), plot=plot.arrange, width=8, height=18, dpi=300)  
  return(plot.arrange)
}


# Filter the data set
data_filtering <- function(df, remove.cols, remove.row.head, remove.row.tail) {
  # remove the unwanted columns
  df.filt <- df[, !(names(df) %in% c('Time'))]
  df.filt <- df.filt[, !(names(df.filt) %in% remove.cols)]
  
  # remove first and last entries
  df.filt <- df.filt[1:remove.row.tail, ]
  df.filt <- df.filt[remove.row.head:nrow(df.filt), ]  
  
  # log my data.
  # FluorescenceIntensity/ProteinConcentration is a sigmoid curve, but the central 
  # part, which represents our data, is nearly-linear. Therefore, we do not need to 
  # log our data. 
  # df.filt <- log10(df.filt)
  
  # apply min max (DISCARD min-max rescaling)
  #df.filt <- data.frame(Time=10 * 0:(nrow(df.filt)-1),
  #                      apply(df.filt, 2, normalise), 
  #                      check.names=FALSE)
  
  df.filt <- data.frame(Time=10 * 0:(nrow(df.filt)-1),
                        df.filt, 
                        check.names=FALSE)  
  
  return (df.filt)
}





# plot the table as heatmap.Rows are sorted by maximum increasing
# :param filename: the filename containing the table to plot
tc_heatmap <- function(folder, filename, labCol, title="Time courses", df.thres=17) {
  df <- read.table(file.path(folder, paste0(filename)), header=TRUE, 
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
  repeats <- nrow(df)
  
  # reduce the data
  if(ncol(df) > df.thres) {
    df <- df[1:df.thres,]
  }
  
  # plot the heatmap
  
  # creates a 5 x 5 inch image
  png(paste(file_path_sans_ext(filename), ".png", sep=""),    # create PNG for the heat map
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  # normalise each row within [0,1]
  df <- t(apply(df, 1, normalise))

  # replace NA with 0
  df[is.na(df)] <- 0
  
  par(cex.main=1.5, cex.axis=0.9)
  g <- heatmap.2(as.matrix(df),
                 main = title, # heat map title
                 col=matlab.like(256),
                 scale="none",
                 dendrogram="none", Rowv = FALSE, Colv = FALSE,
                 density.info="none",  # turn off density plot inside color legend
                 key=FALSE, symkey=FALSE, # turn off key                 
                 trace="none",        # turn off trace lines inside the heat map
                 labRow=FALSE,        # turn off row labels
                 margins = c(7, 4),   # reduce the right margin
                 lwid=c(4,25), lhei=c(4,25), # adjust the margin when key=FALSE
                 
                 labCol=labCol,
                 srtCol=45,
                 cexCol=2.5,
                 cexRow=1.6
                 # don't use this now as I haven't found a way to increase this fonts..
                 #xlab="Time (s)", ylab="Repeats by peak time"
  )
  # This works but is not elegant..
  mtext("Time [s]", side=1, line=4, cex=2.8)
  mtext(paste0("normalised repeats (n=",repeats, ")"), side=4, line=1, cex=2.8)
  
  # close the PNG device  
  dev.off()
  return(g)
}



### STATISTICS ###

# https://uk.mathworks.com/help/stats/lognrnd.html?requestedDomain=true 

logspace_mu <- function(mu, v) {
  log((mu^2)/sqrt(v+mu^2))
}


logspace_sd <- function(mu, v) {
  sqrt(log(v/(mu^2)+1))
}


skewness <- function(x, na.rm = FALSE, ...) {
  if (na.rm) x = x[!is.na(x)]
  return(sum((x-mean(x))^3/length(x))/sqrt(var(x))^3)
}


kurtosis <- function(x, na.rm = FALSE, ...) {
  if (na.rm) x = x[!is.na(x)]
  return(sum((x-mean(x))^4/length(x))/sqrt(var(x))^4)
}



# Basic scatter plot
scatter_plot <- function(X, Y, se=TRUE, annot.size=5.5, annot.x=6.5, annot.y=1) {
  df <- data.frame(X, Y)
  # linear model for my data
  fit <- lm(Y ~ X, data = df)
  digits <- 2
  g <- ggplot(df, aes(x=X, y=Y, group=X)) +
    geom_point() +
    geom_smooth(method=lm, se=se, aes(group=1)) +
    annotate("text", x=annot.x, y=annot.y, label = equation(fit, digits=digits), parse=TRUE, size=annot.size) +
    annotate("text", x=annot.x, y=annot.y-0.04, label = corr.coef(X, Y, digits=digits), parse=TRUE, size=annot.size) +
    theme_basic(24)

  return(g)
}


# Basic violin plot
violin_plot <- function(X, Y) {
  df <- data.frame(X, Y)
  g <- ggplot(df, aes(x=factor(X), y=Y, group=X)) +
    geom_violin(trim = FALSE, size=0.3) + 
    geom_jitter(height = 0, width = 0.1, size=0.5) +
    theme_basic(24)
  return(g)
}

# box violin plot
box_plot <- function(X, Y) {
  df <- data.frame(X, Y)
  g <- ggplot(df, aes(x=factor(X), y=Y, group=X)) +
    geom_boxplot() + 
    theme_basic(24)
  return(g)
}


###  for normal distribution test

# QQ plot
qq_plot <- function(vec, distribution=stats::qnorm) {
  d <- data.frame(residuals=vec)
  g <- ggplot(d, aes(sample=residuals)) + 
    stat_qq(distribution=distribution) + 
    stat_qq_line(distribution=distribution) +
    labs(title='Normal Q-Q') +
    theme_basic(24)
  return(g)
} 


# density plot
density_plot <- function(vec) {
  df <- data.frame(residuals=vec)
  mu <- data.frame(mu=mean(vec))
  g <- ggplot(df, aes(x=residuals)) +
    geom_density(colour = "black", fill = "#56B4E9", alpha=0.2) + 
    geom_vline(data=mu, aes(xintercept=mu), linetype="dashed") +
    theme_basic(24)
  return(g)
}


# density plot with colour
density_plot_wcolour <- function(vec, colour) {
  df <- data.frame(value=vec, variable=factor(colour))
  g <- ggplot(df, aes(x=value, group=variable, col=variable, fill=variable)) +
    stat_density(alpha=0.25) + 
    theme_basic(24)
  return(g)
}


# test whether the kurtosis is significantly different from zero.
# not normal if > 1
kurtosis.test <- function (x) {
  m4 <- sum((x-mean(x))^4)/length(x)
  s4 <- var(x)^2
  kurt <- (m4/s4) - 3
  sek <- sqrt(24/length(x))
  totest <- kurt/sek
  pvalue <- pt(totest,(length(x)-1))
  pvalue 
}

# test whether the skewness is significantly different from zero.
# not normal if > 1
skew.test <- function (x) {
  m3 <- sum((x-mean(x))^3)/length(x)
  s3 <- sqrt(var(x))^3
  skew <- m3/s3
  ses <- sqrt(6/length(x))
  totest <- skew/ses
  pt(totest,(length(x)-1))
  pval <- pt(totest,(length(x)-1))
  pval
}