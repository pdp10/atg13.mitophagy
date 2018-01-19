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



# https://uk.mathworks.com/help/stats/lognrnd.html?requestedDomain=true 

meanlog <- function(mu, v) {
  log((mu^2)/sqrt(v+mu^2))
}


sdlog <- function(mu, v) {
  sqrt(log(v/(mu^2)+1))
}


# Return the skewness. For a normal distribution, skewness is 0.
skewness <- function(x, na.rm = FALSE, ...) {
  if (na.rm) x = x[!is.na(x)]
  return(sum((x-mean(x))^3/length(x))/sqrt(var(x))^3)
}


# Return the kurtosis excess (kurtosis - 3). For a normal distribution, kurtosis excess is 0.
kurtosis <- function(x, na.rm = FALSE, ...) {
  if (na.rm) x = x[!is.na(x)]
  kurt <- sum((x-mean(x))^4/length(x))/sqrt(var(x))^4 
  return(kurt - 3)
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
  pval <- pt(totest,(length(x)-1))
  pval
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


# min max normalisation
# call: apply(df, 2, normalise) 
# to normalise each column within [0,1]
normalise <- function(x, na.rm = TRUE) {
  ranx <- range(x, na.rm = na.rm)
  (x - ranx[1]) / diff(ranx)
}









####################
## Statistical plots
####################

# QQ plot of a vector against a distribution 
qqplot <- function(vec, distribution=stats::qnorm) {
  d <- data.frame(resids=vec)
  g <- ggplot(d, aes(sample=resids)) + 
    stat_qq(distribution=distribution) + 
    stat_qq_line(distribution=distribution)
}


# overlay histogram, normal distribution
hist_w_norm <- function(vec) {
  df <- data.frame(x=vec)
  g <- ggplot(df, aes(x)) +
    geom_histogram(bins=length(df$x)/3, colour = "black", fill = "white", alpha=0.2) + 
    stat_function(fun = dnorm, 
                  args = list(mean = mean(df$x), sd = sd(df$x)), 
                  lwd = 1,
                  colour='red')
}


# overlay histogram, normal and lognormal distributions
hist_w_dist <- function(vec) {
  df <- data.frame(x=vec)
  g <- ggplot(df, aes(x)) +
    geom_histogram(aes(y = ..density..), bins=length(df$x), color = "black", fill = 'white') +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(df$x), sd = sd(df$x)), 
                  lwd = 1,
                  aes(colour='norm')) +
    stat_function(fun = dlnorm, 
                  args = list(meanlog = meanlog(mean(df$x), var(df$x)), sdlog = sdlog(mean(df$x), var(df$x))), 
                  lwd = 1, 
                  aes(colour='lnorm')) +
    scale_colour_manual('distrib', values=c("blue", "red"))
}



# overlay histogram, coloured density
hist_w_coldensity <- function(vec, colour) {
  df <- data.frame(value=vec, variable=factor(colour))
  g <- ggplot(df, aes(x=value, group=variable, col=variable, fill=variable)) +
    geom_histogram(bins=length(df$value)/6, colour = "black", fill = "white", alpha=0.2) +
    stat_density(alpha=0.25)
}


plot_correlation <- function(df, xlab, ylab) {
  #print(df)
  pcorr <- cor(df[,c(xlab)], df[,c(ylab)])  
  #print(pcorr)
  g <- ggplot(df, aes_string(xlab, ylab)) +
    geom_point() +     # Use hollow circles
    geom_smooth(method=lm) +  # Add linear regression line (by default includes 95% confidence region)
    ggtitle(paste0("R=", round(pcorr, digits = 3)))
}


# box plot
box_plot <- function(X, Y) {
  df <- data.frame(X, Y)
  g <- ggplot(df, aes(x=factor(X), y=Y, group=X)) +
    geom_boxplot()
}
