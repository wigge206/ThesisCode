## Function scatterplot generate scatter plots for Mean, SD, CV and MAD
## Input data is probe level stat
library(ggplot2)
library(reshape2)
library(gridExtra)
scatterplot = function(data, stat = c("sd", "mad", "cv", "mean"), model = "loess", Legend = T, ...){
  data = data[,grepl(stat, colnames(data))]
  plot(data[,1], data[,2], col =rgb(102,102,102,50,maxColorValue=255),pch=16, bg ='white', yaxs = "i", xaxs = "i", ...)
  abline(0,1, col = "red", lwd=2)
  if (model == "lm"){
    abline(lm(data[,2]~data[,1]), col = 'blue', lwd=2)
    model.type = "linear regression"
  } else if (model == "lowess"){
    lines(lowess(data[,1], data[,2]), col = "blue", lwd=2)
    model.type = "LOWESS"
  } else if (model == 'loess'){
    lines(lowess(data[,1], data[,2]), col = "blue", lwd=2)
    model.type = "LOESS"
  } else {
    lines(lowess(data[,1], data[,2]), col = "skyblue", lwd=2)
    abline(lm(data[,2]~data[,1]), col = 'blue', lty = 5, lwd=2)
    model.type = c("LOESS", "linear regression")
  }
  if(Legend & length(model.type) == 1 ){
    legend('top', lty = 1, col= c("red", "blue"), legend = c("No variance", model.type))
  } else if (Legend){
    legend('top', lty = c(1,1,5), col= c("red", "blue", "blue"), legend = c("No variance", model.type))
  }
}

## Load data
load("C://UOC/gwiggins/Projects/PhD/Bioinformatics/All_globalStats.RData")

all.df <- mget(ls(pattern="global"))

outfile <- gsub("\\.global\\.","_", names(all.df))
outfile <- gsub("\\.","_", outfile)
outfile <-paste("C://UOC/gwiggins/Projects/PhD/Bioinformatics/Plots/Global_var/",outfile,sep="")
xlab <- c("Basal", "Basal", "BRCA1", "BRCA1", "BRCA2","BRCA2", "Basal", "BRCA1","BRCA2", "Basal", "BRCA1", "BRCA2")
ylab <- c("Non-basal", "Non-basal", "BRCAx", "Spordaic", "BRCAx","Sporadic", "Non-basal", "BRCAx","BRCAx", "Non-basal", "BRCAx", "BRCAx")

for(i in 1:length(all.df)){
  xlab.p <- paste("Mean (", xlab[i], ")", sep="")
  ylab.p <- paste("Mean (", ylab[i], ")", sep="" )
  png(paste(outfile[i],"mean.png",sep="_"), res=300, width = 2.5, height =2.5,  units ='in', pointsize=9)
  scatterplot(all.df[[i]], stat="mean", Legend = F, model = 'both', xlab=xlab.p, ylab=ylab.p, cex=.5)
  dev.off()
  
  xlab.p <- paste("SD (", xlab[i], ")", sep="")
  ylab.p <- paste("SD (", ylab[i], ")", sep="" )
  png(paste(outfile[i],"SD.png",sep="_"), res=300, width = 2.5, height =2.5,  units ='in', pointsize=9)
  scatterplot(all.df[[i]], stat="sd", Legend = F, model = 'both',  xlab=xlab.p, ylab=ylab.p, cex=.5)
  dev.off()
  
  xlab.p <- paste("MAD (", xlab[i], ")", sep="")
  ylab.p <- paste("MAD (", ylab[i], ")", sep="" )
  png(paste(outfile[i],"MAD.png",sep="_"), res=300, width = 2.5, height =2.5,  units ='in', pointsize=9)
  scatterplot(all.df[[i]], stat="mad", Legend = F, model = 'both',  xlab=xlab.p, ylab=ylab.p, cex=.5)
  dev.off()
  
  xlab.p <- paste("CV (", xlab[i], ")", sep="")
  ylab.p <- paste("CV (", ylab[i], ")", sep="" )
  png(paste(outfile[i],"cv.png",sep="_"), res=300, width = 2.5, height =2.5,  units ='in', pointsize=9)
  scatterplot(all.df[[i]], stat="cv", Legend = F, model = 'both',  xlab=xlab.p, ylab=ylab.p, cex=.5)
  dev.off()
 
}  
  

