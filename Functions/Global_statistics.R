Global.stats = function(data, classifer, ID1 = "Group1", ID2 = "Group2"){
  sd.g1 = apply(data[,classifer == 1],1,sd)
  sd.g2 = apply(data[,classifer == 0],1,sd)
  mean.g1 = apply(data[,classifer == 1],1,mean)
  mean.g2 = apply(data[,classifer == 0],1,mean)
  mad.g1 = apply(data[,classifer == 1],1,mad)
  mad.g2 = apply(data[,classifer == 0],1,mad)
  cv.g1 = apply(data[,classifer == 1],1, function(x) sd(x)/mean(x)*100)
  cv.g2 = apply(data[,classifer == 0],1, function(x) sd(x)/mean(x)*100)
  
  df =data.frame(row.names = rownames(data), mean.g1, sd.g1, mad.g1, cv.g1, mean.g2, sd.g2, mad.g2, cv.g2)
  colnames(df) = c(paste(c("mean", "sd", "mad", "cv"), ID1, sep = "_"), paste(c("mean", "sd", "mad", "cv"), ID2, sep = "_"))
  return(df)
}
