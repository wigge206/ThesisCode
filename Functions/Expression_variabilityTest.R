library(lawstat)
library(limma)
library(hgu133plus2.db)

## Function to preform dvar and DE, pvalue correction is FDR
ExpVar.test <- function(data, classifier, p.adj = "BH"){
  Ftest = apply(data,1, function(x) var.test(x[classifier == 1],x[classifier == 0]))
  LevTest = apply(data, 1, function(x) levene.test(x, classifier, location="median", correction.method="zero.correction")) 
  
  df = data.frame(row.names = rownames(data), F = unlist(lapply(Ftest, function(x) x$statistic)),
                  p.value = unlist(lapply(Ftest, function(x) x$p.value)),
                  lev.t.statistic = unlist(lapply(LevTest, function(x) x$statistic)),
                  lev.p.value = unlist(lapply(LevTest, function(x) x$p.value)))
  df$lev.adj.p.value = p.adjust(df$lev.p.value, method = p.adj, n = nrow(df))
  
  
  classifier = ifelse(classifier == 1, "G1","G2")
  fgroup = factor(classifier)
  
  design = model.matrix(~ fgroup + 0, data)
  colnames(design) = levels(fgroup)
  rownames(design) = colnames(data)
  
  fit = lmFit(data, design)
  cont.matrix = makeContrasts(G1-G2, levels=design)
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2, 0.01)
  x = topTable(fit2, adjust=p.adj, sort.by="B", number=nrow(df),confint = T)
  x = x[match(row.names(df), row.names(x)),]
  colnames(x)[5:8] = paste("limma", colnames(x)[3:6], sep =".")
  
  df = cbind(df, x)
  
  return(df)
}