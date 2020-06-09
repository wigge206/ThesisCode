# Code calls a function which runs var.test, levene test (lawstat) and DE (limma) based on a binary definer
# e.g. A vector is made where by a 1 symbols the first group and 0 the second. order matches samples in df
library(lawstat)
library(limma)
## Function to preform dvar and DE, pvalue correction is FDR by default
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
  x = topTable(fit2, adjust=p.adj, sort.by="B", number=nrow(df))
  x = x[match(row.names(df), row.names(x)),]
  colnames(x)[3:6] = paste("limma", colnames(x)[3:6], sep =".")
  
  df = cbind(df, x)
  
  return(df)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Waddell expression stats (DVAR AND DE)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##### To reduce multiple testing only transcripts with symbols
## comparison made
# Basal vs non-basal
# BRCA1 vs BRCAx
# BRCA2 vs BRCAx


load("Analysis_GEO/Expression_variation_WaddellData/Data/Geo_data/CompleteData.RData")
rm(GSE19177_expr)

## pick dataset quantile normalised (Log2), probes filtered
Waddell = GSE19177.filter

## remove samples that is "unknown genotype from expr and charateristic
GSE19177_ch = GSE19177_ch[!grepl("unknown", GSE19177_ch)]
Waddell = Waddell[,match(names(GSE19177_ch), colnames(Waddell))]
Waddell = as.data.frame(Waddell)
GPL6106 <- GPL6106[match(rownames(Waddell), GPL6106$ID),]
Waddell <- Waddell[GPL6106$Symbol != "",]

normalised.exp.b1 = Waddell[, grepl("BRCA1", GSE19177_ch)]
normalised.exp.basal.b1 = Waddell[, grepl("Basal", GSE19177_ch)]
normalised.exp.b2 = Waddell[, grepl("BRCA2|non", GSE19177_ch)]
normalised.exp.b1v2 = Waddell[, grepl(": BRCA1|BRCA2", GSE19177_ch)]


classifer.basal = ifelse(grepl("Basal", GSE19177_ch),1,0)
classifer.b1 = ifelse(grepl(": BRCA1", GSE19177_ch[grepl("BRCA1", GSE19177_ch)]),1,0)
classifer.basal.b1 = ifelse(grepl(": BRCA1", GSE19177_ch[grepl("Basal", GSE19177_ch)]),1,0) # B1.basal vs non-b1.basal (include b2)
classifer.b2 = ifelse(grepl("BRCA2", GSE19177_ch[grepl("BRCA2|non", GSE19177_ch)]),1,0)
classifer.b1.2 = ifelse(grepl(": BRCA1|BRCA2", GSE19177_ch),1,0) # B1/2 vs non-b1/2
classifer.b1v2 = ifelse(grepl("BRCA1", GSE19177_ch[grepl(": BRCA1|BRCA2", GSE19177_ch)]),1,0) # B1vB2


Waddell.basal = ExpVar.test(Waddell, classifer.basal)
Waddell.brca1 = ExpVar.test(normalised.exp.b1, classifer.b1)
Waddell.basal.brca1 = ExpVar.test(normalised.exp.basal.b1, classifer.basal.b1)
Waddell.brca2 = ExpVar.test(normalised.exp.b2, classifer.b2)
Waddell.brca1.2 = ExpVar.test(Waddell, classifer.b1.2)
Waddell.brca1v2 = ExpVar.test(normalised.exp.b1v2, classifer.b1v2)

Waddell.basal = data.frame(GPL6106[match(rownames(Waddell.basal), GPL6106$ID), c("Search_key", "IlluminaID","Transcript", "Symbol")], Waddell.basal, row.names = rownames(Waddell.basal))
Waddell.brca1 = data.frame(GPL6106[match(rownames(Waddell.brca1), GPL6106$ID), c("Search_key", "IlluminaID", "Transcript", "Symbol")], Waddell.brca1, row.names = rownames(Waddell.basal))
Waddell.basal.brca1 = data.frame(GPL6106[match(rownames(Waddell.basal.brca1), GPL6106$ID), c("Search_key", "IlluminaID", "Transcript", "Symbol")], Waddell.basal.brca1, row.names = rownames(Waddell.basal))
Waddell.brca2 = data.frame(GPL6106[match(rownames(Waddell.brca2), GPL6106$ID), c("Search_key", "IlluminaID", "Transcript", "Symbol")], Waddell.brca2, row.names = rownames(Waddell.basal))
Waddell.brca1.2 = data.frame(GPL6106[match(rownames(Waddell.brca1.2), GPL6106$ID), c("Search_key", "IlluminaID", "Transcript", "Symbol")], Waddell.brca1.2, row.names = rownames(Waddell.basal))
Waddell.brca1v2 = data.frame(GPL6106[match(rownames(Waddell.brca1v2), GPL6106$ID), c("Search_key", "IlluminaID", "Transcript", "Symbol")], Waddell.brca1v2, row.names = rownames(Waddell.basal))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Nagel expression stats (DVAR AND DE)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##### To reduce multiple testing only transcripts with symbols
## comparison made
# Basal vs non-basal
# BRCA1 vs BRCAx
# BRCA2 vs BRCAx

load("~/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_Nagel/Geo_data/CompleteData.RData")
affy.entrez = as.list(hgu133plus2SYMBOL)
#Nagel =  read.delim("~/Projects/Microarray data/Nagel/RMA_Nagel.txt", row.names =1)
Nagel_subtype = read.csv("~/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_Nagel/Nagel_pam50classification.csv", row.names= 1)

## Pick which dataset to use, RMA, mas5 or RMA with probes filtered.
Nagel = GSE27830.filter

## remove CHEK2
Nagel = Nagel[,!grepl("CHEK2",GSE27830_ch)]
Nagel_subtype = Nagel_subtype[!grepl("CHEK2", Nagel_subtype$X),]
GSE27830_ch = GSE27830_ch[!grepl("CHEK2",GSE27830_ch)]

affy.entrez <- affy.entrez[match(rownames(Nagel), names(affy.entrez))]
Nagel <- Nagel[!is.na(affy.entrez),]

Nagel.b1.bx = Nagel[,grepl("Non|BRCA1", GSE27830_ch)]
Nagel.basal.b1 = Nagel[,grepl("Basal", Nagel_subtype$pam50.classification.subtype)]
Nagel.b2.bx = Nagel[,grepl("Non|BRCA2", GSE27830_ch)] ## only 6 brca2 samples
Nagel.b1vb2 = Nagel[,grepl("BRCA1|BRCA2", GSE27830_ch)] ## only 6 brca2 samples

classifer.basal = ifelse(Nagel_subtype$pam50.classification.subtype == "Basal", 1, 0)
classifer.b1 = ifelse(grepl("BRCA1", GSE27830_ch[grepl("Non|BRCA1", GSE27830_ch)]), 1, 0)
classifer.basal.b1 = ifelse(grepl("BRCA1", GSE27830_ch[grepl("Basal", Nagel_subtype$pam50.classification.subtype)]), 1, 0)
classifer.b2 = ifelse(grepl("BRCA2", GSE27830_ch[grepl("Non|BRCA2", GSE27830_ch)]), 1, 0)
classifer.b1.2 = ifelse(grepl("BRCA1|BRCA2", GSE27830_ch), 1, 0)
classifer.b1v2 = ifelse(grepl("BRCA1", GSE27830_ch[grepl("BRCA1|BRCA2", GSE27830_ch)]), 1, 0)

Nagel.basal = ExpVar.test(Nagel, classifer.basal)
Nagel.brca1 = ExpVar.test(Nagel.b1.bx, classifer.b1)
Nagel.basal.brca1 = ExpVar.test(Nagel.basal.b1, classifer.basal.b1)
Nagel.brca2 = ExpVar.test(Nagel.b2.bx, classifer.b2)
Nagel.brca1.brca2 = ExpVar.test(Nagel, classifer.b1.2)
Nagel.brca1vbrca2 = ExpVar.test(Nagel.b1vb2, classifer.b1v2)

Nagel.basal = cbind(unlist(affy.entrez[match(rownames(Nagel.basal), names(affy.entrez))]),Nagel.basal)
Nagel.brca1 = cbind(unlist(affy.entrez[match(rownames(Nagel.brca1), names(affy.entrez))]),Nagel.brca1)
Nagel.basal.brca1 = cbind(unlist(affy.entrez[match(rownames(Nagel.basal.brca1), names(affy.entrez))]),Nagel.basal.brca1)
Nagel.brca2 = cbind(unlist(affy.entrez[match(rownames(Nagel.brca2), names(affy.entrez))]),Nagel.brca2)
Nagel.brca1.brca2 = cbind(unlist(affy.entrez[match(rownames(Nagel.brca1.brca2), names(affy.entrez))]),Nagel.brca1.brca2)
Nagel.brca1vbrca2 = cbind(unlist(affy.entrez[match(rownames(Nagel.brca1vbrca2), names(affy.entrez))]),Nagel.brca1vbrca2)

colnames(Nagel.basal)[1] = "Gene.symbol"
colnames(Nagel.brca1)[1] = "Gene.symbol"
colnames(Nagel.basal.brca1)[1] = "Gene.symbol"
colnames(Nagel.brca2)[1] = "Gene.symbol"
colnames(Nagel.brca1.brca2)[1] = "Gene.symbol"
colnames(Nagel.brca1vbrca2)[1] = "Gene.symbol"


#write.csv(Nagel.brca1, "~/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_Nagel/BRCA1_expression_Stats.csv")
#write.csv(Nagel.brca2, "~/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_Nagel/BRCA2_expression_Stats.csv")
#write.csv(Nagel.basal, "~/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_Nagel/Basal_expression_Stats.csv")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Larsen expression stats (DVAR AND DE)
## There is two lots of data, the original (GSE40115) and the later which includes BRCAx samples (GSE49481). See GSE40115vs GSE49481.Rmd
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## comparison made
# Basal vs non-basal
# BRCA1 vs BRCAx
# BRCA2 vs BRCAx
# BRCA1 vs Sporadic
# BRCA2 vs Sporadic

## Load either or both datasets
load("Analysis_GEO/Expression_variation_Larsen/Geo_data/GSE49481/CompleteData.Rdata")
Probe_ID = read.delim("~/Projects/Microarray data/Larsen/GSE40115/Discovery_GPL15931/Probe_ID.txt")
PAM50 <-unlist(lapply(GSE49481_ch, function(x) lapply(strsplit(as.character(x[grepl("pam50agilent:",x)]), ": "), "[", 2)))
BRCA <- unlist(lapply(GSE49481_ch, function(x) lapply(strsplit(as.character(x[grepl("group:",x)]), ": "), "[", 2)))

## Select which dataframe to use (GSE40115_CY5_quantile_filtered or GSE40115_exprs)
# GSE40115_exprs - Larsen et al normalised, NB: may contain NA's for probes that were filtered out
# GSE40115_CY5_quantile_filtered - CY5 intensties, background corrected, qunatile normalised, highest mean expressing probe selected for duplicated probes 
Larsen <- GSE49481_CY5quantile_filtered
Larsen <- Larsen[, match(names(GSE49481_ch), colnames(Larsen))]

## Subset df so only columns that are required for tests are included
Larsen.b1.bx = Larsen[, c(grepl("BRCA1$|non-BRCA1/2", BRCA))]
Larsen.basal.b1 = Larsen[, grepl("Basal", PAM50)]
Larsen.b2.bx = Larsen[,c(grepl("BRCA2$|non-BRCA1/2",  BRCA))]
Larsen.b1.2.bx = Larsen[, c(!grepl("Sporadic", BRCA))]
Larsen.b1v2 = Larsen[, c(grepl("BRCA1$|BRCA2$", BRCA))]
Larsen.b1.sporadic = Larsen[, c(grepl("BRCA1$|Sporadic", BRCA))]
Larsen.b2.sporadic = Larsen[,c(grepl("BRCA2$|Sporadic",  BRCA))]
Larsen.b1.2.sporadic = Larsen[,c(!grepl("non-BRCA1/2",  BRCA))]


classifer.basal = ifelse(grepl("Basal", PAM50), 1, 0)
classifer.b1.bx = ifelse(grepl("BRCA1$", BRCA[grepl("BRCA1$|non-BRCA1/2", BRCA)]), 1, 0)
classifer.basal.b1 = ifelse(grepl("BRCA1$", BRCA[grepl("Basal", PAM50)]), 1, 0)
classifer.b2.bx = ifelse(grepl("BRCA2$", BRCA[grepl("BRCA2$|non-BRCA1/2", BRCA)]), 1, 0)
classifer.b1.2.bx = ifelse(grepl("BRCA1$|BRCA2", BRCA[!grepl("Sporadic", BRCA)]), 1, 0)
classifer.b1v2 = ifelse(grepl("BRCA1$", BRCA[grepl("BRCA1$|BRCA2$", BRCA)]), 1, 0)
classifer.b1.sporadic = ifelse(grepl("BRCA1$", BRCA[grepl("BRCA1$|Sporadic", BRCA)]), 1, 0)
classifer.b2.sporadic = ifelse(grepl("BRCA2$", BRCA[grepl("BRCA2$|Sporadic", BRCA)]), 1, 0)
classifer.b1.2.sporadic = ifelse(grepl("BRCA1$|BRCA2$", BRCA[!grepl("non-BRCA1/2", BRCA)]), 1, 0)

Larsen.basal = ExpVar.test(Larsen, classifer.basal)
Larsen.brca1.bx = ExpVar.test(Larsen.b1.bx, classifer.b1.bx)
Larsen.basal.b1 = ExpVar.test(Larsen.basal.b1, classifer.basal.b1) #
Larsen.brca2.bx = ExpVar.test(Larsen.b2.bx, classifer.b2.bx)
Larsen.b1.2.bx = ExpVar.test(Larsen.b1.2.bx, classifer.b1.2.bx) #
Larsen.b1v2 = ExpVar.test(Larsen.b1v2, classifer.b1v2) #
Larsen.brca1.sporadic = ExpVar.test(Larsen.b1.sporadic, classifer.b1.sporadic)
Larsen.brca2.sporadic = ExpVar.test(Larsen.b2.sporadic, classifer.b2.sporadic)
Larsen.b1.b2.sporadic = ExpVar.test(Larsen.b1.2.sporadic, classifer.b1.2.sporadic)

Larsen.basal = cbind(Probe_ID[match(row.names(Larsen.basal), Probe_ID$ID),c("ID", "GENE_SYMBOL")], Larsen.basal)
Larsen.brca1.bx = cbind(Probe_ID[match(row.names(Larsen.brca1.bx), Probe_ID$ID),c("ID", "GENE_SYMBOL")], Larsen.brca1.bx)
Larsen.basal.b1 = cbind(Probe_ID[match(row.names(Larsen.basal.b1), Probe_ID$ID),c("ID", "GENE_SYMBOL")], Larsen.basal.b1) 
Larsen.brca2.bx = cbind(Probe_ID[match(row.names(Larsen.brca2.bx), Probe_ID$ID),c("ID", "GENE_SYMBOL")], Larsen.brca2.bx)
Larsen.b1.2.bx = cbind(Probe_ID[match(row.names(Larsen.b1.2.bx), Probe_ID$ID),c("ID", "GENE_SYMBOL")], Larsen.b1.2.bx)
Larsen.b1v2 = cbind(Probe_ID[match(row.names(Larsen.b1v2), Probe_ID$ID),c("ID", "GENE_SYMBOL")], Larsen.b1v2)
Larsen.brca1.sporadic = cbind(Probe_ID[match(row.names(Larsen.brca1.sporadic), Probe_ID$ID),c("ID", "GENE_SYMBOL")], Larsen.brca1.sporadic)
Larsen.brca2.sporadic = cbind(Probe_ID[match(row.names(Larsen.brca2.sporadic), Probe_ID$ID),c("ID", "GENE_SYMBOL")], Larsen.brca2.sporadic)
Larsen.b1.b2.sporadic = cbind(Probe_ID[match(row.names(Larsen.b1.b2.sporadic), Probe_ID$ID),c("ID", "GENE_SYMBOL")], Larsen.b1.b2.sporadic)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### bc2116 metacohort expression stats (DVAR AND DE)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##### To reduce multiple testing only transcripts with symbols
## comparison made
# Basal vs non-basal

affy.entrez = as.list(hgu133plus2SYMBOL)
bc2116clin <- read.csv('~/Projects/PhD/Bioinformatics/BRCA2116/BRCA2116data_v2/BC2116-clinicaldata.csv')
load('~/Projects/PhD/Bioinformatics/BRCA2116/BRCA2116data_v2/BC2116-rma-noqnorm-COMBAT.RData')

classifer.basal = ifelse(bc2116clin$Subtype == "Basal", 1, 0)

affy.entrez <- affy.entrez[match(rownames(bc2116.combat), names(affy.entrez))]
bc2116.combat <- bc2116.combat[!is.na(affy.entrez),]

bc2116.basal = ExpVar.test(bc2116.combat, classifer.basal)
bc2116.basal = cbind(unlist(affy.entrez[match(rownames(bc2116.basal), names(affy.entrez))]),bc2116.basal)

colnames(bc2116.basal)[1] = "Gene.symbol"


#write.csv(bc2116.basal, "~/Projects/PhD/Bioinformatics/BRCA2116/basal_variation.txt")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Save all expression data to one dataframe
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

save(list = c("Waddell.basal","Waddell.brca1","Waddell.basal.brca1","Waddell.brca2", "Waddell.brca1.2", "Waddell.brca1v2",
              "Nagel.basal", "Nagel.brca1", "Nagel.basal.brca1","Nagel.brca2", "Nagel.brca1.brca2", "Nagel.brca1vbrca2",
              "Larsen.basal", "Larsen.brca1.bx","Larsen.basal.b1", "Larsen.brca2.bx","Larsen.b1.2.bx","Larsen.b1v2", "Larsen.brca1.sporadic", "Larsen.brca2.sporadic","Larsen.b1.b2.sporadic",
              "bc2116.basal"), file = "All_expressionAnalysis.RData")
