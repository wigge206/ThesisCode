# Code calls a function (Global.stats) which calculates mean, sd, cv and mean for each group (defined by binary vector)
# e.g. A vector is required where by a 1 symbols the first group and 0 the second. order matches samples in df

source("~/Projects/PhD/Bioinformatics/Scripts/Functions/Global_statistics.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Waddell global stats
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

load("Analysis_GEO/Expression_variation_WaddellData/Data/Geo_data/CompleteData.RData")
# select dataset
Waddell = GSE19177.filter
GSE19177_ch = GSE19177_ch[!grepl("unknown", GSE19177_ch)]

Waddell = Waddell[,match(names(GSE19177_ch), colnames(Waddell))]
Waddell = as.data.frame(Waddell)

normalised.exp.b1 = Waddell[, grepl("BRCA1|non-BRCA", GSE19177_ch)]
normalised.exp.b2 = Waddell[, grepl("BRCA2|non-BRCA", GSE19177_ch)]

classifer.basal = ifelse(grepl("Basal", GSE19177_ch), 1, 0)
classifer.b1 = GSE19177_ch[grepl("BRCA1|non-BRCA", GSE19177_ch)]
classifer.b1 = ifelse(grepl(": BRCA1", classifer.b1), 1, 0)
classifer.b2 = GSE19177_ch[grepl("BRCA2|non-BRCA", GSE19177_ch)]
classifer.b2 = ifelse(grepl("BRCA2", classifer.b2), 1, 0)

Waddell.global.basal = Global.stats(Waddell, classifer.basal, ID1 = "Basal", ID2 = "Non.basal")
Waddell.global.brca1 = Global.stats(normalised.exp.b1, classifer.b1, ID1 = "BRCA1", ID2 = "BRCAx")
Waddell.global.brca2 = Global.stats(normalised.exp.b2, classifer.b2, ID1 = "BRCA2", ID2 = "BRCAx")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Nagel gloabl stats 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
load("~/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_Nagel/Geo_data/CompleteData.RData")

affy.entrez = as.list(hgu133plus2SYMBOL)

## Pick which dataset to use, RMA, mas5 or RMA with probes filtered.
Nagel = GSE27830.filter
Nagel_subtype = read.csv("~/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_Nagel/Nagel_pam50classification.csv")

Nagel.b1.bx = Nagel[,grepl("Non|BRCA1", GSE27830_ch)]
Nagel.b2.bx = Nagel[,grepl("Non|BRCA2", GSE27830_ch)] ## only 6 brca2 samples

classifer.basal = ifelse(Nagel_subtype$pam50.classification.subtype == "Basal", 1, 0)
classifer.b1 = ifelse(grepl("BRCA1", GSE27830_ch[grepl("Non|BRCA1", GSE27830_ch)]), 1, 0)
classifer.b2 = ifelse(grepl("BRCA2", GSE27830_ch[grepl("Non|BRCA2", GSE27830_ch)]), 1, 0)

Nagel.global.basal = Global.stats(Nagel, classifer.basal, ID1 = "Basal", ID2 = "Non.basal")
Nagel.global.brca1 = Global.stats(Nagel.b1.bx, classifer.b1, ID1 = "BRCA1", ID2 = "BRCAx")
Nagel.global.brca2 = Global.stats(Nagel.b2.bx, classifer.b2, ID1 = "BRCA2", ID2 = "BRCAx")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Larsen global stats 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

load("Analysis_GEO/Expression_variation_Larsen/Geo_data/GSE49481/CompleteData.Rdata")
Probe_ID = read.delim("~/Projects/Microarray data/Larsen/GSE40115/Discovery_GPL15931/Probe_ID.txt")
PAM50 <-unlist(lapply(GSE49481_ch, function(x) lapply(strsplit(as.character(x[grepl("pam50agilent:",x)]), ": "), "[", 2)))
BRCA <- unlist(lapply(GSE49481_ch, function(x) lapply(strsplit(as.character(x[grepl("group:",x)]), ": "), "[", 2)))

Larsen = GSE49481_CY5quantile_filtered[,match(names(GSE49481_ch), colnames(GSE49481_CY5quantile_filtered))]

## Subset df so only columns that are required for tests are included
Larsen.b1.bx = Larsen[, c(grepl("BRCA1$|non-BRCA1/2", BRCA))]
Larsen.b2.bx = Larsen[,c(grepl("BRCA2$|non-BRCA1/2",  BRCA))]
Larsen.b1.sporadic = Larsen[, c(grepl("BRCA1$|Sporadic", BRCA))]
Larsen.b2.sporadic = Larsen[,c(grepl("BRCA2$|Sporadic",  BRCA))]

classifer.basal = ifelse(grepl("Basal", PAM50), 1, 0)
classifer.b1.bx = ifelse(grepl("BRCA1$", BRCA[grepl("BRCA1$|non-BRCA1/2", BRCA)]), 1, 0)
classifer.b2.bx = ifelse(grepl("BRCA2$",BRCA[grepl("BRCA2$|non-BRCA1/2", BRCA)]), 1, 0)
classifer.b1.sporadic = ifelse(grepl("BRCA1$", BRCA[grepl("BRCA1$|Sporadic", BRCA)]), 1, 0)
classifer.b2.sporadic = ifelse(grepl("BRCA2$", BRCA[grepl("BRCA2$|Sporadic", BRCA)]), 1, 0)

Larsen.global.basal = Global.stats(Larsen, classifer.basal, ID1 = "Basal", ID2 = "Non.basal")
Larsen.global.brca1.bx = Global.stats(Larsen.b1.bx, classifer.b1.bx, ID1 = "BRCA1", ID2 = "BRCAx")
Larsen.global.brca2.bx = Global.stats(Larsen.b2.bx, classifer.b2.bx, ID1 = "BRCA2", ID2 = "BRCAx")
Larsen.global.brca1.sporadic = Global.stats(Larsen.b1.sporadic, classifer.b1.sporadic, ID1 = "BRCA1", ID2 = "Sporadic")
Larsen.global.brca2.sporadic = Global.stats(Larsen.b2.sporadic, classifer.b2.sporadic, ID1 = "BRCA2", ID2 = "Sporadic")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### bc2116 metacohort global stats
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
affy.entrez = as.list(hgu133plus2SYMBOL)
bc2116clin <- read.csv('~/Projects/PhD/Bioinformatics/BRCA2116/BRCA2116data_v2/BC2116-clinicaldata.csv')
load('~/Projects/PhD/Bioinformatics/BRCA2116/BRCA2116data_v2/BC2116-rma-noqnorm-COMBAT.RData')

classifer.basal = ifelse(bc2116clin$Subtype == "Basal", 1, 0)

bc2116.global.basal = Global.stats(bc2116.combat, classifer.basal, ID1 = "Basal", ID2 = "Non.basal")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Save all expression data to one dataframe
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

save(list = c("Waddell.global.basal","Waddell.global.brca1","Waddell.global.brca2","Nagel.global.basal", "Nagel.global.brca1","Nagel.global.brca2", "Larsen.global.basal", "Larsen.global.brca1.sporadic", "Larsen.global.brca1.bx", "Larsen.global.brca2.sporadic", "Larsen.global.brca2.bx", "bc2116.global.basal"), file = "All_globalStats.RData")
