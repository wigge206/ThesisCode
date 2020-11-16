library(hgu133plus2.db)
load("C:/UOC/gwiggins/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_WaddellData/Data/Geo_data/CompleteData.RData")
rm(GSE19177_expr)
## pick dataset quantile normalised (Log2), probes filtered
Waddell = GSE19177.filter

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

SDR.b1Bx <- apply(normalised.exp.b1[classifer.b1==1],1,sd)/apply(normalised.exp.b1[classifer.b1==0],1,sd)
SDR.b2Bx <- apply(normalised.exp.b2[classifer.b2==1],1,sd)/apply(normalised.exp.b2[classifer.b2==0],1,sd)
SDR.basal <- apply(Waddell[classifer.basal==1],1,sd)/apply(Waddell[classifer.basal==0],1,sd)

WaddellSDR = data.frame(GPL6106[match(names(SDR.basal), GPL6106$ID), c("Search_key", "IlluminaID","Transcript", "Symbol")], BRCA1.Bx = SDR.b1Bx, 
                           BRCA2.bx = SDR.b2Bx, Basal=SDR.basal)

colnames(WaddellSDR)[4] <- "Gene.symbol"

##

load("C:/UOC/gwiggins/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_Nagel/Geo_data/CompleteData.RData")
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


Nag_SDR.b1Bx <- apply(Nagel.b1.bx[classifer.b1==1],1,sd)/apply(Nagel.b1.bx[classifer.b1==0],1,sd)
Nag_SDR.b2Bx <- apply(Nagel.b2.bx[classifer.b2==1],1,sd)/apply(Nagel.b2.bx[classifer.b2==0],1,sd)
Nag_SDR.basal <- apply(Nagel[classifer.basal==1],1,sd)/apply(Nagel[classifer.basal==0],1,sd)

NagelSDR = data.frame(Gene.symbol= unlist(affy.entrez[match(names(Nag_SDR.b1Bx), names(affy.entrez))]),BRCA1.Bx = Nag_SDR.b1Bx, 
                         BRCA2.bx = Nag_SDR.b2Bx, Basal=Nag_SDR.basal)


##
## Load either or both datasets
load("C:/UOC/gwiggins/Projects/PhD/Bioinformatics/Analysis_GEO/Expression_variation_Larsen/Geo_data/GSE49481/CompleteData.Rdata")
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
Larsen.b2.bx = Larsen[,c(grepl("BRCA2$|non-BRCA1/2",  BRCA))]
Larsen.b1.sporadic = Larsen[, c(grepl("BRCA1$|Sporadic", BRCA))]
Larsen.b2.sporadic = Larsen[,c(grepl("BRCA2$|Sporadic",  BRCA))]

classifer.basal = ifelse(grepl("Basal", PAM50), 1, 0)
classifer.b1.bx = ifelse(grepl("BRCA1$", BRCA[grepl("BRCA1$|non-BRCA1/2", BRCA)]), 1, 0)
classifer.b2.bx = ifelse(grepl("BRCA2$", BRCA[grepl("BRCA2$|non-BRCA1/2", BRCA)]), 1, 0)
classifer.b1.sporadic = ifelse(grepl("BRCA1$", BRCA[grepl("BRCA1$|Sporadic", BRCA)]), 1, 0)
classifer.b2.sporadic = ifelse(grepl("BRCA2$", BRCA[grepl("BRCA2$|Sporadic", BRCA)]), 1, 0)


lar_SDR.b1Bx <- apply(Larsen.b1.bx[classifer.b1.bx==1],1,sd)/apply(Larsen.b1.bx[classifer.b1.bx==0],1,sd)
lar_SDR.b1spor <- apply(Larsen.b1.sporadic[classifer.b1.sporadic==1],1,sd)/apply(Larsen.b1.sporadic[classifer.b1.sporadic==0],1,sd)
lar_SDR.b2Bx <- apply(Larsen.b2.bx[classifer.b2.bx==1],1,sd)/apply(Larsen.b2.bx[classifer.b2.bx==0],1,sd)
lar_SDR.b2spor <- apply(Larsen.b2.sporadic[classifer.b2.sporadic==1],1,sd)/apply(Larsen.b2.sporadic[classifer.b2.sporadic==0],1,sd)
lar_SDR.basal <- apply(Larsen[classifer.basal==1],1,sd)/apply(Larsen[classifer.basal==0],1,sd)

LarsenSDR = cbind(Probe_ID[match(names(lar_SDR.b1Bx), Probe_ID$ID),c("ID", "GENE_SYMBOL")], 
                  data.frame(BRCA1.Bx = lar_SDR.b1Bx, BRCA2.bx = lar_SDR.b2Bx, Basal=lar_SDR.basal, BRCA1.spor = lar_SDR.b1spor,
                             BRCA2.spor=lar_SDR.b2spor))
colnames(LarsenSDR)[2] <- "Gene.symbol"


## meta-cohort
affy.entrez = as.list(hgu133plus2SYMBOL)
bc2116clin <- read.csv('C:/UOC/gwiggins/Projects/PhD/Bioinformatics/BRCA2116/BRCA2116data_v2/BC2116-clinicaldata.csv')
load('C:/UOC/gwiggins/Projects/PhD/Bioinformatics/BRCA2116/BRCA2116data_v2/BC2116-rma-noqnorm-COMBAT.RData')

classifer.basal = ifelse(bc2116clin$Subtype == "Basal", 1, 0)

affy.entrez <- affy.entrez[match(rownames(bc2116.combat), names(affy.entrez))]
bc2116.combat <- bc2116.combat[!is.na(affy.entrez),]

meta.basal <- apply(bc2116.combat[classifer.basal==1],1,sd)/apply(bc2116.combat[classifer.basal==0],1,sd)
metaSDR = data.frame(Gene.symbol= unlist(affy.entrez[match(names(meta.basal), names(affy.entrez))]),Basal=meta.basal)

SDR.list <- list(WaddellSDR, NagelSDR, LarsenSDR, metaSDR)
names(SDR.list) <- c("Waddell.SDR", "Nagel.SDR", "Larsen.SDR", "Meta-cohort.SDR")

save(SDR.list, file= "C:/UOC/gwiggins/Projects/PhD/Bioinformatics/SDRs.Rdata")
