library(DOSE)
library(clusterProfiler)
library(illuminaHumanv2.db)
## Load data
load("~/Projects/PhD/Bioinformatics/ExpressionAnalysis_listform.Rdata")

Lar.probe = read.delim("~/Projects/Microarray data/Larsen/GPL15931-probeEntrez_annotation.txt")

entrez = as.list(hgu133plus2ENTREZID)
ilmn.entrez = as.list(illuminaHumanv2ENTREZID)

## Create background genes list
bc2116.universal = unlist(entrez[names(entrez) %in% rownames(basal$bc2116.basal)])
wad.universal = unlist(ilmn.entrez[names(ilmn.entrez) %in% basal$Waddell.basal$IlluminaID])
Lar.universal = as.character(Lar.probe$entrez) ## imported Lar probe from GEO
nag.universal = unlist(entrez[names(entrez) %in% rownames(basal$Nagel.basal)])

## Pathway enrichment for basal vs non-basal breast tumours
basal.sig<-lapply(basal, function(i) i[i$lev.adj.p.value < 0.05,])


signif_entrez.bc2116 = entrez[names(entrez) %in% rownames(basal.sig$bc2116.basal)]
GO.bc2116 = enrichGO(signif_entrez.bc2116, 'org.Hs.eg.db', readable =T, universe = bc2116.universal,ont ="BP", pvalueCutoff = 1, qvalueCutoff = 1)

signif_entrez.wad = ilmn.entrez[names(ilmn.entrez) %in% basal.sig$Waddell.basal$IlluminaID]
GO.wad = enrichGO(signif_entrez.wad, 'org.Hs.eg.db', readable =T, universe = wad.universal,ont ="BP", pvalueCutoff = 1, qvalueCutoff = 1)

signif_entrez.lar = Lar.probe$entrez[match(basal.sig$Nagel.basal$ID,Lar.probe$ID)]
GO.lar = enrichGO(signif_entrez.lar, 'org.Hs.eg.db', readable =T, universe = Lar.universal,ont ="BP", pvalueCutoff = 1, qvalueCutoff = 1)

signif_entrez.nag = entrez[names(entrez) %in% rownames(basal.sig$Nagel_Dvar)]
GO.nag = enrichGO(signif_entrez.nag, 'org.Hs.eg.db', readable =T, universe = nag.universal,ont ="BP", pvalueCutoff = 1, qvalueCutoff = 1)
