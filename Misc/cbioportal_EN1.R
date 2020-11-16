library("cgdsr")
library('ggplot2')
library(reshape2)
library(gridExtra)
library(lawstat)
mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)

genes <- c("IGF2BP3","EN1","DSC3","ESR1","BRCA1")
##Get available case lists (collection of samples) for a given cancer study
mycancerstudy =  'brca_metabric'
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]

# Get available genetic profiles
getGeneticProfiles(mycgds,mycancerstudy)
#select genetic profile
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[c(1),1]
mygeneticprofile.rna = getGeneticProfiles(mycgds,mycancerstudy)[c(3),1]

# Get data slices for a specified list of genes, genetic profile and case list
METBRIC_CNA<-getProfileData(mycgds,c(genes),mygeneticprofile,mycaselist)
METBRIC_RNA<-getProfileData(mycgds,c(genes),mygeneticprofile.rna,mycaselist)
METABRIC_clinical <- getClinicalData(mycgds,mycaselist)
METABRIC_clinical<- METABRIC_clinical[match(rownames(METBRIC_RNA), rownames(METABRIC_clinical)),]

## exclude DSC3 as no expression
p.df <- cbind(METBRIC_RNA[,c("EN1", "IGF2BP3")],METABRIC_clinical)
tmp <- melt(p.df, id.vars = colnames(METABRIC_clinical))

##calculate levene.test pvalue
basal.t.f <- ifelse(p.df$CLAUDIN_SUBTYPE[!is.na(p.df$IGF2BP3)]== "Basal",1,0)
EN1.basal.p <- levene.test(p.df$EN1[!is.na(p.df$EN1)], ifelse(p.df$CLAUDIN_SUBTYPE[!is.na(p.df$EN1)]== "Basal",1,0), location='median')
IGF2BP3.basal.p <- levene.test(p.df$IGF2BP3[!is.na(p.df$IGF2BP3)], ifelse(p.df$CLAUDIN_SUBTYPE[!is.na(p.df$IGF2BP3)]== "Basal",1,0), location='median')


subtype_meta<- ggplot(tmp[!is.na(tmp$value),], aes(x=CLAUDIN_SUBTYPE, y=value))+
  geom_boxplot(fill='#B0C4De', alpha=0.7, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha=0.5, shape=16, size=0.8)+theme_bw(base_size = 10) + 
  facet_grid(~variable) +ylab("RNA expression (Z-Scores)") + xlab(NULL)+ theme(strip.text = element_text(face = "italic"))

png(filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/EN1_IGF2BP3_subtype.png", width=7, height=2,units='in', res=300)
subtype_meta +theme(axis.title = element_text(size=8))
dev.off()

EN1_erStatus <- ggplot(tmp[!is.na(tmp$value)& tmp$variable=="EN1",], aes(x=ER_STATUS, y=value))+
  geom_boxplot(fill='#B0C4De', alpha=0.7, outlier.shape = NA, size=0.1)+ xlab("ER status") +
  geom_jitter(position=position_jitter(0.2), alpha=0.5, shape=16, size = 0.1)+
  theme_bw(base_size = 6) +ylab(expression(paste(italic("EN1 "), "expression")))
levene.test(p.df$EN1[!is.na(p.df$EN1)], p.df$ER_STATUS[!is.na(p.df$EN1)], location='median')
#2.18726e-106

png(filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/EN1_ER_bx.png", width=1, height=1,units='in', res=300)
EN1_erStatus +theme(axis.title = element_text(size=4))
dev.off()


IGF2BP3_erStatus <- ggplot(tmp[!is.na(tmp$value)& tmp$variable=="IGF2BP3",], aes(x=ER_STATUS, y=value))+
  geom_boxplot(fill='#B0C4De', alpha=0.7, outlier.shape = NA, size=0.1)+ xlab("ER status") +
  geom_jitter(position=position_jitter(0.2), alpha=0.5, shape=16, size=0.1)+theme_bw(base_size = 6) +
  ylab(expression(paste(italic("IGF2BP3 "), "expression")))
levene.test(p.df$IGF2BP3[!is.na(p.df$IGF2BP3)], p.df$ER_STATUS[!is.na(p.df$IGF2BP3)], location='median')
#2.868254e-87
png(filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/IGF2BP3_ER_bx.png", width=1, height=1,units='in', res=300)
IGF2BP3_erStatus +theme(axis.title = element_text(size=4))
dev.off()


#ggsave(subtype_meta, filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/EN1_IGF2BP3_subtype.png", width=7, height=3)
#ggsave(EN1_erStatus, filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/EN1_ER_bx.png", width=4, height=4)
#ggsave(IGF2BP3_erStatus, filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/IGF2BP3_ER_bx.png", width=4, height=4)





EN1.esr1<-ggplot(METBRIC_RNA, aes(x=ESR1, y=EN1))+geom_point(col='#B0C4De', size=0.8, alpha=0.8) +theme_bw(base_size = 10) + xlab(expression(paste(italic("ESR1 "), "expression (Z-Scores)"))) + 
ylab(expression(paste(italic("EN1 "), "expression (Z-Scores)"))) +theme(axis.title = element_text(size=8))

png("C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/EN1_esr1_scatter.png",width=3.5, height=2.5, units='in',res =300 )
EN1.esr1
dev.off()

METBRIC_RNA$median <- ifelse(METBRIC_RNA$ESR1 > median(METBRIC_RNA$ESR1, na.rm=T), "Greater","Less")
META.complete.obs <- METBRIC_RNA[!is.na(METBRIC_RNA$EN1)& !is.na(METBRIC_RNA$ESR1),]
META.complete.obs$median<-factor(META.complete.obs$median, levels=c("Less", "Greater"))
med.bx<-ggplot(META.complete.obs, aes(x=median, y=EN1))+geom_boxplot(fill='#B0C4De', alpha=0.7, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha=0.5, shape=16)+theme_bw(base_size = 16) + ylab("EN1 expression")+
  xlab("ESR1 expression (split by median)")

IGF2BP3.esr1 <- ggplot(METBRIC_RNA, aes(x=ESR1, y=IGF2BP3))+geom_point(col='#B0C4De', size=0.8, alpha=0.8) +theme_bw(base_size = 10) + xlab("ESR1 expression") + 
  ylab(expression(paste(italic("IGF2BP3 "), "expression (Z-Scores)"))) +theme(axis.title = element_text(size=8))
png("C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/IGF2BP3_esr1_scatter.png",width=3.5, height=2.5, units='in',res =300 )
IGF2BP3.esr1
dev.off()

META.complete.obs <- METBRIC_RNA[!is.na(METBRIC_RNA$IGF2BP3)& !is.na(METBRIC_RNA$ESR1),]
META.complete.obs$median<-factor(META.complete.obs$median, levels=c("Less", "Greater"))
igf.med.bx<-ggplot(META.complete.obs, aes(x=median, y=IGF2BP3))+geom_boxplot(fill='#B0C4De', alpha=0.7, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha=0.5, shape=16)+theme_bw(base_size = 16) + ylab("IGF2BP3 expression")+
  xlab("ESR1 expression (split by median)")

ggsave(igf.med.bx, filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/IGF2BP3_esr1_medi.png", width=5, height=5)
ggsave(med.bx, filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/EN1_esr1_medi.png", width=5, height=5)
ggsave(EN1.esr1, filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/EN1_esr1_scatter.png", width=5.5, height=4)
ggsave(IGF2BP3.esr1, filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/IGF2BP3_esr1_scatter.png", width=5.5, height=4)



## Boxplot of CN number
EN1 <- data.frame(row.names=rownames(METBRIC_CNA),CN=METBRIC_CNA$EN1, RNA=METBRIC_RNA$EN1)
EN1 <- EN1[!is.na(EN1$RNA),]
EN1$CN[EN1$CN == -1] <-"Deletion"
EN1$CN[EN1$CN == 0] <-"Diploid"
EN1$CN[EN1$CN == 1] <-"Gain"
EN1$CN[EN1$CN == 2] <-"Amplification"

EN1$CN <- factor(EN1$CN, levels= c("Deletion", "Diploid", "Gain", "Amplification"))

EN1_cn<-ggplot(EN1, aes(y=RNA, x=CN,group=CN)) +geom_boxplot(fill='#B0C4De', alpha=0.7, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha=0.5, shape=16, size=0.8)+theme_bw(base_size = 10) +xlab("Copy number variation") +
  ylab(expression(paste(italic("EN1 "), "RNA expression (Z-Scores)"))) +theme(axis.title = element_text(size=8))

IGF2BP3 <- data.frame(row.names=rownames(METBRIC_CNA),CN=METBRIC_CNA$IGF2BP3, RNA=METBRIC_RNA$IGF2BP3)
IGF2BP3 <- IGF2BP3[!is.na(IGF2BP3$RNA),]
IGF2BP3$CN[IGF2BP3$CN == -1] <-"Deletion"
IGF2BP3$CN[IGF2BP3$CN == 0] <-"Diploid"
IGF2BP3$CN[IGF2BP3$CN == 1] <-"Gain"
IGF2BP3$CN[IGF2BP3$CN == 2] <-"Amplification"

IGF2BP3$CN <- factor(IGF2BP3$CN, levels= c("Deletion", "Diploid", "Gain", "Amplification"))

IGF2BP3_cn<-ggplot(IGF2BP3, aes(y=RNA, x=CN,group=CN)) +geom_boxplot(fill='#B0C4De', alpha=0.7, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha=0.5, shape=16, size=0.8)+theme_bw(base_size = 10) +xlab("Copy number variation") +
  ylab(expression(paste(italic("IGF2BP3 "), "RNA expression (Z-Scores)"))) +theme(axis.title = element_text(size=8))

png(filename = "C:/UOC/gwiggins/Projects/PhD/Thesis Figures/ExpressionVar/EN1_IGF2BP3_CNV.png", width=7, height=4,units='in', res=300)
grid.arrange(EN1_cn,IGF2BP3_cn, nrow=2)
dev.off()



subtype_meta<- ggplot(tmp[!is.na(tmp$value),], aes(x=CLAUDIN_SUBTYPE, y=value))+
  geom_boxplot(fill='#B0C4De', alpha=0.7, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha=0.5, shape=16, size=0.8)+theme_bw(base_size = 10) + 
  facet_grid(~variable) +ylab("RNA expression (Z-Scores)") + xlab(NULL)+ theme(strip.text = element_text(face = "italic"))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### TCGA nature
mycancerstudy <- "brca_tcga_pub"
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
#select genetic profile
mygeneticprofile.cna = getGeneticProfiles(mycgds,mycancerstudy)[c(2),1]
mygeneticprofile.rna = getGeneticProfiles(mycgds,mycancerstudy)[c(10),1]

TCGA_pub_cna<-getProfileData(mycgds,c(genes),mygeneticprofile,mycaselist)
TCGA_pub_rna<-getProfileData(mycgds,c(genes),mygeneticprofile.rna,mycaselist)

ggplot(TCGA_pub_rna, aes(x=ESR1, y=EN1))+geom_point(col='#B0C4De') +theme_bw(base_size = 16) + xlab("ESR1 expression") + 
  ylab("EN1 expression")


ggplot(TCGA_pub_rna, aes(x=ESR1, y=IGF2BP3))+geom_point(col='#B0C4De') +theme_bw(base_size = 16) + xlab("ESR1 expression") + 
  ylab("IGF2BP3 expression")

ggplot(TCGA_pub_rna, aes(x=ESR1, y=DSC3))+geom_point(col='#B0C4De') +theme_bw(base_size = 16) + xlab("ESR1 expression") + 
  ylab("DSC3 expression")
