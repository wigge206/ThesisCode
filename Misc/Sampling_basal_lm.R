library(ggplot)
set.seed(123)
nonVar <- rnorm(1000, mean=10, sd=1)
V <- rnorm(1000, mean =10, sd=5)

predicted.sdr <- sd(V)/sd(nonVar)

dat <- list()
for(i in 10:1000){
  
  dat[[i-9]] <- replicate(n=100, expr = mean(sd(sample(V,i))/mean(sd(sample(nonVar,i)))))
  
}

df.sdr <- data.frame(sdr=unlist(lapply(dat,mean)), sdr.sd=unlist(lapply(dat,sd)), samples =10:1000)

ggplot(df.sdr, aes(samples, sdr)) +geom_point() +theme_bw() +geom_smooth(method="loess", se=F)+
  ylab("Standard deviation ratio (SDR)") + xlab("Sample size") + 
  geom_hline(yintercept=predicted.sdr, color='red', size =1, linetype='dashed')

## Save as PDF for vector graphic for line
ggsave("SDR_simulation.pdf", height=4, width=8)


bc2116clin <- read.csv('C:/UOC/gwiggins/projects/PhD/Bioinformatics/BRCA2116/BRCA2116data_v2/BC2116-clinicaldata.csv')
load('C:/UOC/gwiggins/projects/PhD/Bioinformatics/BRCA2116/BRCA2116data_v2/BC2116-rma-noqnorm-COMBAT.RData')

basal.id <- bc2116clin$SampleID[bc2116clin$Subtype == "Basal"]
Nonbasal.id <- bc2116clin$SampleID[bc2116clin$Subtype != "Basal"]
basal <- bc2116.combat[, as.character(basal.id)]
non.basal <- bc2116.combat[, as.character(Nonbasal.id)]

predicted.lm <- lm(apply(non.basal,1,sd)~apply(basal,1,sd))[[1]][2]

dat.lm <- list()
dat.lm.c <- list()
for(i in 10:490){
  print(i-9)
  dat.lm[[i-9]] <- replicate(n=100, expr = {
    basal.sd <- apply(basal[,c(sample(1:490,i))],1,sd) 
    nonbasal.sd <- apply(non.basal[,c(sample(1:1626,i))],1,sd) 
    
    lm(nonbasal.sd~basal.sd)[[1]][2]
  })
}

my_y <- expression(paste("Linear model (", beta, ")"))
df <- data.frame(lm.coeff=unlist(lapply(dat.lm,mean)), samples =10:490)
ggplot(df, aes(samples, lm.coeff)) +geom_point() +theme_bw() +geom_smooth(method="loess")+
  ylab(my_y) + xlab("Sample size") + 
  geom_hline(yintercept=predicted.lm, color='red', size =1, linetype='dashed') +
  theme(legend.position= c(0.8,0.2))
ggsave("LM_sampling.pdf", height=4, width=8)
