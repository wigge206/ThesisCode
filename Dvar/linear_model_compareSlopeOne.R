### Comparing linear model to one
## Generate a function that calculates T stat from lm beta compared to 1 and s.e.of beta
## P-value is two sided
ttest <- function(reg, coefnum, val){
  co <- coef(summary(reg))
  tstat <- (co[coefnum,1]-val)/co[coefnum,2]
  2 * pt(abs(tstat), reg$df.residual, lower.tail = FALSE)
}

load("C:/UOC/gwiggins/Projects/PhD/Bioinformatics/All_globalStats.RData")

slope <- c()
CI.lower <- c()
CI.upper <- c()

df<-mget(ls(patter="global"))
models <- list()
for(dat in 1:length(df)){
  reg <- list()
  for (i in 1:((ncol(df[[dat]]))/2)){
    reg[[i]] <- lm(df[[dat]][,i+4]~df[[dat]][,i])
    
  }
  names(reg) <- c("Mean","sd", "mad","cv")
  models[[dat]] <- reg
}
names(models) <- names(df)


## Calculate p-value
model.p <- list()
for(dat in 1:length(models)){
  reg.p <- list()
  for(i in 1:length(models[[dat]])){
    reg.p[[i]] <- ttest(models[[dat]][[i]], 2, 1)
  }
  names(reg.p) <- c("Mean","sd", "mad","cv")
  model.p[[dat]] <- unlist(reg.p)
}
names(model.p) <- names(df)

betas <- lapply(models, function(i) unlist(lapply(i, function(i) coef(summary(i))[2,1])))
lower.CI <- lapply(models, function(i) unlist(lapply(i, function(i) confint(i)[2,1])))
upper.CI <- lapply(models, function(i) unlist(lapply(i, function(i) confint(i)[2,2])))
beta.se <- lapply(models, function(i) unlist(lapply(i, function(i) coef(summary(i))[2,2])))


## Make df
df <- list()
for (i in 1:4){
  df[[i]] <- data.frame(t(data.frame(betas))[,i],t(data.frame(lower.CI))[,i],t(data.frame(upper.CI))[,i],t(data.frame(beta.se))[,i],t(data.frame(model.p))[,i])
}
df<- do.call(cbind, df)
colnames(df) <- paste(rep(c("Mean", "SD","MAD","CV"), each=5), c("beta", "lower.CI","upper.CI", "SE", "p"), sep="_")
