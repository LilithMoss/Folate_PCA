###############################
# Calculate PC's
# Look at validity of PC's
# Test those PC's against
# Case/Control Status
###############################
#Load Libraries
library(survival)
pca.output.path <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Folate/DATA_CLEANING/Isaac's_Data/Data/PCA Analysis/"

#Read in Matched case/control data
dat <- read.csv("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Folate/DATA_CLEANING/Isaac's_Data/Data/dat.epi.csv")

#Identify which columns are actually the variables
names(dat) #18 - 4242

#Isolate the variables
vars <- dat[,18:4242] 
vars.pca <- prcomp(vars,center = TRUE,scale. = TRUE)


#How much of the variance is captures by the PC's
jpeg("PCAs.jpg")
plot(vars.pca, type = "l")
dev.off()

summary(vars.pca) #Over 95% is captures by the first 96 PCs
write.csv(summary(vars.pca)$importance,file=paste0(pca.output.path,"pca.summay.csv"))

#Visualize PCA's
biplot(vars.pca,scale=0)

#How to regress the PCs on outcome? Logistic?
pcas <- vars.pca$x
dat.pc <- cbind(dat,pcas)
write.csv(dat.pc, "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Folate/DATA_CLEANING/Isaac's_Data/Data/dat.epi.pca.csv",row.names=F)

#Analyze via adjusted regression
d.pc <- dat.pc[,c(1:17,4243:4454)]
#Regress individual PCs 1st 96 of them
#Adjusted Model ############################################
results.adjusted <- do.call(rbind,lapply(20:115,function(i){
  reg <- clogit(case ~ d.pc[,i] + GFR + strata(matchID), data=d.pc)
  #reg <- clogit(case ~ d.pc[,i] + GFR + strata(matchID), data=d.pc)
  beta <- summary(reg)$coefficients[1,]
  or <- exp(cbind(OR = coef(reg), confint(reg)))[1,]
  #out <- cbind(colnames(dat)[i],beta,or)
  out <- cbind(t(beta),t(or))
}))
res.temp.adjusted <- as.data.frame(results.adjusted)
res.temp.adjusted$PC <- colnames(d.pc[,20:115])
res.df.adjusted <- res.temp.adjusted[,c(9,1:8)] #Re-adjust for easy viewing (PCs first)
res.df.sort.adjusted <- res.df.adjusted[order(res.df.adjusted[,6]),] #Sort by p-value
write.csv(res.df.sort.adjusted,file=paste0(pca.output.path,"pca.results.csv"),row.names=F)
############################################################

#Regress all PCs in one model
PCs <- formula( d.pc[,c(2,18:115)] ) 
X <- d.pc[,c(2,18:115)]
reg.all <- clogit(case~GFR+strata(matchID)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40+PC41+PC42+PC43+PC44+PC45+PC46+PC47+PC48+PC49+PC50+PC51+PC52+PC53+PC54+PC55+PC56+PC57+PC58+PC59+PC60+PC61+PC62+PC63+PC64+PC65+PC66+PC67+PC68+PC69+PC70+PC71+PC72+PC73+PC74+PC75+PC76+PC77+PC78+PC79+PC80+PC81+PC82+PC83+PC84+PC85+PC86+PC87+PC88+PC89+PC90+PC91+PC92+PC93+PC94+PC95+PC96,data=d.pc)
reg.all <- clogit(case~GFR+strata(matchID)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35,data=d.pc)

############################################################
#Try smaller number of variables
#Isolate the variables
vars <- dat[,18:25] 
vars.pca <- prcomp(vars,center = TRUE,scale. = TRUE)

#How much of the variance is captures by the PC's
jpeg("PCAs.jpg")
plot(vars.pca, type = "l")
dev.off()

#Use dat.pc to evaluate PC1,PC2,PC36, and PC5
#PC1
#Adjusted Model ############################################
results.adjusted <- do.call(rbind,lapply(18:4242,function(i){
  reg <- lm(PC1 ~ dat.pc[,i], data=dat.pc)
  beta <- summary(reg)$coefficients[1,]
}))
res.temp.adjusted <- as.data.frame(results.adjusted)
res.temp.adjusted$biomarker <- colnames(dat[,18:4242])
res.df.adjusted <- res.temp.adjusted[,c(5,1:4)]
res.df.sort.adjusted <- res.df.adjusted[order(res.df.adjusted[,5]),]
write.csv(res.df.sort.adjusted,file=paste0(pca.output.path,"PC1.analysis.csv"),row.names=F)


#Use dat.pc to evaluate PC1,PC2,PC36, and PC5
#PC36
#Adjusted Model ############################################
results.adjusted <- do.call(rbind,lapply(18:4242,function(i){
  reg <- lm(PC36 ~ dat.pc[,i], data=dat.pc)
  beta <- summary(reg)$coefficients[1,]
}))
res.temp.adjusted <- as.data.frame(results.adjusted)
res.temp.adjusted$biomarker <- colnames(dat[,18:4242])
res.df.adjusted <- res.temp.adjusted[,c(5,1:4)]
res.df.sort.adjusted <- res.df.adjusted[order(res.df.adjusted[,5]),]
write.csv(res.df.sort.adjusted,file=paste0(pca.output.path,"PC36.analysis.csv"),row.names=F)


