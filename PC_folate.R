###############################
# Calculate PC's
# Look at validity of PC's
# Test those PC's against
# Case/Control Status
###############################
#Load Libraries

#Read in Matched case/control data
dat <- read.csv("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/Folate/DATA_CLEANING/Isaac's_Data/Data/dat.epi.csv")

#Identify which columns are actually the variables
names(dat) #18 - 4242

#Isolate the variables
vars <- dat[,18:4242] 
vars.pca <- prcomp(vars,center = TRUE,scale. = TRUE)

#HOw much of the variance is captures by the PCA's
jpeg("PCAs.jpg")
plot(vars.pca, type = "l")
dev.off()

summary(vars.pca) #Over 95% is captures by the first 96 PCAs

#Visualize PCA's
biplot(vars.pca)
