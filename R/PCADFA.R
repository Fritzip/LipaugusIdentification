setwd("C:/Users/Bernadette/Desktop/Maxime/Lipaugus/LipaugusIdentification/")
library(ade4)
library(MASS)
library(class)
library(scatterplot3d)

### Import data & Omit missing values ###
data = read.table(file="measures.dat",header=F,sep=' ')
out = which(data$V2=='X')
data = data[-out,]
data = na.omit(data)
rawmeasures = data[,c(-1,-2)]

### Do PCA ###
n = 7 # Keep 7 axis (all)
acp1 = dudi.pca(rawmeasures, center = T, scale = T, scannf = F, nf=7)

### Plot PCA results ###
scatter(acp1, posieig = "none", xax = 1, yax = 2)
s.label(acp1$li, xax=1, yax=2)
s.class(acp1$li,data$V2)

### Create variables ###
groups = factor(data$V2)
measures = acp1$l1

### Train sample ###
len = dim(measures)[1]
train <- sample(1:len, len) # 100% to train - 0% to test

### Do LDA ###
temp <- lda(groups ~ . , data=cbind(groups,measures), scale=TRUE, subset=train)

### New parameter space ###
temp$scaling

### Measures in the new parameter space ###
transformed <- as.matrix(measures)%*%temp$scaling

### Plot LDA results ###
palette(c("black" ,  "red",     "green3",  "blue" ,   "cyan",    "magenta" ,"yellow",  "gray", "purple"))
par(mfrow=c(1,3))
plot(transformed, col=groups, pch = 16)
plot(transformed[,2], transformed[,3], col=groups, pch = 16)
plot(transformed[,1], transformed[,3], col=groups, pch = 16)
par(mfrow=c(1,1))
scatterplot3d(x=transformed[,3],y=transformed[,2],z=transformed[,1],pch=16,color=as.numeric(groups))

