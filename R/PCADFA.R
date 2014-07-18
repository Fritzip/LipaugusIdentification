setwd("C:/Users/Bernadette/Desktop/Maxime/Lipaugus/LipaugusIdentification/output/")
library(ade4)
library(MASS)
library(class)
library(scatterplot3d)
data = read.table(file="measures3.dat",header=F,sep=' ')
out = which(data$V2=='Y' | data$V2=='X' | data$V2=='b' | data$V2=='d' | data$V2=='c' | data$V2=='j')
data = data[-out,]
which(data$V2=='Y')
data = data[c(-36,-100),]
measures = data[,c(-1,-2)]
measures[is.na(measures)] <- 0
acp1 = dudi.pca(measures, center = T, scale = T)


scatter(acp1, posieig = "none", xax = 1, yax = 2)
s.label(acp1$li, xax=1, yax=2)
s.class(acp1$li,data$V2)

###
Classes = factor(data$V2)
RawData = acp1$l1

for(i in 1:13) {
  Temp <- RawData[,i]
  a <- mean(Temp)
  b <- var(Temp)
  RawData[,i] <- (Temp-a)/b
}
head(RawData)

train <- sample(1:dim(RawData)[1],70)

result1 <- knn(RawData[train,],RawData[-train,],Classes[train],k=1)
table1 <- table(old=Classes[-train],new=result1)
table1
chisq.test(table1)


temp <- lda(Classes ~ . , data=cbind(Classes,RawData), scale=TRUE, subset=train)
temp$scaling
temp
palette(c("black" ,  "red",     "green3",  "blue" ,   "cyan",    "magenta" ,"yellow",  "gray", "purple"))
transformed <- as.matrix(RawData)%*%temp$scaling
par(mfrow=c(1,3))
plot(transformed,pch = 16)
plot(transformed[,2],transformed[,3],pch = 16)
plot(transformed[,1],transformed[,3],pch = 16)
par(mfrow=c(1,1))
scatterplot3d(x=transformed[,3],y=transformed[,2],z=transformed[,1],pch=16,color=as.numeric(Classes))

result2 <- knn(transformed[train,],transformed[-train,],Classes[train],k=1)
table2 <- table(old=Classes[-train],new=result2)
table2
chisq.test(table2)
