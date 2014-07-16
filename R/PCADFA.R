setwd("C:/Users/Bernadette/Desktop/Maxime/Lipaugus/LipaugusIdentification/")
library(ade4)
library(MASS)
library(class)
data = read.table(file="measures2.dat",header=F,sep=' ')
which(data$V2=='Y')
which(data$V2=='X')
data = data[c(-36,-100),]
measures = data[,c(-1,-2)]
n = 30
acp1 = dudi.pca(measures, center = T, scale = T)
30

scatter(acp1, posieig = "none", xax = 1, yax = 2)
s.label(acp1$li, xax=1, yax=2)
s.class(acp1$li,data$V2)

###
Classes = factor(data$V2)
RawData = acp1$l1

for(i in 1:n) {
  Temp <- RawData[,i]
  a <- mean(Temp)
  b <- var(Temp)
  RawData[,i] <- (Temp-a)/b
}
head(RawData)

train <- sample(1:dim(RawData)[1],80)

result1 <- knn(RawData[train,],RawData[-train,],Classes[train],k=1)
table1 <- table(old=Classes[-train],new=result1)
table1
chisq.test(table1)


temp <- lda(Classes ~ . , data=cbind(Classes,RawData), scale=TRUE, subset=train)
temp$scaling
temp

transformed <- as.matrix(RawData)%*%temp$scaling
plot(transformed,col=Classes,lwd = 5)
plot(transformed[,2],transformed[,3],col=Classes)

result2 <- knn(transformed[train,],transformed[-train,],Classes[train],k=1)
table2 <- table(old=Classes[-train],new=result2)
table2
chisq.test(table2)
