setwd("C:/Users/Bernadette/Desktop/Maxime/Lipaugus/LipaugusIdentification/")
library(ade4)
data = read.table(file="measures2.dat",header=F,sep=' ')
which(data$V2=='Y')
which(data$V2=='X')
data = data[c(-36,-100),]
measures = data[,c(-1,-2)]
head(measures)
#colnames(data) = c('a0','a1','b1','a2','b2','w','max','length','energy')
acp1 = dudi.pca(measures, center = T, scale = T)
5
#par(mfrow=c(1,3))
scatter(acp1, posieig = "none", xax = 1, yax = 2)
scatter(acp1, posieig = "none", xax = 2, yax = 3)
scatter(acp1, posieig = "none", xax = 1, yax = 3)
s.label(acp1$li, xax=1, yax=2)
s.label(acp1$li, xax=2, yax=3)
s.label(acp1$li, xax=1, yax=3)
#par(mfrow=c(1,1))

s.class(acp1$li,data$V2)

###
Classes = factor(data$V2)
RawData = acp1$l1

### LDA

# Linear Discriminant Analysis with Jacknifed Prediction
library(MASS)
fit <- lda(as.factor(data$V2) ~ x1 + x2 + x3, data=newdata,
           na.action="na.omit", CV=TRUE)
fit # show results 



Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                   Sp = rep(c("s","c","v"), rep(50,3)))
train <- sample(1:150, 75)
table(Iris$Sp[train])
## your answer may differ
##  c  s  v
## 22 23 30
z <- lda(Sp ~ ., Iris, prior = c(1,1,1)/3, subset = train)
predict(z, Iris[-train, ])$class
##  [1] s s s s s s s s s s s s s s s s s s s s s s s s s s s c c c
## [31] c c c c c c c v c c c c v c c c c c c c c c c c c v v v v v
## [61] v v v v v v v v v v v v v v v
(z1 <- update(z, . ~ . - Petal.W.))

library(class)

Raw <- read.csv(url("http://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data"),header=FALSE)
head(Raw)

Classes <- Raw[,1]
RawData <- Raw[,2:14]
for(i in 1:5) {
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
plot(transformed,col=Classes)
plot(transformed[,2],transformed[,3],col=Classes)

result2 <- knn(transformed[train,],transformed[-train,],Classes[train],k=1)
table2 <- table(old=Classes[-train],new=result2)
table2
chisq.test(table2)









tr <- sample(1:50, 25)
train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
test <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
cl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
z <- qda(train, cl)
table(predict(z,test)$class,cl)


