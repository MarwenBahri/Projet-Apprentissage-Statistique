rm(list=objects());graphics.off()
setwd("D:/ENSTA/Cours/STA203/Projet/")

#Partie 2: Analyse exploratoire
#Part2.Question 1
load("cookie.app.RData")
load("cookie.val.RData")

xtrain <- as.matrix(cookie.app)
xtest <- as.matrix(cookie.val)

ytrain <- xtrain[,1]
ytest <- xtest[,1]


xtrain <- xtrain[,-1]
xtest <- xtest[,-1]


plot.new()
library(car)

boxplot(xtrain)
par(mfrow=c(2,2),mar=c(2,2,2,1))
boxplot(xtrain[,1:50])
boxplot(xtrain[,51:100])
boxplot(xtrain[,101:150])
boxplot(xtrain[,151:200])

boxplot(xtrain[,201:250])
boxplot(xtrain[,251:300])
boxplot(xtrain[,301:350])
boxplot(xtrain[,351:400])

boxplot(xtrain[,401:450])
boxplot(xtrain[,451:500])
boxplot(xtrain[,501:550])
boxplot(xtrain[,551:600])

par(mfrow=c(1,2),mar=c(6,1,6,1))
boxplot(xtrain[,601:650])
boxplot(xtrain[,651:700])

par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
matplot(x=1:700,y=t(xtrain),type="l",xlab="variables (fréquences)", ylab="absorption",lty=c(1,1))

C <- cor(xtrain)
C.eig <- eigen(C)
plot(C.eig$values,xlab="valeur propre",ylab="")
points(C.eig$values,col="red")
nb.vp <- length(which(C.eig$values >= 0.001)) #nombre de valeurs propres non nulles = 31 << p=700

#Part2.Question 2
library(FactoMineR)
res=PCA(xtrain)
nbr.vp <- nrow(res$eig) #nombre de valeurs propres
#graphe des valeurs propres
par(mfrow=c(1,1))
plot(res$eig[,1],xlab="valeur propre",ylab="")
points(res$eig[,1],pch=20,col="red")

barplot(res$eig[,2],main="% inertie",names=paste("Dim",1:nrow(res$eig)))
res <- PCA(xtrain,ncp=nbr.vp)
par(mfrow=c(1,2))
#premier plan principal
plot(res,axes = c(1,2), choix = "var",graph.type = "classic")
plot(res,axes = c(1,2), choix = "ind",graph.type = "classic")
#deuxième plan principal
plot(res,axes = c(2,3), choix = "var",graph.type = "classic")
plot(res,axes = c(2,3), choix = "ind",graph.type = "classic")
#troisième plan principal
plot(res,axes = c(3,4), choix = "var",graph.type = "classic")
plot(res,axes = c(3,4), choix = "ind",graph.type = "classic")
#quatrième plan principal
plot(res,axes = c(4,5), choix = "var",graph.type = "classic")
plot(res,axes = c(4,5), choix = "ind",graph.type = "classic")
#cinquième plan principal
plot(res,axes = c(5,6), choix = "var",graph.type = "classic")
plot(res,axes = c(5,6), choix = "ind",graph.type = "classic")

#Part2.Question 3 
reconstruct<-function(res,nr,Xm,Xsd)
{
  coord.var = as.matrix(res$var$coord)[,1:nr,drop=F]
  coord.ind = as.matrix(res$ind$coord)[,1:nr,drop=F]
  x_rec = coord.ind%*%t(sweep(coord.var,2,sqrt(res$eig[1:nr,1]),FUN="/"))
  x_rec = sweep(x_rec,2,Xsd,FUN="*")
  x_rec = sweep(x_rec,2,Xm,FUN="+")
  return(x_rec)
}
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)

res <- PCA(xtrain,ncp=39)
library("Metrics")
mat=as.matrix(xtrain)
Xm <- apply(mat,2,mean)
Xsd <- apply(mat,2,sd)

Xrec1=reconstruct(res,1,Xm,Xsd)
Xrec2=reconstruct(res,2,Xm,Xsd)
Xrec3=reconstruct(res,3,Xm,Xsd)
Xrec4=reconstruct(res,4,Xm,Xsd)
Xrec5=reconstruct(res,5,Xm,Xsd)
Xrec39=reconstruct(res,39,Xm,Xsd)

par(mfrow=c(3,2),mar=c(2,4,4,2)+0.1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 1",type="l",main=paste("RMSE=",round(rmse(mat,Xrec1),5),"MAE=",round(mae(mat,Xrec1),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec2),xlab="",ylab="nr = 2",type="l",main=paste("RMSE=",round(rmse(mat,Xrec2),5),"MAE=",round(mae(mat,Xrec2),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 3",type="l",main=paste("RMSE=",round(rmse(mat,Xrec3),5),"MAE=",round(mae(mat,Xrec3),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 4",type="l",main=paste("RMSE=",round(rmse(mat,Xrec4),5),"MAE=",round(mae(mat,Xrec4),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 5",type="l",main=paste("RMSE=",round(rmse(mat,Xrec5),5),"MAE=",round(mae(mat,Xrec5),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 39",type="l",main=paste("RMSE=",round(rmse(mat,Xrec39),5),"MAE=",round(mae(mat,Xrec39),5),collapse = " "),cex.main=1)

par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
matplot(Xrec1[,24],type="l",col="gray0",ylab="X24")
matlines(Xrec2[,24],col="red")
matlines(Xrec3[,24],col="darkblue")
matlines(Xrec4[,24],col="chartreuse")
matlines(Xrec5[,24],col="gray")
matlines(Xrec39[,24],col="darkmagenta")
legend("topleft",c("Xrec1","Xrec2","Xrec3","Xrec4","Xrec5","Xrec39"),
       col=c("gray0","red","darkblue","chartreuse","gray","darkmagenta"),lty=c(1:6),lwd = 0.5 , xpd = T,cex=0.7 )


#Partie 3: Régression pénalisée
#Part3.Question 1
library(glmnet)

grid=10^seq(6,-10,length=100) # la grille de lambda
ridge.glm=glmnet(xtrain,ytrain,alpha=0,lambda=grid)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)#default
plot(ridge.glm,xlab="L2 Norm") #coefficients de la regression ridge

coef(ridge.glm)[1,] #l'intercept
plot(log(grid),coef(ridge.glm)[1,],xlab="log(k)",ylab="intercept")
points(log(grid),coef(ridge.glm)[1,],col="red",pch=20)
abline(h=mean(ytrain),lwd=2,col="blue")
text(-20, mean(ytrain)+2, "mean(y)", col = "blue")

#recalcule de l'intercept
theta0 <- mean(ytrain)*rep(1,length(grid)) - colMeans(xtrain%*%coef(ridge.glm)[-1,])
plot(log(grid),theta0)
points(log(grid),theta0,lwd=2,col="red")

#si on centre ytrain seulement
ridge.glm=glmnet(xtrain,scale(ytrain,scale=FALSE),alpha=0,lambda=grid)
plot(log(grid),coef(ridge.glm)[1,],xlab="log(k)",ylab="intercept")
points(log(grid),coef(ridge.glm)[1,],col="red",pch=20)
abline(h=0,lwd=2,col="blue")
text(-20, 2, "mean(y)", col = "blue")

#si on centre xtrain seulement
ridge.glm=glmnet(scale(xtrain,scale=FALSE),ytrain,alpha=0,lambda=grid)
plot(log(grid),coef(ridge.glm)[1,],xlab="log(k)",ylab="intercept")
points(log(grid),coef(ridge.glm)[1,],col="red",pch=20)
abline(h=mean(ytrain),lwd=2,col="blue")
text(-20, mean(ytrain)+0.02, "mean(y)", col = "blue")

#si on centre tous les deux
ridge.glm=glmnet(scale(xtrain,scale=FALSE),scale(ytrain,scale=FALSE),alpha=0,lambda=grid)
plot(log(grid),coef(ridge.glm)[1,],xlab="log(k)",ylab="intercept")
points(log(grid),coef(ridge.glm)[1,],col="red",pch=20)
abline(h=0,lwd=2,col="blue")
text(-20, 0.02, "mean(y)", col = "blue")

#estimation de theta quand lambda tend vers 0
xscaled <- scale(xtrain)
yscaled <- scale(ytrain)

C <- eigen(t(xscaled)%*%xscaled)
ind <- which(C$values > 0.0001) #length(ind)=rang de xscaled
eig.values <- C$values[ind]
v <- C$vectors[ind,]
u <- eigen(xscaled%*%t(xscaled))$vectors[ind,]

mat <- matrix(0,ncol(xtrain),nrow(xtrain))
for (i in 1:length(ind)){
  m <- matrix(0,ncol(xtrain),nrow(xtrain))
  for (j in 1:nrow(xtrain)){
    m[,j] <- t(v[i,])*u[i,j]
  }
  m <- m*(1/sqrt(C$values[ind[i]]))
  mat <- mat+m
}
theta_l <- mat%*%yscaled
plot(theta_l)
points(theta_l,col="red",pch=20)

#Part3.Question 2
library(MASS)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)

ridge.lm <- lm.ridge(ytrain~xtrain,lambda=grid)
plot(ridge.lm)
plot(lm.ridge(ytrain~xtrain,lambda=grid[40:80]))

all(t(coef(ridge.lm))==coef(ridge.glm))#FALSE, pas les mêmes coéfficients

#calcule de theta tilde
theta.tilde <- matrix(0,701,length(grid))
xcentre <- scale(xtrain,scale=FALSE)
ycentre <- scale(ytrain,scale=FALSE)
for (i in 1:length(grid)){
  theta.tilde[-1,i] <- solve(t(xcentre)%*%xcentre + grid[i]*diag(1,700,700))%*%(t(xcentre)%*%ycentre)
}
#calcule de theta
xbar <- matrix(apply(xtrain,2,mean),nrow(xtrain),ncol(xtrain),byrow=T)
ybar <- rep(mean(ytrain),nrow(xtrain))
theta <- matrix(0,701,length(grid))
for (i in 1:length(grid)){
  theta[-1,i] <- solve(t(xtrain)%*%xtrain + grid[i]*diag(1,700,700))%*%((t(xcentre)%*%xcentre+grid[i]*diag(1,700,700))%*%theta.tilde[-1,i]+t(xbar)%*%ycentre+t(xtrain)%*%ybar)
}
#calcul de l'intercept
for (i in 1:length(grid)){
  theta[1,i] <- mean(ytrain-xtrain%*%theta[-1,i])
}
plot(log(grid),theta[1,],ylab="intercept",xlab="log(k)")
points(log(grid),theta[1,],col="red",pch=20)

#Part3.Question 3
library(pls)
set.seed(87)
cv.seg <- cvsegments(nrow(xtrain),4,type="random")
cv.errors = vector()
grid=10^seq(4,-5.5,length=100)

cv.models  = matrix(0,length(grid),2)

for (i in 1:length(grid)){
  k = grid[i]
  errorsk <- rep(0,length(cv.seg))
  for (j in 1:length(cv.seg)){
    xval <- xtrain[unlist(cv.seg[j]),]
    yval <- ytrain[unlist(cv.seg[j])]
    xt <- xtrain[unlist(cv.seg[-j]),]
    yt <- ytrain[unlist(cv.seg[-j])]
    reg <- glmnet(xt,yt,alpha=0,lambda=k)
    ypred <- predict(reg,xval)
    err <- mean((yval-ypred)^2)
    errorsk[j] = err
  }
  cv.models[i,1] = k
  cv.models[i,2] = mean(errorsk)
}
#kappa associé à l'erreur moyenne minimale
ind <- which(cv.models[,2]==min(cv.models[,2]))
kappa <- cv.models[ind,1]
print(kappa)
#plot EQM
par(mfrow=c(1,1))
plot(log(cv.models[,1]),cv.models[,2],xlab="log(k)",ylab="EQM")
points(log(cv.models[,1]),cv.models[,2],col="red",pch=20)
abline(h=cv.models[ind,2],col="green")
text(-10, cv.models[ind,2]+0.5,round(cv.models[ind,2],4) , col = "green")

set.seed(87)
cv.out=cv.glmnet(xtrain,ytrain,alpha=0,nfolds = 4,lambda=grid) 

plot(cv.out)
abline(h=min(cv.out$cvm),col="green")
text(-10, min(cv.out$cvm)+0.5,round(min(cv.out$cvm),4) , col = "green")

print(cv.out$lambda.min)
# différent de kappa que l'on a trouvé

# on va utiliser le paramètre kappa puisque il a la valeur de l'EQM la plus petite  
reg.ridge <- glmnet(xtrain,ytrain,alpha=0,lambda=kappa)
pred <- predict(reg.ridge,xtest)
err.gener <- mean((ytest-pred)^2)#0.4903069
print(err.gener)
#le kappa que nous avons trouvé donne une erreur de généralisation grande
reg.ridge <- glmnet(xtrain,ytrain,alpha=0,lambda=cv.out$lambda.min)
pred <- predict(reg.ridge,xtest)
err.gener <- mean((ytest-pred)^2)#0.6997558
print(err.gener)
#=> le kappa que nous avons trouvé est mieux


#Partie 4: Régression logistique pénalisée
#Question 1
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
z=sapply(ytrain,function(x){if(x<18){x<-0}
  else{x<-1}})

ztest=sapply(ytest,function(x){if(x<18){x<-0}
  else{x<-1}})
sum(z[which(z==1)])
sum(ztest[which(ztest==1)])
#Question 2
logridge<-cv.glmnet(as.matrix(xtrain),z,alpha=0,lambda=grid,family="binomial")
loglasso<-cv.glmnet(as.matrix(xtrain),z,alpha=1,lambda=grid,family="binomial")

plot(logridge,ylab="Binomial Deviance-Ridge Regression")
plot(loglasso,ylab="Binomial Deviance-Lasso Regression")

#Question 3
ridgemodel=glmnet(as.matrix(xtrain),z,alpha=0,family = "binomial",lambda = logridge$lambda.min)
x.test <- model.matrix(sucres ~., cookie.val)[,-1]
library(magrittr)
probabilities <- ridgemodel %>% predict(newx = x.test)
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
observed.classes <- ztest
1-mean(predicted.classes == observed.classes)
library(pROC)
res.roc <- roc(observed.classes, probabilities)
plot(1-res.roc$specificities,res.roc$sensitivities,type="l",main="ROC régression ridge",xlab="1-specificité"
     ,ylab="sensitivité")
lines(c(0,1),c(0,1),col=2,lty=2)  
segments(c(0,0),c(0,1),c(0,1),c(1,1),lty=3,lwd=2,col=2) 


lassomodel=glmnet(as.matrix(xtrain),z,alpha=1,family = "binomial",lambda = loglasso$lambda.min)
x.test <- model.matrix(sucres ~., cookie.val)[,-1]
probabilities <- lassomodel %>% predict(newx = x.test)
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
observed.classes <- ztest
1-mean(predicted.classes == observed.classes)
res.roc <- roc(observed.classes, probabilities)
plot(1-res.roc$specificities,res.roc$sensitivities,type="l",main="ROC régression lasso",xlab="1-specificité"
     ,ylab="sensitivité")
lines(c(0,1),c(0,1),col=2,lty=2)  
segments(c(0,0),c(0,1),c(0,1),c(1,1),lty=3,lwd=2,col=2)