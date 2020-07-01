# Liwei Wang

# This code is exactly for my undergraduate thesis (major in statistics), the file name is my thesis's title.
# Unfortaintly the thesis was wrote in Chinese, but I will attach that and you could just refer to the 
# English abstract (copied below) and its tables and graphs 

#######################################################

# Abstract

# From the perspective of classification problems and variable selection, this thesis attempts
# to identify key mutation genes that cause two groups of single-cell lesions. Based on the expression
# data of unbalanced, high-dimensional, zero-expansion with small sample single-cell
# RNA gene sequencing data, the following steps were taken for analysis. First, this thesis adopts
# the method of normal distribution test to pre-screen out genetic variables with excessive zero
# expansion and no biological significance. Further, this thesis uses the shrinkage and Logistics
# binary classification method, i.e, Lasso and its variants, to select key genetic variables. At the
# same time, this thesis chooses to use the synthetic sampling method, SMOTE, to generate new
# training data to improve the balance of previous unbalanced data, and effectively improve the
# AUC value of the classifiers and the variable selection results. After finding potential key genetic
# variables that significantly affect cytopathic changes, the candidate gene variables selected in
# this thesis were subjected to T-distribution mean hypothesis test and multiple hypothesis tests
# for secondary verification to obtain the final key gene results. Finally, this article refers to the
# background knowledge of biology and genes online, and explains the final results.

#Keywords: Unbalanced Data; SMOTE; Variable Selection on High-Dimensional Data; Lasso;
#Multiple Hypothesis Testing

########################################################

rm(list=ls())

path <- getwd()
rawdata<-readRDS(paste(path, "/immature_neutrophil.rds", sep=""))

condition<-rawdata[['condition']]
condition<-data.frame(condition)
data<-data.frame(rawdata[['data']])

data<-data.frame(t(data))
condition<-condition$condition
data<-cbind(condition,data)

########################################################

Sum<-apply(data[,-1], 1, sum)
Mean<-apply(data[,-1], 1, mean)
Std<-apply(data[,-1], 1, sd)
Median<-apply(data[,-1], 1, median)
datadiscribe<-cbind(data,Sum,Mean,Std,Median)
head(datadiscribe[,(ncol(datadiscribe)-3):ncol(datadiscribe)])

sum(datadiscribe$condition=='sham')
sum(datadiscribe$condition=='CLP')
plot(Mean,Std,pch=ifelse(factor(condition)=="sham",2, 1),col=
       'black',lty=1,lwd=ifelse(factor(condition)=="sham",3, 1),
     ylab='Standard Deviation', main='Mean-SD Scatter plot')

##############

# Check No missing value

hang<-which(rowSums(is.na(data)) > 0)
hang

###########################################################
library(glmnet)
library(ncvreg)
library(MASS)
library(parallel)
library(caret)
library(pROC)
library(kernlab)
library(ROCR)

################
# Data preprocessing

set.seed(1)
Y<-data$condition
Ytestnum<-c(sample(which(Y=='sham'),length(which(Y=='sham'))/2,
                   replace = F),sample(which(Y=='CLP'),length(which(Y=='CLP'))/2,
                                       replace = F))
Ytestnum<-sample(Ytestnum,length(Ytestnum),replace = F)
Ytest<-Y[Ytestnum]                                      # 1/3 for test data
Ytrain<-Y[-Ytestnum]
Y01<-rep(1,length(Y))
Y01[Y=='sham']<-0
Ytrain01<-rep(1,length(Ytrain))
Ytrain01[Ytrain=='sham']<-0                             #sham (normal data), set to 0，minority part!!!
Ytest01<-rep(1,length(Ytest))
Ytest01[Ytest=='sham']<-0

XXX<-as.matrix(data[,2:ncol(data)])
# Remove the columns that only contains 0
XX<-XXX[,-which(apply(abs(XXX),2,sum)==0)]

############################## S-W normality testing

# First focus on normality testing，P > 0.05 meets the hypothesis that it is normally distribution

P1<-rep(0,ncol(XX))
for (i in 1:ncol(XX)) {
  P1[i]<-shapiro.test(XX[,i])$p.value
}
SigP1<-sort(P1,decreasing = T)[1:635]
SigPorder1<-order(P1,decreasing = T)[1:635]
hist(XX[,SigPorder1[1]],breaks = 20)
hist(XX[,SigPorder1[635]],breaks = 20)

#####  Compared with K-S testing

P2<-rep(0,ncol(XX))
for (i in 1:ncol(XX)) {
  P2[i]<-ks.test(XX[,i],"pnorm",mean(XX[,i]), sd(XX[,i]))$p.value
}
SigP2<-sort(P2,decreasing = T)[1:635]
SigPorder2<-order(P2,decreasing = T)[1:635]
hist(XX[,SigPorder2[1]],breaks = 20)
hist(XX[,SigPorder2[635]],breaks = 20)

########################################################

# Only leave n genes attributes that meets normality and remove another
# Lasso now! (Compared with different kinds of shrinkage methods)

n<-634
SigPorder3<-order(P2,decreasing = T)[1:n]
X<-XX[,SigPorder3]
Xscale<-scale(X)
Xtest<-X[Ytestnum,]
Xtrain<-X[-Ytestnum,]
Xscale<-scale(X)
Xtestscale<-Xscale[Ytestnum,]
Xtrainscale<-Xscale[-Ytestnum,]
dim(Xtrainscale)
dim(Xtestscale)

######################### Phase I: Unbalanced Data Modeling and Results  ########################

#[1]MCP ########################

set.seed(1)
scad.mod<-cv.ncvreg(Xtrainscale,Ytrain01,maxit=1000000000,
                    family='binomial',penalty='MCP')
#plot(scad.mod)
coeffscad<-coef(scad.mod)
sum(coeffscad[-1]!=0)
pred.scad<-predict(scad.mod,Xtestscale,type='response',
                   family='binomial')
pred.scad2<-prediction(pred.scad,Ytest01)
perfscad<-performance(pred.scad2,'sens','spec')
performance(pred.scad2,'auc')@y.values 
#0.8145384
auc.scad<-roc(Ytest01,pred.scad)
print(auc.scad)
plot(auc.scad,ylim=c(0,1),print.thres=T,main=paste
     ('AUC.SCAD',round(auc.scad$auc[[1]],3)))


#[2] adaptiveLasso ####################

set.seed(1)
ridge.mod<-cv.glmnet(Xtrainscale,Ytrain01,alpha=0,
                     family='binomial',type.measure = 'auc')
coeffridge<-coef(ridge.mod)
w<-1/abs(coeffridge[-1])
set.seed(1)
adalasso.mod<-cv.glmnet(Xtrainscale,Ytrain01,alpha=1,
                        penalty.factor=w,family='binomial',type.measure = 'auc')
coeffada2<-coef(adalasso.mod)
length(w)
sum(coeffada2[-1]!=0)
pred.ada<-predict(adalasso.mod,Xtestscale,type='response',
                  family='binomial')
pred.ada2<-prediction(pred.ada,Ytest01)
perfada<-performance(pred.ada2,'sens','spec')
performance(pred.ada2,'auc')@y.values
auc.ada<-roc(Ytest01,pred.ada)
print(auc.ada)
plot(auc.ada,ylim=c(0,1),print.thres=T,main=paste('AUC.ADA',
                                                  round(auc.ada$auc[[1]],3)))


#[3]elastic net ##############################

set.seed(1)
elastic.mod<-cv.glmnet(Xtrainscale,Ytrain01,maxit=1000000,
                       alpha=0.6,family='binomial',type.measure = 'auc')
coeffela<-coef(elastic.mod)
sum(coeffela[-1]!=0)
pred.elas<-predict(elastic.mod,Xtestscale,type='response'
                   ,family='binomial')
pred.elas2<-prediction(pred.elas,Ytest01)
perfelas<-performance(pred.elas2,'sens','spec')
performance(pred.elas2,'auc')@y.values
auc.elas<-roc(Ytest01,pred.elas)
print(auc.elas)
plot(auc.elas,ylim=c(0,1),print.thres=T,main=paste
     ('AUC.elastic',round(auc.elas$auc[[1]],3)))

for (i in seq(0,1,0.1)) {
  set.seed(1)
  elastic.mod<-cv.glmnet(Xtrainscale,Ytrain01,maxit=1000000
                         ,alpha=i,family='binomial',type.measure = 'auc')
  coeffela<-coef(elastic.mod)
  print(i)
  print(sum(coeffela[-1]!=0))
  pred.elas<-predict(elastic.mod,Xtestscale,type='response'
                     ,family='binomial')
  pred.elas2<-prediction(pred.elas,Ytest01)
  perfelas<-performance(pred.elas2,'sens','spec')
  print(performance(pred.elas2,'auc')@y.values)
  print('--------------------')
}


#[4]AdaptiveLasso + elastic net ##############################

set.seed(1)
elastic.mod2<-cv.glmnet(Xtrainscale,Ytrain01,maxit=1000000,
                        alpha=0.8,penalty.factor=w,family='binomial',type.measure = 'auc')
coeffela2<-coef(elastic.mod2)
sum(coeffela2[-1]!=0)
pred.elas<-predict(elastic.mod2,Xtestscale,type='response'
                   ,family='binomial')
pred.elas2<-prediction(pred.elas,Ytest01)
perfelas<-performance(pred.elas2,'sens','spec')
performance(pred.elas2,'auc')@y.values
auc.elas<-roc(Ytest01,pred.elas)
print(auc.elas)
plot(auc.elas,ylim=c(0,1),print.thres=T,main=paste('AUC.elastic
                                                   ',round(auc.elas$auc[[1]],3)))

for (i in seq(0,1,0.1)) {
  set.seed(1)
  elastic.mod<-cv.glmnet(Xtrainscale,Ytrain01,maxit=1000000,
                         alpha=i,penalty.factor=w,family='binomial',type.measure = 'auc')
  coeffela<-coef(elastic.mod)
  print(i)
  print(sum(coeffela[-1]!=0))
  pred.elas<-predict(elastic.mod,Xtestscale,type='response'
                     ,family='binomial')
  pred.elas2<-prediction(pred.elas,Ytest01)
  perfelas<-performance(pred.elas2,'sens','spec')
  print(performance(pred.elas2,'auc')@y.values)
  print('--------------------')
}


#############################

Ytrain01<-as.factor(Ytrain01)
Ytest01<-as.factor(Ytest01)
train<-cbind(Ytrain01,data.frame(Xtrain))
test<-cbind(Ytest01,data.frame(Xtest))

# Data Visualization Using PCA and TSNE #############################

# TSNE

library(Rtsne)
library(rgl)
library(scatterplot3d)

color<-rep('black', length(Y))
color[which(Y=='sham')]<-'green'      # Green represents 'sham' data
set.seed(1)
tsne2<-Rtsne(X,dims = 2)
plot(tsne2$Y[,1],tsne2$Y[,2],main='TSNE',col=color,pch=20)
set.seed(1)
tsne3<-Rtsne(X,dims = 3,perplexity = 30,verbose = T,
             max_iter = 1000)
scatterplot3d(tsne3$Y[,1],tsne3$Y[,2],tsne3$Y[,3],color=color,
              highlight.3d=F,pch = 20,angle=80,main='TSNE',
              type = 'p',grid = F,col.grid = "lightblue",
              col.axis = 'lightblue',col.lab = T,
              scale.y =2,xlim = c(-20,20),ylim = c(-25,25),
              zlim = c(-30,30))
# PCA
df.pr <- princomp(t(X))
plot(df.pr$loadings[,1],df.pr$loadings[,2],main='PCA',
     col=color,pch=20,xlab = 'Comp.1',ylab = 'Comp.2')
scatterplot3d(df.pr$loadings[,1],df.pr$loadings[,2],
              df.pr$loadings[,3],color=color,
              highlight.3d=F,pch = 20,angle=80,main='PCA',
              type = 'p',grid = F,col.grid = "lightblue",
              col.axis = 'lightblue',col.lab = T,
              scale.y =5,xlab = 'Comp.1',ylab = 'Comp.2',zlab = 'Comp.3')

######################### Phase II: Balanced Data Modeling with SMOTE and Results.  ########################

library(DMwR)

#SMOTE
#now using SMOTE to create a more "balanced problem"

set.seed(1)
dataSMOTE <- SMOTE(Ytrain01 ~ ., train,perc.over = 2500
                   ,perc.under = 105,k=5)
table(dataSMOTE$Ytrain01)
dim(dataSMOTE)
YSMOTE<-dataSMOTE$Ytrain01
Y01.SMOTE<-rep(0,length(YSMOTE))
Y01.SMOTE[YSMOTE=='1']<-1                # Same, 0 represents 'sham'
X.SMOTE<-dataSMOTE[,2:ncol(dataSMOTE)]
Xscale.SMOTE<-scale(X.SMOTE)
Ytrain01.SMOTE<-Y01.SMOTE
Xtrain.SMOTE<-X.SMOTE
Xtrainscale.SMOTE<-Xscale.SMOTE
dim(Xtrainscale.SMOTE)


# Data Visualization Using PCA and TSNE #############################

# TSNE

color.SMOTE<-rep('black', length(YSMOTE))
color.SMOTE[which(YSMOTE=='0')]<-'green'  
set.seed(1)
tsne2<-Rtsne(X,dims = 2)
plot(tsne2$Y[,1],tsne2$Y[,2],main='TSNE',col=color.SMOTE,pch=20)
set.seed(1)
tsne.SMOTE<-Rtsne(X.SMOTE,dims = 3,perplexity = 30,
                  check_duplicates = FALSE,verbose=T,max_iter=1000)
scatterplot3d(tsne.SMOTE$Y[,1],tsne.SMOTE$Y[,2],
              tsne.SMOTE$Y[,3],color=color.SMOTE,
              highlight.3d=F,pch = 20,angle=80,main='TSNE',
              type = 'p',grid = F,col.grid = "lightblue",
              col.axis = 'lightblue',col.lab = T,
              scale.y =5,xlim = c(-20,20),ylim = c(-25,25)
              ,zlim = c(-30,30))
# PCA
set.seed(1)
df.pr <- princomp(t(X.SMOTE))
plot(df.pr$loadings[,1],df.pr$loadings[,2],main='PCA',
     col=color.SMOTE,pch=20,xlab = 'Comp.1',ylab = 'Comp.2')
scatterplot3d(df.pr$loadings[,1],df.pr$loadings[,2],
              df.pr$loadings[,3],color=color.SMOTE,
              highlight.3d=F,pch = 20,angle=80,main='PCA',
              type = 'p',grid = F,col.grid = "lightblue",
              col.axis = 'lightblue',col.lab = T,
              scale.y =5,xlab = 'Comp.1',ylab = 'Comp.2',zlab = 'Comp.3')


# Data Modeling after SMOTE

#1.MCP ###########################################

set.seed(1)
scad.mod.SMOTE<-cv.ncvreg(Xtrainscale.SMOTE,
                          Ytrain01.SMOTE,family='binomial',penalty='SCAD')
coeffscad.SMOTE<-coef(scad.mod.SMOTE)
sum(coeffscad.SMOTE[-1]!=0)
pred.scad.SMOTE<-predict(scad.mod.SMOTE,Xtestscale,
                         type='response',family='binomial')
pred.scad.SMOTE2<-prediction(pred.scad.SMOTE,Ytest01)
perfscad.SMOTE<-performance(pred.scad.SMOTE2,'sens','spec')
performance(pred.scad.SMOTE2,'auc')@y.values
auc.scad.SMOTE<-roc(Ytest01,pred.scad.SMOTE)
print(auc.scad.SMOTE)
plot(auc.scad.SMOTE,ylim=c(0,1),print.thres=T,main=paste
     ('AUC.SCAD.SMOTE',round(auc.scad.SMOTE$auc[[1]],3)))


#2.Adaptive ###########################################

set.seed(1)
ridge.mod2.SMOTE<-cv.glmnet(Xtrainscale.SMOTE,Ytrain01.SMOTE,
                            alpha=0,family='binomial',type.measure = 'auc')
w.SMOTE<-1/abs(coef(ridge.mod2.SMOTE)[-1])
set.seed(1)
elastic.mod.SMOTE<-cv.glmnet(Xtrainscale.SMOTE,Ytrain01.SMOTE
                             ,alpha=1,penalty.factor=w.SMOTE,family='binomial',
                             type.measure = 'auc')
coeffela.SMOTE<-coef(elastic.mod.SMOTE)
sum(coeffela.SMOTE[-1]!=0)
pred.elas.SMOTE<-predict(elastic.mod.SMOTE,Xtestscale,
                         type='response',family='binomial')
pred.elas2.SMOTE<-prediction(pred.elas.SMOTE,Ytest01)
perfelas.SMOTE<-performance(pred.elas2.SMOTE,'sens','spec')
performance(pred.elas2.SMOTE,'auc')@y.values
auc.elas.SMOTE<-roc(Ytest01,pred.elas.SMOTE)
print(auc.elas.SMOTE)
plot(auc.elas.SMOTE,ylim=c(0,1),print.thres=T,main=paste
     ('AUC.ELAS.SMOTE',round(auc.elas.SMOTE$auc[[1]],3)))

#3.elastic net ###########################################

set.seed(1)
elastic.mod.SMOTE<-cv.glmnet(Xtrainscale.SMOTE,Ytrain01.SMOTE,alpha=0.8,family='binomial',type.measure = 'auc')
coeffela.SMOTE<-coef(elastic.mod.SMOTE)
sum(coeffela.SMOTE[-1]!=0)
pred.elas.SMOTE<-predict(elastic.mod.SMOTE,Xtestscale,
                         type='response',family='binomial')
pred.elas2.SMOTE<-prediction(pred.elas.SMOTE,Ytest01)
perfelas.SMOTE<-performance(pred.elas2.SMOTE,'sens','spec')
performance(pred.elas2.SMOTE,'auc')@y.values
auc.elas.SMOTE<-roc(Ytest01,pred.elas.SMOTE)
print(auc.elas.SMOTE)
plot(auc.elas.SMOTE,ylim=c(0,1),print.thres=T,main=paste
     ('AUC.ELAS.SMOTE',round(auc.elas.SMOTE$auc[[1]],3)))

for (i in seq(0,1,0.1)) {
  set.seed(1)
  elastic.mod.SMOTE<-cv.glmnet(Xtrainscale.SMOTE,Ytrain01.SMOTE,alpha=i,family='binomial',type.measure = 'auc')
  coeffela.SMOTE<-coef(elastic.mod.SMOTE)
  print(i)
  print(sum(coeffela.SMOTE[-1]!=0))
  pred.elas.SMOTE<-predict(elastic.mod.SMOTE,Xtestscale
                           ,type='response',family='binomial')
  pred.elas.SMOTE2<-prediction(pred.elas.SMOTE,Ytest01)
  perfelas<-performance(pred.elas.SMOTE2,'sens','spec')
  print(performance(pred.elas.SMOTE2,'auc')@y.values)
  print('--------------------')
}

#4.adaptiveLasso + elastic net ###########################################

set.seed(1)
ridge.mod2.SMOTE<-cv.glmnet(Xtrainscale.SMOTE,Ytrain01.SMOTE,alpha=0,family='binomial',type.measure = 'auc')
w.SMOTE<-1/abs(coef(ridge.mod2.SMOTE)[-1])
set.seed(1)
elastic.mod.SMOTE<-cv.glmnet(Xtrainscale.SMOTE,Ytrain01.SMOTE,alpha=0.9,penalty.factor=w.SMOTE,family='binomial',
                             type.measure = 'auc')
coeffela.SMOTE<-coef(elastic.mod.SMOTE)
sum(coeffela.SMOTE[-1]!=0)
pred.elas.SMOTE<-predict(elastic.mod.SMOTE,Xtestscale,type='response',family='binomial')
pred.elas2.SMOTE<-prediction(pred.elas.SMOTE,Ytest01)
perfelas.SMOTE<-performance(pred.elas2.SMOTE,'sens','spec')
performance(pred.elas2.SMOTE,'auc')@y.values#0.873
auc.elas.SMOTE<-roc(Ytest01,pred.elas.SMOTE)
print(auc.elas.SMOTE)
plot(auc.elas.SMOTE,ylim=c(0,1),print.thres=T,main=paste('
                                                         AUC.ELAS.SMOTE',round(auc.elas.SMOTE$auc[[1]],3)))

set.seed(1)
ridge.mod2.SMOTE<-cv.glmnet(Xtrainscale.SMOTE,Ytrain01.SMOTE,
                            alpha=0,family='binomial',type.measure = 'auc')
w.SMOTE<-1/abs(coef(ridge.mod2.SMOTE)[-1])

for (i in seq(0,1,0.1)) {
  set.seed(1)
  elastic.mod.SMOTE<-cv.glmnet(Xtrainscale.SMOTE,Ytrain01.SMOTE
                               ,alpha=i,penalty.factor=w.SMOTE,family='binomial',
                               type.measure = 'auc')
  coeffela.SMOTE<-coef(elastic.mod.SMOTE)
  print(i)
  print(sum(coeffela.SMOTE[-1]!=0))
  pred.elas.SMOTE<-predict(elastic.mod.SMOTE,Xtestscale,
                           type='response',family='binomial')
  pred.elas.SMOTE2<-prediction(pred.elas.SMOTE,Ytest01)
  perfelas<-performance(pred.elas.SMOTE2,'sens','spec')
  print(performance(pred.elas.SMOTE2,'auc')@y.values)
  print('--------------------')
}

######################### Phase III: Multiple Hypothesis Test for potiential genes found by previous models  ########################

# T-Testing：
# Calr，Tspo，Ubc，Naca，Serinc3，Rasgrp4，Hsp90ab1，Sp3，Nt5c3，Sf3b2，Dock2，Lilrb4a

# 1.Calr
Calr0<-X[,'Calr'][which(Y01==0)] # Select sham, normal part
Calr1<-X[,'Calr'][which(Y01==1)] # Select CLP
t.test(Calr0,Calr1,var.equal = TRUE,alternative = "less") # sham's mean siginicantly < CLP by t-test
pCalr<-t.test(Calr0,Calr1,var.equal = TRUE,alternative 
              = "less")$p.value

# 2.Tspo
Tspo0<-X[,'Tspo'][which(Y01==0)]
Tspo1<-X[,'Tspo'][which(Y01==1)]
t.test(Tspo0,Tspo1,var.equal = TRUE,alternative = "less") # sham's mean siginicantly < CLP by t-test
pTspo<-t.test(Tspo0,Tspo1,var.equal = TRUE,alternative
              = "less")$p.value

# 3.Ubc
Ubc0<-X[,'Ubc'][which(Y01==0)]
Ubc1<-X[,'Ubc'][which(Y01==1)]
t.test(Ubc0,Ubc1,var.equal = TRUE,alternative = "greater") # sham's mean siginicantly > CLP by t-test
pUbc<-t.test(Ubc0,Ubc1,var.equal = TRUE,alternative 
             = "greater")$p.value

# 4.Naca
Naca0<-X[,'Naca'][which(Y01==0)]
Naca1<-X[,'Naca'][which(Y01==1)]
t.test(Naca0,Naca1,var.equal = TRUE,alternative = "greater") # sham's mean siginicantly > CLP by t-test
pNaca<-t.test(Naca0,Naca1,var.equal = TRUE,alternative
              = "greater")$p.value

# 5.Serinc3
Serinc30<-X[,'Serinc3'][which(Y01==0)]
Serinc31<-X[,'Serinc3'][which(Y01==1)]
t.test(Serinc30,Serinc31,var.equal = TRUE,alternative = "greater") # sham's mean siginicantly > CLP by t-test
pSerinc3<-t.test(Serinc30,Serinc31,var.equal = TRUE,alternative 
                 = "greater")$p.value

# 6.Rasgrp4
Rasgrp40<-X[,'Rasgrp4'][which(Y01==0)]
Rasgrp41<-X[,'Rasgrp4'][which(Y01==1)]
t.test(Rasgrp40,Rasgrp41,var.equal = TRUE,alternative = "greater") # sham's mean siginicantly > CLP by t-test
pRasgrp4<-t.test(Rasgrp40,Rasgrp41,var.equal = TRUE,alternative
                 = "greater")$p.value

# 7.Hsp90ab1
Hsp90ab10<-X[,'Hsp90ab1'][which(Y01==0)]
Hsp90ab11<-X[,'Hsp90ab1'][which(Y01==1)]
t.test(Hsp90ab10,Hsp90ab11,var.equal = TRUE,alternative = "less") # sham's mean siginicantly < CLP by t-test
pHsp90ab1<-t.test(Hsp90ab10,Hsp90ab11,var.equal = TRUE,alternative
                  = "less")$p.value

# 8.Sp3
Sp30<-X[,'Sp3'][which(Y01==0)]
Sp31<-X[,'Sp3'][which(Y01==1)]
t.test(Sp30,Sp31,var.equal = TRUE,alternative = "two.sided")
t.test(Sp30,Sp31,var.equal = TRUE,alternative = "less")
t.test(Sp30,Sp31,var.equal = TRUE,alternative = "greater") # sham's mean is not siginicantly different from CLP by t-test
pSp3<-t.test(Sp30,Sp31,var.equal = TRUE,alternative = "less")$p.value

# 9.Nt5c3
Nt5c30<-X[,'Nt5c3'][which(Y01==0)]
Nt5c31<-X[,'Nt5c3'][which(Y01==1)]
t.test(Nt5c30,Nt5c31,var.equal = TRUE,alternative = "less") # sham's mean siginicantly < CLP by t-test
pNt5c3<-t.test(Nt5c30,Nt5c31,var.equal = TRUE,alternative
               = "less")$p.value

# 10.Sf3b2
Sf3b20<-X[,'Sf3b2'][which(Y01==0)]
Sf3b21<-X[,'Sf3b2'][which(Y01==1)]
t.test(Sf3b20,Sf3b21,var.equal = TRUE,alternative = "two.sided")
t.test(Sf3b20,Sf3b21,var.equal = TRUE,alternative = "less")
t.test(Sf3b20,Sf3b21,var.equal = TRUE,alternative = "greater") # sham's mean is not siginicantly different from CLP by t-test
pSf3b2<-t.test(Sf3b20,Sf3b21,var.equal = TRUE,alternative
               = "greater")$p.value

# 11.Dock2
Dock20<-X[,'Dock2'][which(Y01==0)]
Dock21<-X[,'Dock2'][which(Y01==1)]
t.test(Dock20,Dock21,var.equal = TRUE,alternative = "two.sided")
t.test(Dock20,Dock21,var.equal = TRUE,alternative = "less")
t.test(Dock20,Dock21,var.equal = TRUE,alternative = "greater") # sham's mean is not siginicantly different from CLP by t-test
pDock2<-t.test(Dock20,Dock21,var.equal = TRUE,alternative
               = "greater")$p.value

# 12.Lilrb4a
Lilrb4a0<-X[,'Lilrb4a'][which(Y01==0)]
Lilrb4a1<-X[,'Lilrb4a'][which(Y01==1)]
t.test(Lilrb4a0,Lilrb4a1,var.equal = TRUE,alternative = "less") # sham's mean siginicantly < CLP by t-test
pLilrb4a<-t.test(Lilrb4a0,Lilrb4a1,var.equal = TRUE,alternative
                 = "less")$p.valu


# Multiple Hypothesis Test -- FDR

Pless<-cbind(pCalr,pTspo,pHsp90ab1,pNt5c3,pLilrb4a,pSp3)
Pgreater<-c(pUbc,pNaca,pSerinc3,pRasgrp4,pSf3b2,pDock2)
PFDRless<-p.adjust(Pless,method = 'fdr',n=length(Pless))
PFDRgreater<-p.adjust(Pgreater,method = 'fdr',n=length(Pgreater))

