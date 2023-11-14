###################################################
# Title: MaMD Analytical: A Computer Program for 
# Macromorphoscopic Trait Analysis in Forensic Anthropology
# Purpose: Publication and validation of MaMD Analytical v. 1.0.0
# Version: 1.0.0
# Author: jth, rrd
# Last updated: 17 Mar 2022(jth)
# Updates: 
# 1. 17 Mar 2022: jth initiated code
# 2.
# Notes:
#
#
###################################################
# START WITH A CLEAN WORKING DIRECTORY
rm(list=ls())
###################################################
# Install libraries and dependencies
###################################################
# par
set.seed(1234)
digits=4
options(scipen = 999)
###################################################
library(ModelMetrics)
library(nnet)
library(dplyr)
library(knitr)
library(reshape2)
library(caret)
library(stats)
library(ggplot2)
library(tidyverse)
###################################################
# install data
data<-read.csv("mamd.csv", sep=',', header = T)
data<-as.data.frame(data)
data.1<-read.csv("mamd.csv", sep=',', header = T)
head(data)
#test<-read.csv("test.csv", sep=',', header = T )
mamd<-data
mamd.1<-data.1
###################################################
head(mamd)
names(mamd)<-c("ancestry", 'ANS', 'INA', 'IOB', 'MT', 'NAW', 'NBC', 'NO', 'PBD', 'PZT', 'ZS')
names(mamd.1)<-c("ancestry", 'ANS', 'INA', 'IOB', 'MT', 'NAW', 'NBC', 'NO', 'PBD', 'PZT', 'ZS')
str(mamd)
mamd$ancestry<-as.factor(mamd$ancestry)
mamd.1$ancestry<-as.factor(mamd.1$ancestry)
mamd$ANS<-as.factor(mamd$ANS)
mamd$INA<-as.factor(mamd$INA)
mamd$IOB<-as.factor(mamd$IOB)
mamd$MT<-as.factor(mamd$MT)
mamd$NAW<-as.factor(mamd$NAW)
mamd$NBC<-as.factor(mamd$NBC)
mamd$NO<-as.factor(mamd$NO)
mamd$PBD<-as.factor(mamd$PBD)
mamd$PZT<-as.factor(mamd$PZT)
mamd$ZS<-as.factor(mamd$ZS)
###################################################
mamd<-na.omit(mamd) %>% droplevels()
mamd.1<-na.omit(mamd.1) %>% droplevels()
###################################################
#data selection for tuning and training
## 75% of the sample size
smp_size <- floor(0.75 * nrow(mamd))
train_ind <- sample(seq_len(nrow(mamd)), size = smp_size)

train <- mamd[train_ind, ]
test <- mamd[-train_ind, ]
###################################################

#MaMDnnet<-function(data){}
fit8<-nnet(ancestry~., data = train, size = 15, rang = 0.1, 
          decay = 5e-4, maxit = 20000)

f8<-fitted(fit8)                                                                                           
f.8<-melt(f8)
names(f.8)<-c("Entry","Group","Posterior")
f.8 <- f.8[order(f.8$Posterior),] 
f8<-f.8[1,]
mod8 <- predict(fit8, type="class")

ctab8<-table(train$ancestry, mod8)
ctab.prop8<-round(prop.table(ctab8,1)*100,digits=2)
diag(prop.table(ctab8,1))
sum(diag(prop.table(ctab8)))

#test model on hold-out sample

raw.post3<-predict(fit8, newdata = test,type=c("raw"))
pred.o<-predict(fit8, newdata = test,type=c("class"))
post.prob<-(pred.o)
ctab<-table(test$ancestry,post.prob)
ctab.prop<-round(prop.table(ctab,1)*100,digits=2)
diag(prop.table(ctab,1))
sum(diag(prop.table(ctab)))

# For graphic
nnet.1<-neuralnet(ancestry~., train, linear.output = FALSE,
                  algorithm = "rprop+", hidden=c(15),
                  threshold = 0.1,rep=5, stepmax = 1e+06,
                  learningrate = 0.001, lifesign="full")

plot(nnet.1)

# (list(post.prob, raw.post3))
# 
# MaMDnnet(mamd)
###################################################


###################################################
# end of program
###################################################



