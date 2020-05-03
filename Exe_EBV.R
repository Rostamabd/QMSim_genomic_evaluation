#!/usr/bin/env Rscript

rm(list=ls())
library(BGLR)

######make a mrker file using QMSim marker file########
a <- read.table("./r_toy_par/p1_mrk_001.txt",skip = 1,colClasses = c("numeric","character"))

cat("dimention_geno = ", dim(a), "\n")

#head(a)
ID <- a[,1]

X <- data.frame(do.call(rbind, strsplit(a[,2],"")))
X[] <- as.numeric(factor(as.matrix(X)))
rm(a)

cat("dimention_New_geno = ", dim(X), "\n")

### recode hetrozygots to 1
X[X==3] <- 1
X[X==4] <- 1

### Quality control for MAF

      snpremove<-(colMeans(X)/2)>0.05 & (colMeans(X)/2)<0.95
      X <- X[,snpremove]

##############make Phen#############
Data <- read.table("data.tmp",header=T)

### Just the recent 2000 animals will be picked up for the prediction
if(nrow(Data)> 2000){
X <- X[((nrow(Data)-2000)+1):nrow(Data),]
ID <- ID[((nrow(Data)-2000)+1):nrow(Data)]
Data <- Data[((nrow(Data)-2000)+1):nrow(Data),]
}

cat("dimention_Data = ", dim(Data), "\n")

yNA <- Data[(1:nrow(Data)),11]


### Scale genotype (X) matrix
X <- scale(X,center=TRUE,scale=FALSE)


nIter <- 10000 ; burnIn <- 2000

ETA<-list(list(X=X,model='BayesB'))

fmBB<-BGLR(y=yNA,ETA=ETA, nIter=nIter, burnIn=burnIn,verbose = FALSE,saveAt='BB_')


GEBV <- cbind(ID,fmBB$yHat)

colnames(GEBV) <- c("ID","EBV")


write.table(GEBV,file="my_bv.txt",row.names=F,quote=F)

