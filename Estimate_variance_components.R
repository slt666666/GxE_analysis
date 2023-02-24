# load library
library(rrBLUP)
library(BGLR)

# read genotype_data
geno <- read.csv(paste0("../data/N12_genotype.csv"))
# read phenotype_data
pheno <- read.csv("../data/all_data_2year_for_R.csv")

# make Additive relationship matrix using A.mat function in rrBLUP
colnames(geno)[4] <- paste0(RIL, "_000")
colnames(geno)[5] <- paste0(RIL, "_999")
geno <- cbind(geno[, 1:3], geno[,4:dim(geno)[2]]-1)
M <- t(geno[, -c(1:3)])
K <- A.mat(M,shrink=FALSE)
    
pheno <- pheno[pheno$RIL=="RIL1",]
pheno <- pheno[,c("id", "Location", trait)]
    
# generate matrix for the effects of locations
ZE <- model.matrix(~factor(pheno$Location)-1)
ZE2 <- A.mat(ZE,shrink=FALSE)

# generate genomic relationship matrix
K2 <- K[pheno$id, pheno$id]

# generate incidence matrix for GxE effect
# followed "Paulino et al., 2015, Crop Science"
ZVAR <- model.matrix(~factor(x=pheno$id, levels=rownames(K), ordered=TRUE)-1)
K3 <- K
diag(K3) <- diag(K3)+1e-4
L <- t(chol(K3/mean(diag(K3))))
ZAZ <- tcrossprod(ZVAR%*%L)
ZZ <- tcrossprod(ZE)
K4 <- ZZ*ZAZ
diag(K4) <- diag(K4)+1e-4
L3 <- t(chol(K4/mean(diag(K4))))

# estimate variance components using BGLR
ETA<-list( PED=list(K=K2,model='RKHS'),
            ENV=list(K=ZE2,model='RKHS'),
            GxE=list(X=L3,model='BRR') )
set.seed(0)
fm<-BGLR(y=pheno[, trait],ETA=ETA,saveAt='M2_',nIter=6000,burnIn=2000,verbose=F)

# variance of genotype effect
Vg <- fm$ETA$PED$varU
# variance of environmental effect
VE <- fm$ETA$ENV$varU
# variance of GxE effect
VgE <- fm$ETA$GxE$varB
# variance of residual effect
Ve <- fm$varE
