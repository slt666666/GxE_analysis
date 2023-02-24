# For variance components of G, E, GxE random effect
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
pheno <- pheno[,c("id", "Location", "GN_2018")]
    
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


# ---------------------------------------------------------------------------------------- #

# For variance components of G, E, SNPxE random effect
# make Additive relationship matrix using A.mat function in rrBLUP
geno <- read.csv(paste0("../data/genotype/", RIL, "_genotype.csv"))
colnames(geno)[4] <- paste0(RIL, "_000")
colnames(geno)[5] <- paste0(RIL, "_999")
geno <- cbind(geno[, 1:3], geno[,4:dim(geno)[2]]-1)
M <- t(geno[, -c(1:3)])
K <- A.mat(M,shrink=FALSE)

pheno <- read.csv("../data/all_data_2year_R.csv")
pheno <- pheno[pheno$RIL=="RIL1",]
pheno <- pheno[,c("id", "Location", "GN_2018")]

# make GxE column
SNP_geno <- t(geno[marker, -c(1:3)])
colnames(SNP_geno) <- "SNP"
SNP_geno <- cbind(id = rownames(SNP_geno), SNP_geno)
tmp_pheno <- merge(pheno, SNP_geno, by="id")

GxE <- c()
for (j in 1:dim(tmp_pheno)[1]) {
    if (is.na(tmp_pheno[j, 4])) {
        GxE <- c(GxE, NA)
    } else if ((tmp_pheno[j, 2] == "Aomori") & (tmp_pheno[j, 4] == "-1")) {
        GxE <- c(GxE, "a")
    } else if ((tmp_pheno[j, 2] == "Iwate") & (tmp_pheno[j, 4] == "-1")) {
        GxE <- c(GxE, "b")
    } else if ((tmp_pheno[j, 2] == "Fukushima") & (tmp_pheno[j, 4] == "-1")) {
        GxE <- c(GxE, "c")
    } else if ((tmp_pheno[j, 2] == "Aomori") & (tmp_pheno[j, 4] == "1")) {
        GxE <- c(GxE, "d")
    } else if ((tmp_pheno[j, 2] == "Iwate") & (tmp_pheno[j, 4] == "1")) {
        GxE <- c(GxE, "e")
    } else if ((tmp_pheno[j, 2] == "Fukushima") & (tmp_pheno[j, 4] == "1")) {
        GxE <- c(GxE, "f")
    } else {
        GxE <- c(GxE, NA)
    } 
}
GxE <- data.frame(GxE)
names(GxE) <- c("GxE")
tmp_pheno <- cbind(tmp_pheno, GxE)
tmp_pheno <-  na.omit(tmp_pheno)

# generate matrix for the effects of locations
ZE <- model.matrix(~factor(tmp_pheno$Location)-1)
ZE2 <- A.mat(ZE,shrink=FALSE)

# generate genomic relationship matrix
K2 <- K[tmp_pheno$id, tmp_pheno$id]

# generate matrix for the effects of SNPxEnv
ME <- model.matrix(~factor(tmp_pheno$GxE)-1)
ME2 <- A.mat(ME,shrink=FALSE)

ETA<-list( SNP=list(~factor(SNP),data=tmp_pheno,model='FIXED'),
            PED=list(K=K2,model='RKHS'),
            ENV=list(K=ZE2,model='RKHS'),
            MxE=list(K=ME2,model='RKHS') )
set.seed(0)
# estimate variance components using BGLR
fm<-BGLR(y=tmp_pheno[, trait],ETA=ETA,saveAt='M2_',nIter=6000,burnIn=2000,verbose=F)

# variance of genotype effect
Vg <- fm$ETA$PED$varU
# variance of enviromental effect
VE <- fm$ETA$ENV$varU
# variance of SNP genotype x E effect
VgE <- fm$ETA$MxE$varU
# variance of residual effect
Ve <- fm$varE
