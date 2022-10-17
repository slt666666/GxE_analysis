# load library
library(rrBLUP)
library(lme4qtl)

# read genotype_data
geno <- read.csv(paste0("../data/N12_genotype.csv"))
# read phenotype_data
pheno <- read.csv("../data/all_data_2year_for_R.csv")

# extract RIL1 data & GN_2018 data as a example
pheno <- pheno[pheno$RIL=="RIL1",]
pheno <- pheno[,c("id", "Location", "GN_2018")]

# make Additive relationship matrix using A.mat function in rrBLUP
colnames(geno)[4] <- "N12_000"
colnames(geno)[5] <- "N12_999"
geno <- cbind(geno[, 1:3], geno[,4:dim(geno)[2]]-1)
M <- t(geno[, -c(1:3)])
K <- A.mat(M,shrink=FALSE)

# statical test for GxE effect using lme4qtl package
score.calc <- function(geno, pheno, trait, K) {
  m <- nrow(geno)
  score <- array(0, m)
  for (i in 1:m) {
    SNP_geno <- t(geno[i, -c(1:3)])
    colnames(SNP_geno) <- "SNP"
    SNP_geno <- cbind(id = rownames(SNP_geno), SNP_geno)
    tmp_pheno <- merge(pheno, SNP_geno, by="id")
    tmp_pheno <-  na.omit(tmp_pheno)
    
    model1_formula <- as.formula(paste0(trait, " ~ SNP + (1|id) + (1|Location)"))
    model1 <- relmatLmer(model1_formula, tmp_pheno, relmat = list(id = K), REML=FALSE)
    model2_formula <- as.formula(paste0(trait, " ~ SNP + (1|id) + (1+SNP|Location)"))
    model2 <- relmatLmer(model2_formula, tmp_pheno, relmat = list(id = K), REML=FALSE)
    
    # maximum likelihood test (model1 vs model2)
    score[i] <- anova(model1, model2)[2, 8]
    print(paste(i, score[i]))
  }
  return(score)
}

# statistical test for all SNPs
m <- nrow(geno)
it <- split(1:m,factor(cut(1:m,n.core,labels=FALSE)))
score <- unlist(parallel::mclapply(it, function(SNPs){score.calc(geno[SNPs, ], pheno, trait, K)},mc.cores=n.core))

result <- within(geno[1:m, 1:3], Score <- score)
# write.table(result, "../data/result/N12_GN_2018.csv", sep=",", row.names=FALSE)
