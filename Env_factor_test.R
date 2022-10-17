# load library
library(rrBLUP)
library(lme4qtl)

# read genotype_data
geno <- read.csv(paste0("../data/N12_genotype.csv"))
# read phenotype and environmental factor data
pheno_and_env <- read.csv("../data/Pheno_and_env_2018.csv")

# make Additive relationship matrix using A.mat function in rrBLUP
colnames(geno)[4] <- "N12_000"
colnames(geno)[5] <- "N12_999"
geno <- cbind(geno[, 1:3], geno[,4:dim(geno)[2]]-1)
M <- t(geno[, -c(1:3)])
K <- A.mat(M,shrink=FALSE)

# check list of environmental factors
env_2018 <- c('mean_temp_0_30','highest_temp_0_30','lowest_temp_0_30','precipitation_0_30','sunshine_duration_0_30','mean_temp_15_45','highest_temp_15_45','lowest_temp_15_45','precipitation_15_45','sunshine_duration_15_45','mean_temp_30_60','highest_temp_30_60','lowest_temp_30_60','precipitation_30_60','sunshine_duration_30_60')

# test statistical significance of environmetnal factors in GxE effect
result_2018 <- data.frame(env=env_2018)
for (trait in c("HD_2018", "GN_2018", "LL_2018")) {
    ps <- c()
    for (env in env_2018) {
        data_2018 <- within(data_2018, tmp_GxE <- data_2018[,env]*data_2018$genotype_chr03_1235072)
        # make models
        model1_formula <- as.formula(paste0(trait, " ~ genotype_chr03_1235072 + ", env, " + (1|id) + (1|Location)"))
        model1 <- relmatLmer(model1_formula, data_2018, relmat = list(id = K), REML=FALSE)
        model3_formula <- as.formula(paste0(trait, " ~ genotype_chr03_1235072*", env, " + (1|id) + (1|Location)"))
        model3 <- relmatLmer(model3_formula, data_2018, relmat = list(id = K), REML=FALSE)

        p <- anova(model1, model3)[2, 8]
        ps <- c(ps, p)
    }
    result_2018 <- cbind(result_2018, ps)
}
colnames(result_2018) <- c("env", "HD_2018", "GN_2018", "LL_2018")
result_2018
# write.table(result_2018, "../data/Env_factor_test_result_2018.csv", sep=",", row.names=FALSE)