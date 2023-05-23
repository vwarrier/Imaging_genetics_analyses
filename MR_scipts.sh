###MR scripts###

###Step 1 - rerun the GWAS by splitting it into two halves
library(data.table)
library(dplyr)
pheno = fread("~/UKB_v2/UKB_WM/WMpheno_scaled.txt")
pheno2 = pheno[,c("FID", "IID")]
alpha = sample_n(pheno2, 15898, replace = FALSE) #15898
beta = pheno2[!pheno2$FID %in% alpha$IID,] #15899

write.table(alpha, file = "~/UKB_v2/UKB_WM/keep_alpha2.txt", row.names = F, col.names = T, quote = F)
write.table(beta, file = "~/UKB_v2/UKB_WM/keep_beta2.txt", row.names = F, col.names = T, quote = F)




# Whole brain Structural T1 and T2 with T1 T2 covariate, R
for i in {1..8}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/Structural_WB.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete_structural.txt --threads 20 --keep ~/UKB_v2/UKB_WM/keep_alpha2.txt --out ./GWAS/Whole_brain/alpha/alpha2_Structural_GWAS_$i; done

for i in {1..8}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/Structural_WB.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete_structural.txt --threads 20 --keep ~/UKB_v2/UKB_WM/keep_beta2.txt --out ./GWAS/Whole_brain/beta/beta2_Structural_GWAS_$i; done



# Whole brain DTI
for i in {1..5}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/WMpheno_scaled.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete.txt --threads 20 --keep ~/UKB_v2/UKB_WM/keep_alpha2.txt --out ./GWAS/Whole_brain/alpha/alpha2_mWM_GWAS$i; done

for i in {1..5}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/WMpheno_scaled.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete.txt --threads 20 --keep ~/UKB_v2/UKB_WM/keep_beta2.txt --out ./GWAS/Whole_brain/beta/beta2_mWM_GWAS$i; done



### Clump, bash (change as needed)

for i in {1..5}; do ./plink \
--noweb \
--bfile fileforclumping \
--clump ~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/alpha/alpha2_mWM_GWAS${i}.fastGWA \
--clump-field P \
--clump-p1 0.00000005 \
--clump-p2 1 \
--clump-r2 0.01 \
--clump-kb 1000 \
--threads 15 \
--out ~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/alpha/Clumped/clumped_alpha2_mWM_GWAS$i; done



###Step 2 MR with twosample MR and MRPResso, R

library(TwoSampleMR)
library(data.table)
library(MRPRESSO)



#### Run two sample MR - forward direction


SNP = fread("~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/beta/Clumped/clumped_beta2_Structural_GWAS_6.clumped")
SNP = SNP[,c("SNP")]

SA = fread("~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/beta/beta2_Structural_GWAS_6.fastGWA")
Exposure = merge(SNP, SA, by = "SNP")

Exposure = subset(Exposure, AF1 > 0.01 & AF1 < 0.99)

setnames(Exposure, old = c('BETA','SE', 'A1', 'A2', "AF1"), new = c('beta','se', 'effect_allele', 'other_allele', "eaf"))


Exposure <- format_data(Exposure, type="exposure")


outcome_dat <- read_outcome_data(
  snps = Exposure$SNP,
  filename = "~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/alpha/alpha2_Structural_GWAS_1.fastGWA",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P"
)

#Change the outcome data as needed


dat = harmonise_data(Exposure, outcome_dat, action = 1)
res = mr(dat, method_list=c("mr_ivw",  "mr_weighted_median", "mr_egger_regression"))
res
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)

dat2 = dat[-c(1),]

res = mr(dat2, method_list=c("mr_ivw",  "mr_weighted_median", "mr_egger_regression"))
res
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat2, NbDistribution = 1000,  SignifThreshold = 0.05)


res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

#--------------------------------
#Steiger analyses
#--------------------------------

dat$samplesize.exposure = 16000
dat$samplesize.outcome =  16000  
directionality_test(dat[1,])




###plots

mr_scatter_plot(res, dat) #mrplot
mr_forest_plot(res_single) #forest plot
mr_leaveoneout_plot(res_loo) # leave one out plot



##------------------------------------
#CAUSE
##------------------------------------

#Step 1: Install packages

devtools::install_github("jean997/cause@v1.2.0")
library(cause)
library(readr)
library(dplyr)
library(data.table)

#Step 2: Read and format data

X1 <- fread("~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/alpha/alpha2_LGI.fastGWA")
X2 <- fread("~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/beta/beta2_Structural_GWAS_6.fastGWA")

X <- gwas_merge(X1, X2, snp_name_cols = c("SNP", "SNP"),
                beta_hat_cols = c("BETA", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("A1", "A1"), 
                A2_cols = c("A2", "A2"))

#Step 3: Calculate nuisance parameters

set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)


##variants -read

SNP1 = fread("~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/alpha/Clumped/clumped_alpha2_LGI_CAUSE.clumped")
SNP2 = fread("~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/alpha/Clumped/clumped_alpha2_LGI.clumped")
top_vars1 = SNP1[,c("SNP")]
top_vars1_2 = top_vars1[top_vars1$SNP %in% X$snp,]

top_vars2= SNP2[,c("SNP")]
top_vars2_2 = top_vars2[top_vars2$SNP %in% X$snp,]



res1 <- cause(X=X, variants = top_vars1_2$SNP, param_ests = params, force = TRUE)
res2 <- cause(X=X, variants = top_vars2_2$SNP, param_ests = params, force = TRUE)

summary(res2, ci_size=0.95)

summary(res1, ci_size=0.95)

