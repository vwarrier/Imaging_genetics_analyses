##Step 1: Generate the weights using PRScs

# 1. In R
list1 = list("SA", "Volume", "CT")

for (i in list1){
  data1 = fread(paste0("~/METAL/Whole_brain_meta/final_GWAS/", i, ".txt"))
  data2 = data1[,c("SNP", "A1", "A2", "BETA", "P")]
  write.table(data2, file = paste0("~/METAL/Whole_brain_meta/final_GWAS/", i, "_forPRSice.txt"), row.names = F, col.names = T, quote = F)}

# 2. in bash
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --sst_file=/mnt/b2/home4/arc/vw260/METAL/Whole_brain_meta/final_GWAS/SA_forPRScs.txt --n_gwas=36843 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_SA  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/files_imputed/SFARImergedall --sst_file=/mnt/b2/home4/arc/vw260/METAL/Whole_brain_meta/final_GWAS/SA_forPRScs.txt --n_gwas=36843 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SSC_SA  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --sst_file=/mnt/b2/home4/arc/vw260/METAL/Whole_brain_meta/final_GWAS/Volume_forPRScs.txt --n_gwas=36843 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_Volume  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/files_imputed/SFARImergedall --sst_file=/mnt/b2/home4/arc/vw260/METAL/Whole_brain_meta/final_GWAS/Volume_forPRScs.txt --n_gwas=36843 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SSC_Volume  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --sst_file=/mnt/b2/home4/arc/vw260/METAL/Whole_brain_meta/final_GWAS/CT_forPRScs.txt --n_gwas=36843 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SPARK_CT  --phi=1e-2 
python PRScs.py --ref_dir=./ldblk_1kg_eur --bim_prefix=/mnt/b2/home4/arc/vw260/SFARI/liftOverPlink/files_imputed/SFARImergedall --sst_file=/mnt/b2/home4/arc/vw260/METAL/Whole_brain_meta/final_GWAS/CT_forPRScs.txt --n_gwas=36843 --out_dir=/mnt/b2/home4/arc/vw260/Autism_heterogeneity/temp_pgs_files/SSC_CT  --phi=1e-2 


##Step 2: Combine PGS 
setwd("~/Autism_heterogeneity/temp_pgs_files/")

data1 = fread("SSC_SA_pst_eff_a1_b0.5_phi1e-02_chr1.txt")

for (i in 2:22){
  data2 = fread(paste0("SSC_SA_pst_eff_a1_b0.5_phi1e-02_chr", i, ".txt"))
  data1 = rbind(data1, data2)
}

write.table(data1, file = "~/Autism_heterogeneity/final_data_forPGS/Autism_SSC_SA_score.txt", row.names = F, col.names = F, quote = F)

rm(list = ls())


## Step 3: Run score in plink (bash)

./plink --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --score ./final_data_forPGS/Autism_SPARK_SA_score.txt 2 4 6 center --out ./PGS/SPARK_SA_finalscore
./plink --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --score ./final_data_forPGS/Autism_SPARK_CT_score.txt 2 4 6 center --out ./PGS/SPARK_CT_finalscore
./plink --bfile ~/SPARK/SPARK_v1v2/SPARK_v1v2_eur_allchrs_hg19 --score ./final_data_forPGS/Autism_SPARK_Volume_score.txt 2 4 6 center --out ./PGS/SPARK_Volume_finalscore
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedall --score ./final_data_forPGS/Autism_SSC_SA_score.txt 2 4 6 center --out ./PGS/SSC_SA_finalscore
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedall --score ./final_data_forPGS/Autism_SSC_CT_score.txt 2 4 6 center --out ./PGS/SSC_CT_finalscore
./plink --bfile ~/SFARI/liftOverPlink/files_imputed/SFARImergedall --score ./final_data_forPGS/Autism_SSC_Volume_score.txt 2 4 6 center --out ./PGS/SSC_Volume_finalscore




## Step 4: Analyse the data, in R
library(MASS)
library(data.table)
library(ggplot2)
library(dplyr)
library(lme4)
library(car)

Mother = fread("~/SPARK/Phenotypes/SPARK_mother.txt", header= T)
setnames(Mother, c("subject_sp_id", "family_id"), c("IID", "FID"))

Father = fread("~/SPARK/Phenotypes/SPARK_father.txt", header= T)
setnames(Father, c("subject_sp_id", "family_id"), c("IID", "FID"))

Cases = fread("~/SPARK/Phenotypes/Autistic_proband.txt", header= T)
setnames(Cases, c("subject_sp_id", "family_id"), c("IID", "FID"))

siblings = fread("~/SPARK/Phenotypes/Nonautistic_sibling.txt", header= T)
setnames(siblings, c("subject_sp_id", "family_id"), c("IID", "FID"))

#Next, read the prs scores

setwd("~/Autism_heterogeneity")

prs_sa = fread("./PGS/SPARK_SA_finalscore.profile", header = TRUE)
setnames(prs_sa,"SCORE", "SA_pgs")
prs_sa = prs_sa[,c("IID", "SA_pgs")]

prs_ct = fread("./PGS/SPARK_CT_finalscore.profile", header = TRUE)
setnames(prs_ct,"SCORE", "CT_pgs")
prs_ct = prs_ct[,c("IID", "CT_pgs")]

prs_vol = fread("./PGS/SPARK_Volume_finalscore.profile", header = TRUE)
setnames(prs_vol,"SCORE", "Volume_pgs")
prs_vol = prs_vol[,c("IID", "Volume_pgs")]



prs_merged = merge(prs_sa, prs_ct, by = "IID")
prs_merged = merge(prs_merged, prs_vol, by = "IID")

PCS = fread("~/Autism_heterogeneity/SPARK_all_PC.txt", header = TRUE)
setnames(PCS, "Sample_name", "IID")

merged = merge(prs_merged, PCS, by = "IID")

merged = unique(merged)

med_hist = fread("~/SPARK/Phenotypes/V5/basic_medical_screening.csv")
setnames(med_hist, "subject_sp_id", "IID")

med_hist2 = med_hist[,c("IID", "family_sf_id", "growth_microceph", "growth_macroceph", "asd", "sex", "age_at_eval_months" )]

med_hist2[is.na(med_hist2)] <- 0 
merged2 = merge(merged, med_hist2, by = "IID")

summary(glm(growth_macroceph ~ scale(Volume_pgs) + scale(CT_pgs) + scale(SA_pgs) + sex + asd + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2, family = binomial))

summary(glm(growth_macroceph ~ scale(CT_pgs) + scale(SA_pgs) + sex + asd + scale(age_at_eval_months) + scale(X1) + scale(X2) + scale(X3) + scale(X4) + scale(X5) + scale(X6) + scale(X7) + scale(X8) + scale(X9) + scale(X10), data = merged2, family = binomial))


