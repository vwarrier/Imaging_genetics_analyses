
##QC data in R
## Files for WM 
# R program ukb22918.tab created 2018-09-17 by ukb2r.cpp Mar 14 2018 14:22:05
library(data.table)
library(plyr)

bd <- fread("~/UKB_v2/ukb24196.tab", header=TRUE, sep="\t")

## subset it to self-identified white europeans
file1 = subset(bd, f.21000.0.0 == "1"| f.21000.0.0 == "1001" | f.21000.0.0 == "1002" | f.21000.0.0 == "1003"| f.21000.0.0 == "1004") 

file1$checksex = ifelse(file1$f.31.0.0 == file1$f.22001.0.0, "correct", "incorrect")
file1 = subset(file1, checksex == "correct") # remove sex mismatches
pheno2 = file1[is.na(file1$f.22027.0.0),] #remove excessive heterozygosity



PC1_mean = mean(na.omit(pheno2$f.22009.0.1))
PC1_sd = sd(na.omit(pheno2$f.22009.0.1))
PC2_mean = mean(na.omit(pheno2$f.22009.0.2))
PC2_sd = sd(na.omit(pheno2$f.22009.0.2))

pheno3 = subset(pheno2, f.22009.0.1 < (PC1_mean + 5*PC1_sd) & f.22009.0.1 > (PC1_mean - 5*PC1_sd))
pheno2 = subset(pheno3, f.22009.0.2 < (PC2_mean + 5*PC2_sd) & f.22009.0.2 > (PC2_mean - 5*PC2_sd))

pheno_covar = fread("~/UKB_v2/UKB_WM/UKB_whitematterdata.txt")

setnames(pheno_covar, "ID", "f.eid")

pheno_covar = pheno_covar[(pheno_covar$f.eid %in% pheno2$f.eid),]


####FILES FOR fastGWA###

pheno = pheno_covar[,c("f.eid", "f.eid", "meanNODDI_ICVF", "meanNODDI_ISOVF", "meanNODDI_OD", "meanFA", "meanMD")]
setnames(pheno, 1, "FID")
setnames(pheno, 2, "IID")

pheno_numeric = pheno[,c(3:7)]
pheno_names = pheno[,c(1,2)]

pheno_numeric = data.frame(scale(pheno_numeric))
pheno_numeric[pheno_numeric > 5] <- NA   
pheno_numeric[pheno_numeric < -5] <- NA  

pheno = cbind(pheno_names, pheno_numeric)

write.table(pheno, file = "~/UKB_v2/UKB_WM/WMpheno_scaled.txt", row.names = F, col.names = T, quote = F)

covar = pheno_covar[,c("f.eid", "f.eid","Scan.Site", "Age", "Agesquared", "Sex", "Age_Sex", 
                       "Agesquared_sex", "X22009.0.1", "X22009.0.2", "X22009.0.3", "X22009.0.4", "X22009.0.5",
                       "X22009.0.6","X22009.0.7","X22009.0.8", "X22009.0.9", "X22009.0.10", "X22009.0.11",
                       "X22009.0.12", "X22009.0.13", "X22009.0.14", "X22009.0.15", "X22009.0.16", "X22009.0.17",
                       "X22009.0.18", "X22009.0.19","X22009.0.20" , "X22009.0.21", "X22009.0.22", "X22009.0.23",
                       "X22009.0.24", "X22009.0.25", "X22009.0.26", "X22009.0.27", "X22009.0.28", "X22009.0.29",
                       "X22009.0.30", "X22009.0.31", "X22009.0.32", "X22009.0.33", "X22009.0.34", "X22009.0.35",
                       "X22009.0.36", "X22009.0.37", "X22009.0.38", "X22009.0.39", "X22009.0.40", "eulerleftplusright",
                       "fd", "fd_max")]

setnames(covar, 1, "FID")
setnames(covar, 2, "IID")

write.table(covar, file = "~/UKB_v2/UKB_WM/WMcovarbolt.txt", row.names = F, col.names = T, quote = F)

covar_discrete = covar[,c("Scan.Site", "Sex")]
qcovar = covar[,-c("Scan.Site", "Sex")]
covar_discrete = covar[,c("FID", "IID","Scan.Site", "Sex")]

qcovar_numeric = qcovar[,c(3:49)]
qcovar_names = pheno[,c(1,2)]

qcovar_numeric = data.frame(scale(qcovar_numeric))
qcovar = cbind(qcovar_names, qcovar_numeric)

write.table(covar_discrete, file = "~/UKB_v2/UKB_WM/WMcovardiscrete.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar, file = "~/UKB_v2/UKB_WM/WMqcovar.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar2, file = "~/UKB_v2/UKB_WM/WMqcovar2.txt", row.names = F, col.names = T, quote = F)


### Create Plink files and then GRM
# First create a list of SNPs that have an MAF > 0.1% and imputation r2 > 0.4

##R

maf = fread("~/UKB_v2/maf_info/maf_imp_OKfile.txt")
extract_snps = maf[,c("V2")]
write.table(extract_snps, file = "~/UKB_v2/UKB_WM/extract_snps.txt", row.names = F, col.names = F, quote = F)

##bash
for i in {1..22}; do ./plink2 --bgen ukb_imp_chr${i}_v3.bgen --keep ~/UKB_v2/UKB_WM/gctakeepfile.txt --maf 0.001 --hwe 0.000001 --make-bed  --extract ~/UKB_v2/UKB_WM/extract_snps.txt --out ukbchr_v2_r2correct_v2_${i} --sample ukb20904_imp_chr1_v3_s487334.sample --threads 15; done
for i in {1..22}; do ./plink2 --bfile ukbchr_v2_${i} --extract ~/UKB_v2/UKB_WM/extract_snps.txt --make-bed --out ukbchr_v2_r2correct_${i} --threads 15; done

for i in {1..22}; do ./gcta64 --bfile ukbchr_v2_${i} --make-grm --out GRM_chr${i} --thread-num 10; done
./gcta64  --mgrm GRM_merge.txt  --make-grm  --out full_grm --thread-num 15
./gcta64  --grm full_grm --grm-cutoff 0.05 --make-grm --out unrelated_grm --thread-num 15

./gcta64 --grm full_grm --make-bK-sparse 0.05 --out sp_grm
./gcta64 --grm unrelated_grm --make-bK-sparse 0.05 --out sp_grm_unrelated



###Run fast GWAS, bash
# Whole brain DTI
for i in {1..5}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/WMpheno_scaled.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete.txt --threads 20 --out ./GWAS/Whole_brain/mWM_GWAS$i; done

# Whole brain Structural T1 and T2 with T1 T2 covariate
for i in {1..8}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/Structural_WB.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete_structural.txt --threads 20 --out ./GWAS/Whole_brain/Structural_GWAS_$i; done


##ROI
# ISOVF (change as needed)
for i in {1..180}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/ISOVF.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete.txt --threads 10 --out ./GWAS/NODDI_ISOVF/GWAS$i; done


### X chromosome analyses, bash
## Create Plink files, bash
./plink2 --bgen ukb_imp_chrX_v3.bgen --geno 0.05 --hwe 0.000001 --maf 0.0001 --make-bed --mind 0.05 --out ukbchrlowmafX --sample ukb20904_imp_chrX_v3_s486676.sample --threads 15
./plink2 --bfile ukbchrlowmafX --keep ~/UKB_v2/UKB_WM/gctakeepfile.txt --maf 0.001 --hwe 0.000001 --make-bed --out ./UKB_WM/ukbchr_v2_X --threads 15
./plink2 --bfile ./UKB_WM/ukbchr_v2_X --extract ~/UKB_v2/UKB_WM/extract_snps.txt --make-bed --out ~/UKB_v2/UKB_WM/GCTA/ukbchr_v2_r2correct_X --threads 15

###There are few samples in the --keep file that are absent from the fam file. Let's find that and create a keep-file, R
fam_file = fread("~/UKB_v2/UKB_WM/GCTA/ukbchr_v2_r2correct_X.fam")
x_keep = fam_file[,c(1:2)]
write.table(x_keep, file = "~/UKB_v2/UKB_WM/chrX.idlist", row.names = F, col.names = F, quote = F)

#Whole brain structural for X chromosome, bash
for i in {1..8}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/Structural_WB.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete_structural.txt --threads 20 --model-only --keep ~/UKB_v2/UKB_WM/chrX.idlist --out ./GWAS/Whole_brain/Structural_GWAS_model_$i; done
for i in {1..}; do ./gcta64 --bfile ~/UKB_v2/UKB_WM/GCTA/ukbchr_v2_r2correct_X --load-model ./GWAS/Whole_brain/Structural_GWAS_model_$i.fastGWA  --geno 0.1 --out ./GWAS/Whole_brain/Structural_GWAS_X_$i --threads 10; done

#Whole brain WM for X chromosome, bash
for i in {1..5}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/WMpheno_scaled.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete.txt --threads 20 --model-only --keep ~/UKB_v2/UKB_WM/chrX.idlist --out ./GWAS/Whole_brain/WM_GWAS_model_$i; done
for i in {1..5}; do ./gcta64 --bfile ~/UKB_v2/UKB_WM/GCTA/ukbchr_v2_r2correct_X --load-model ./GWAS/Whole_brain/WM_GWAS_model_$i.fastGWA  --geno 0.1 --out ./GWAS/Whole_brain/WM_GWAS_X_$i --threads 10; done


##ROI X chromosome
for i in {1..180}; do ./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm  --pheno ~/UKB_v2/UKB_WM/ICVF.txt --mpheno $i --qcovar ~/UKB_v2/UKB_WM/WMqcovar.txt --covar ~/UKB_v2/UKB_WM/WMcovardiscrete.txt --threads 15 --model-only --keep ~/UKB_v2/UKB_WM/chrX.idlist --out ./GWAS/NODDI_ICVF/GWAS_model_X$i; done
for i in {1..180}; do ./gcta64 --bfile ~/UKB_v2/UKB_WM/GCTA/ukbchr_v2_r2correct_X --load-model ./GWAS/NODDI_ICVF/GWAS_model_X$i.fastGWA  --geno 0.1 --out ./GWAS/NODDI_ICVF/WM_GWAS_X_$i --threads 15; done



# Combine X chromosome and autosomes in R
setwd()

for (i in 1:180){
  autosome = fread(paste0("GWAS",i, ".fastGWA"))
  x_chromosome = fread(paste0("WM_GWAS_X_",i, ".fastGWA"))
  GWAS = rbind(autosome, x_chromosome)
  write.table(GWAS, file = paste0("./Final_GWAS/GWAS",i,".txt"), row.names = F, col.names = T, quote = F)
}



### Meta-analysis, bash

for i in {1..180}; do ./plink  --meta-analysis /mnt/b2/home4/arc/vw260/UKB_v2/UKB_WM/GCTA/GWAS/LGI/Final_GWAS/GWAS${i}.txt /mnt/b2/home4/arc/vw260/ABCD/ABCDgenotype/GWAS/Cortical_GWAS/LGI/Final_GWAS/GWAS${i}.txt + qt report-all   --meta-analysis-bp-field POS   --out ./LGI_meta/LGI_plinkmeta${i}; done


### Create manhattan, R

library(qqman)


GWAS = fread("~/UKB_v2/UKB_WM/GCTA/GWAS/Whole_brain/Structural_GWAS_X_7.fastGWA")

nom = subset(GWAS, P < 0.001)

manhattan(nom, chr = "CHR", bp = "POS", p = "P", snp = "SNP",
          col = c("gray10", "gray60"), chrlabs = NULL,
          suggestiveline = -log10(1e-06), genomewideline = -log10(5e-08),
          highlight = NULL, logp = TRUE, annotatePval = NULL,
          annotateTop = TRUE)




