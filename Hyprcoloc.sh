###Hyprcoloc scripts###
###############################################
###1. Script to identify blocks for Hyprcoloc##
##############################################

data1 = read.delim("clipboard") ##EWS significant SNPs by phenotype (i.e, all global, or all SA regional)
#SNP CHR BP

data2 = read.delim("clipboard") #LDbins

merged = NULL

for (row in 1:nrow(data1)){
  data3 = data1[row,]
  data4 = subset(data2, CHR == data3$CHR)
  data4 = subset(data4, start < data3$BP & stop > data3$BP)
  data5 = cbind(data3, data4)
  merged = rbind(merged, data5)
}

merged2 = merged[!duplicated(merged$Block_number),]

View(merged2)

write.table(merged2, file = "Hypercoloc_data_whole_brain.txt", row.names = F, col.names = T, quote = F)

###############################################
###2. Global analyses: read all the files##
###############################################


setwd("~/METAL/Whole_brain_meta/final_GWAS")

file_names = list("CT", "Volume", "Foldingindex", "Meancurvature",
                  "Gaussiancurvature", "Intrinsic", "LGI", 
                  "ISOVF", "ICVF", "OD", "FA", "MD")
data1 = fread("SA_meta.txt")


hypercoloc_data = fread("~/Hyprcoloc/Global/Hypercoloc_data_whole_brain.txt")
hypercoloc_data = hypercoloc_data2

merged = NULL

for (row in 1:nrow(hypercoloc_data)){
  data3 = hypercoloc_data[row,]
  data4 = subset(data1, CHR == data3$CHR)
  data4 = subset(data4, BP > data3$start & BP < data3$stop)
  data4$Block_number = data3$Block_number
  merged = rbind(merged, data4)
}


merged = merged[!duplicated(merged$SNP),]
merged = data1[,c("SNP", "BETA", "SE", "Block_number")]
setnames(merged, "BETA", "SA_BETA")
setnames(merged, "SE", "SA_SE")



for (i in file_names){
  data1 = fread(paste0(i, "_meta.txt"))
  data1 = data1[!duplicated(data1$SNP),]
  merged1 = data1[,c("SNP", "BETA", "SE")]
  setnames(merged1, "SE", paste0(i, "_SE"))
  setnames(merged1, "BETA", paste0(i, "_BETA"))
  merged = merge(merged, merged1, by = "SNP")
}

merged = na.omit(merged)

save(merged, file  = "~/Hyprcoloc/Global/data_forcoloc.RDa")


###############################################
###3. Run Hyprcoloc##
###############################################

load("~/Hyprcoloc/Global/data_forcoloc.RDa")
omega = unique(merged$Block_number)
alpha = list(omega)
alpha = list(alpha)

for (i in omega){
  beta = paste0(i)
  data1 = subset(merged, Block_number == beta)
  
  BETA3 = data1[,c("SA_BETA", "CT_BETA", "Volume_BETA", "Foldingindex_BETA", "Meancurvature_BETA", "Gaussiancurvature_BETA",
                    "Intrinsic_BETA", "LGI_BETA", "ISOVF_BETA", "ICVF_BETA", "OD_BETA", "FA_BETA", "MD_BETA")]
  BETA4 = as.matrix(BETA3)
  rownames(BETA4) = data1$SNP
  
  SE3 = data1[,c("SA_SE", "CT_SE", "Volume_SE", "Foldingindex_SE", "Meancurvature_SE", "Gaussiancurvature_SE",
                  "Intrinsic_SE", "LGI_SE", "ISOVF_SE", "ICVF_SE", "OD_SE", "FA_SE", "MD_SE")]
  SE4 = as.matrix(SE3)
  rownames(SE4) = data1$SNP
  
  BETA4[BETA4==0] <- 0.00000001
  SE4[SE4==0]<- 0.00000001
  
  traits <- names(BETA3)
  rsid = rownames(BETA4)
  
  res = hyprcoloc(BETA4, SE4, trait.names=traits, snp.id=rsid, uniform.priors = FALSE)
  
  save(res, file = paste0("~/Hyprcoloc/Global/", i, "_resultsglobal.RDa"))
}


###############################################
###4. Regional analyses: read all the files##
###############################################


library(data.table)
library(hyprcoloc)
setwd("~/METAL/Folding_meta/final_GWAS/")

data1 = fread("~/METAL/Whole_brain_meta/final_GWAS/Foldingindex_meta.txt")
hypercoloc_data = fread("~/Hyprcoloc/Regional_FI/Hypercoloc_data_FI.txt")

merged = NULL

for (row in 1:nrow(hypercoloc_data)){
  data3 = hypercoloc_data[row,]
  data4 = subset(data1, CHR == data3$CHR)
  data4 = subset(data4, BP > data3$start & BP < data3$stop)
  data4$Block_number = data3$Block_number
  merged = rbind(merged, data4)
}


merged = merged[!duplicated(merged$SNP),]
merged = merged[,c("SNP", "BETA", "SE", "Block_number")]
setnames(merged, "BETA", "Global_BETA")
setnames(merged, "SE", "Global_SE")



for (i in 1:180){
  data1 = fread(paste0("Folding_plinkmeta", i, ".txt"))
  data1 = data1[!duplicated(data1$SNP),]
  merged1 = data1[,c("SNP", "BETA", "SE")]
  setnames(merged1, "SE", paste0(i, "_SE"))
  setnames(merged1, "BETA", paste0(i, "_BETA"))
  merged = merge(merged, merged1, by = "SNP")
}

merged = na.omit(merged)

save(merged, file  = "~/Hyprcoloc/Regional_FI/data_forcoloc.RDa")

rm(list = ls())


###############################################
###5. Run Hyprcoloc##
###############################################

load("~/Hyprcoloc/Regional_FI/data_forcoloc.RDa")
omega = unique(merged$Block_number)


for (i in omega){
  beta = paste0(i)
  data1 = subset(merged, Block_number == beta)
  
  BETA3 =  data1 %>% dplyr:: select(contains("BETA"))
  
  BETA4 = as.matrix(BETA3)
  rownames(BETA4) = data1$SNP
  
  SE3 = data1 %>% dplyr:: select(contains("SE"))
  
  SE4 = as.matrix(SE3)
  rownames(SE4) = data1$SNP
  
  BETA4[BETA4==0] <- 0.00000001
  SE4[SE4==0]<- 0.00000001
  
  traits <- names(BETA3)
  rsid = rownames(BETA4)
  
  res = hyprcoloc(BETA4, SE4, trait.names=traits, snp.id=rsid, uniform.priors = FALSE)
  
  save(res, file = paste0("~/Hyprcoloc/Regional_FI/", i, "_resultsglobal.RDa"))
}

