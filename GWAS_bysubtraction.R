##Analysis for conditioning regional effects on global effects##



library(data.table)
library(GenomicSEM)

###Step 1: Make the reference file###

data1 = fread("~/METAL/METAANALYSIS1.TBL")
data1_2 = data1[,c("MarkerName", "Allele1", "Allele2", "Freq1")]
setnames(data1_2, old = c("MarkerName", "Allele1", "Allele2", "Freq1"), new = c("SNP", "A1", "A2", "MAF"))
data3 = SA_meta[,c("SNP", "CHR", "BP")]

merged = merge(data3, data1_2, by = "SNP")
merged2 = merged[,c("SNP", "CHR", "BP", "MAF", "A1", "A2")]
merged2$A1 = toupper(merged2$A1)
merged2$A2 = toupper(merged2$A2)

merged2_major = subset(merged2, MAF > 0.5)
merged2_minor = subset(merged2, MAF < 0.5)
merged2_major$MAF = 1 - merged2_major$MAF

setnames(merged2_major, old = c("A1", "A2"), new = c("A2", "A1"))

merged2_major2 = merged2_major[,c("SNP", "CHR", "BP", "MAF", "A1", "A2")]
merged2_all = rbind(merged2_major2, merged2_minor)

write.table(merged2_all, file = "~/ldsc/temp/reference.txt", row.names = F, col.names = T, quote = F)

rm(list = ls())

####Step 2: Read the files and prepare the SNPs###

Significance_list = fread("~/METAL/MD_meta/final_GWAS/Clumped/Significantlists.txt")
setnames(Significance_list, old = c("BETA_meta", "SE_meta", "P_meta"), new = c("BETA", "SE", "P"))
global_meta = fread("~/METAL/Whole_brain_meta/final_GWAS/MD_meta.txt")
global_significant = global_meta[global_meta$SNP %in% Significance_list$SNP,]

v1 = unique(Significance_list$region)

for (i in v1){
  
  p = paste0("~/METAL/MD_meta/final_GWAS/munged/MD_sumstats_", i, ".sumstats")
  
  traits <- c(p,"~/METAL/Whole_brain_meta/final_GWAS/OD.sumstats")
  sample.prev <- c(NA,NA)
  population.prev <- c(NA,NA)
  ld<-"~/ldsc/eur_w_ld_chr/"
  wld <- "~/ldsc/eur_w_ld_chr/"
  trait.names<-c("regional_GWAS", "global_GWAS")
  
  LDSCoutput <- ldsc(traits, 
                     sample.prev, 
                     population.prev, 
                     ld, 
                     wld, 
                     trait.names)
  
  
  model<-'global =~ NA*regional_GWAS + start(0.4)*global_GWAS
        local=~ NA*regional_GWAS
         
         global~~1*global
         local~~1*local
         global~~0*local

         global_GWAS ~~ 0*regional_GWAS
         global_GWAS~~0*global_GWAS
         regional_GWAS~~0*regional_GWAS'
  
  output<-usermodel(LDSCoutput,estimation="DWLS",model=model)
  
  output
  
  sig_region = subset(Significance_list, region == i)
  sig_region_pruned = sig_region[,1:8]
  global_pruned = global_significant[global_significant$SNP %in% sig_region_pruned$SNP,]
  
  write.table(sig_region_pruned, file = "~/ldsc/temp/regional.txt", row.names = F, col.names = T, quote = F)
  write.table(global_pruned, file = "~/ldsc/temp/global.txt", row.names = F, col.names = T, quote = F)
  
  files = c("~/ldsc/temp/regional.txt", "~/ldsc/temp/global.txt")
  ref = "~/ldsc/temp/reference.txt"
  trait.names = c("global","regional")
  se.logit = c(F,F)
  
  
  
  try(p_sumstats<-sumstats(files, ref, trait.names, se.logit, OLS=c(T,T),linprob=NULL,  N=c(36843,36843)))
  
  model<-'global=~ NA*regional_GWAS + start(0.5)*global_GWAS + start(0.2)*regional_GWAS
        local=~ NA*regional_GWAS + start(0.2)*regional_GWAS
        
        global~SNP
        local~SNP
         
         global~~1*global
         local~~1*local
         global~~0*local

         global_GWAS ~~ 0*regional_GWAS
         global_GWAS~~0*global_GWAS
         regional_GWAS~~0*regional_GWAS
         SNP~~SNP'
  
  #Run the Genomic SEM GWAS
  try(outputGWAS<-userGWAS(covstruc=LDSCoutput,SNPs=p_sumstats,estimation="DWLS",model=model,sub =c("global~SNP","local~SNP"), parallel = TRUE))# printwarn = FALSE
  
  
  try(save(outputGWAS, file = paste0("~/ldsc/genomic_SEM_conditional/MD/outputGWAS_", i, ".rds")))
  rm(outputGWAS)
  rm(model)
  rm(LDSCoutput)
  rm(output)
  rm(p_sumstats)
  
}

rm(list = ls())