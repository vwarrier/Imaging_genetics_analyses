##Genetic correlations across regions and SNP heritability using genomic SEM, R###

library(data.table)
library(gdata)
library(GenomicSEM)


setwd("~/METAL/FA_meta/final_GWAS/munged")

for (i in 1:180){
  data1 = fread(paste0("FA_sumstats_", i, ".sumstats"))
  data1$N = 36000
  write.table(data1, file = paste0("FA_sumstats_", i, ".sumstats"), row.names = F, col.names = T, quote = F, sep = "\t")
} 


traits = c('FA_sumstats_1.sumstats', 'FA_sumstats_2.sumstats', 'FA_sumstats_3.sumstats', 'FA_sumstats_4.sumstats', 'FA_sumstats_5.sumstats', 'FA_sumstats_6.sumstats', 'FA_sumstats_7.sumstats', 'FA_sumstats_8.sumstats', 'FA_sumstats_9.sumstats', 'FA_sumstats_10.sumstats', 'FA_sumstats_11.sumstats', 'FA_sumstats_12.sumstats', 'FA_sumstats_13.sumstats', 'FA_sumstats_14.sumstats', 'FA_sumstats_15.sumstats', 'FA_sumstats_16.sumstats', 'FA_sumstats_17.sumstats', 'FA_sumstats_18.sumstats', 'FA_sumstats_19.sumstats', 'FA_sumstats_20.sumstats',
           'FA_sumstats_21.sumstats', 'FA_sumstats_22.sumstats', 'FA_sumstats_23.sumstats', 'FA_sumstats_24.sumstats', 'FA_sumstats_25.sumstats', 'FA_sumstats_26.sumstats', 'FA_sumstats_27.sumstats', 'FA_sumstats_28.sumstats', 'FA_sumstats_29.sumstats', 'FA_sumstats_30.sumstats', 'FA_sumstats_31.sumstats', 'FA_sumstats_32.sumstats', 'FA_sumstats_33.sumstats', 'FA_sumstats_34.sumstats', 'FA_sumstats_35.sumstats', 'FA_sumstats_36.sumstats', 'FA_sumstats_37.sumstats', 'FA_sumstats_38.sumstats', 'FA_sumstats_39.sumstats', 'FA_sumstats_40.sumstats',
           'FA_sumstats_41.sumstats', 'FA_sumstats_42.sumstats', 'FA_sumstats_43.sumstats', 'FA_sumstats_44.sumstats', 'FA_sumstats_45.sumstats', 'FA_sumstats_46.sumstats', 'FA_sumstats_47.sumstats', 'FA_sumstats_48.sumstats', 'FA_sumstats_49.sumstats', 'FA_sumstats_50.sumstats', 'FA_sumstats_51.sumstats', 'FA_sumstats_52.sumstats', 'FA_sumstats_53.sumstats', 'FA_sumstats_54.sumstats', 'FA_sumstats_55.sumstats', 'FA_sumstats_56.sumstats', 'FA_sumstats_57.sumstats', 'FA_sumstats_58.sumstats', 'FA_sumstats_59.sumstats', 'FA_sumstats_60.sumstats',
           'FA_sumstats_61.sumstats', 'FA_sumstats_62.sumstats', 'FA_sumstats_63.sumstats', 'FA_sumstats_64.sumstats', 'FA_sumstats_65.sumstats', 'FA_sumstats_66.sumstats', 'FA_sumstats_67.sumstats', 'FA_sumstats_68.sumstats', 'FA_sumstats_69.sumstats', 'FA_sumstats_70.sumstats', 'FA_sumstats_71.sumstats', 'FA_sumstats_72.sumstats', 'FA_sumstats_73.sumstats', 'FA_sumstats_74.sumstats', 'FA_sumstats_75.sumstats', 'FA_sumstats_76.sumstats', 'FA_sumstats_77.sumstats', 'FA_sumstats_78.sumstats', 'FA_sumstats_79.sumstats', 'FA_sumstats_80.sumstats',
           'FA_sumstats_81.sumstats', 'FA_sumstats_82.sumstats', 'FA_sumstats_83.sumstats', 'FA_sumstats_84.sumstats', 'FA_sumstats_85.sumstats', 'FA_sumstats_86.sumstats', 'FA_sumstats_87.sumstats', 'FA_sumstats_88.sumstats', 'FA_sumstats_89.sumstats', 'FA_sumstats_90.sumstats', 'FA_sumstats_91.sumstats', 'FA_sumstats_92.sumstats', 'FA_sumstats_93.sumstats', 'FA_sumstats_94.sumstats', 'FA_sumstats_95.sumstats', 'FA_sumstats_96.sumstats', 'FA_sumstats_97.sumstats', 'FA_sumstats_98.sumstats', 'FA_sumstats_99.sumstats', 'FA_sumstats_100.sumstats',
           'FA_sumstats_101.sumstats', 'FA_sumstats_102.sumstats', 'FA_sumstats_103.sumstats', 'FA_sumstats_104.sumstats', 'FA_sumstats_105.sumstats', 'FA_sumstats_106.sumstats', 'FA_sumstats_107.sumstats', 'FA_sumstats_108.sumstats', 'FA_sumstats_109.sumstats', 'FA_sumstats_110.sumstats', 'FA_sumstats_111.sumstats', 'FA_sumstats_112.sumstats', 'FA_sumstats_113.sumstats', 'FA_sumstats_114.sumstats', 'FA_sumstats_115.sumstats', 'FA_sumstats_116.sumstats', 'FA_sumstats_117.sumstats', 'FA_sumstats_118.sumstats', 'FA_sumstats_119.sumstats', 'FA_sumstats_120.sumstats',
           'FA_sumstats_121.sumstats', 'FA_sumstats_122.sumstats', 'FA_sumstats_123.sumstats', 'FA_sumstats_124.sumstats', 'FA_sumstats_125.sumstats', 'FA_sumstats_126.sumstats', 'FA_sumstats_127.sumstats', 'FA_sumstats_128.sumstats', 'FA_sumstats_129.sumstats', 'FA_sumstats_130.sumstats', 'FA_sumstats_131.sumstats', 'FA_sumstats_132.sumstats', 'FA_sumstats_133.sumstats', 'FA_sumstats_134.sumstats', 'FA_sumstats_135.sumstats', 'FA_sumstats_136.sumstats', 'FA_sumstats_137.sumstats', 'FA_sumstats_138.sumstats', 'FA_sumstats_139.sumstats', 'FA_sumstats_140.sumstats',
           'FA_sumstats_141.sumstats', 'FA_sumstats_142.sumstats', 'FA_sumstats_143.sumstats', 'FA_sumstats_144.sumstats', 'FA_sumstats_145.sumstats', 'FA_sumstats_146.sumstats', 'FA_sumstats_147.sumstats', 'FA_sumstats_148.sumstats', 'FA_sumstats_149.sumstats', 'FA_sumstats_150.sumstats', 'FA_sumstats_151.sumstats', 'FA_sumstats_152.sumstats', 'FA_sumstats_153.sumstats', 'FA_sumstats_154.sumstats', 'FA_sumstats_155.sumstats', 'FA_sumstats_156.sumstats', 'FA_sumstats_157.sumstats', 'FA_sumstats_158.sumstats', 'FA_sumstats_159.sumstats', 'FA_sumstats_160.sumstats',
           'FA_sumstats_161.sumstats', 'FA_sumstats_162.sumstats', 'FA_sumstats_163.sumstats', 'FA_sumstats_164.sumstats', 'FA_sumstats_165.sumstats', 'FA_sumstats_166.sumstats', 'FA_sumstats_167.sumstats', 'FA_sumstats_168.sumstats', 'FA_sumstats_169.sumstats', 'FA_sumstats_170.sumstats', 'FA_sumstats_171.sumstats', 'FA_sumstats_172.sumstats', 'FA_sumstats_173.sumstats', 'FA_sumstats_174.sumstats', 'FA_sumstats_175.sumstats', 'FA_sumstats_176.sumstats', 'FA_sumstats_177.sumstats', 'FA_sumstats_178.sumstats', 'FA_sumstats_179.sumstats', 'FA_sumstats_180.sumstats')
sample.prev <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
population.prev <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)

ld <- "~/ldsc/eur_w_ld_chr/"

wld <- "~/ldsc/eur_w_ld_chr/"

trait.names<-c('FA_1', 'FA_2', 'FA_3', 'FA_4', 'FA_5', 'FA_6', 'FA_7', 'FA_8', 'FA_9', 'FA_10', 'FA_11', 'FA_12', 'FA_13', 'FA_14', 'FA_15', 'FA_16', 'FA_17', 'FA_18', 'FA_19', 'FA_20', 'FA_21', 'FA_22', 'FA_23', 'FA_24', 'FA_25', 'FA_26', 'FA_27', 'FA_28', 'FA_29', 'FA_30', 'FA_31', 'FA_32', 'FA_33', 'FA_34', 'FA_35', 'FA_36', 'FA_37', 'FA_38', 'FA_39', 'FA_40', 'FA_41', 'FA_42', 'FA_43', 'FA_44', 'FA_45', 'FA_46', 'FA_47', 'FA_48', 'FA_49', 'FA_50', 'FA_51', 'FA_52', 'FA_53', 'FA_54', 'FA_55', 'FA_56', 'FA_57', 'FA_58', 'FA_59', 'FA_60', 'FA_61', 'FA_62', 'FA_63', 'FA_64', 'FA_65', 'FA_66', 'FA_67', 'FA_68', 'FA_69', 'FA_70', 'FA_71', 'FA_72', 'FA_73', 'FA_74', 'FA_75', 'FA_76', 'FA_77', 'FA_78', 'FA_79', 'FA_80', 'FA_81', 'FA_82', 'FA_83', 'FA_84', 'FA_85', 'FA_86', 'FA_87', 'FA_88', 'FA_89', 'FA_90', 'FA_91', 'FA_92', 'FA_93', 'FA_94', 'FA_95', 'FA_96', 'FA_97', 'FA_98', 'FA_99', 'FA_100', 'FA_101', 'FA_102', 'FA_103', 'FA_104', 'FA_105', 'FA_106', 'FA_107', 'FA_108', 'FA_109', 'FA_110', 'FA_111', 'FA_112', 'FA_113', 'FA_114', 'FA_115', 'FA_116', 'FA_117', 'FA_118', 'FA_119', 'FA_120', 'FA_121', 'FA_122', 'FA_123', 'FA_124', 'FA_125', 'FA_126', 'FA_127', 'FA_128', 'FA_129', 'FA_130', 'FA_131', 'FA_132', 'FA_133', 'FA_134', 'FA_135', 'FA_136', 'FA_137', 'FA_138', 'FA_139', 'FA_140', 'FA_141', 'FA_142', 'FA_143', 'FA_144', 'FA_145', 'FA_146', 'FA_147', 'FA_148', 'FA_149', 'FA_150', 'FA_151', 'FA_152', 'FA_153', 'FA_154', 'FA_155', 'FA_156', 'FA_157', 'FA_158', 'FA_159', 'FA_160', 'FA_161', 'FA_162', 'FA_163', 'FA_164', 'FA_165', 'FA_166', 'FA_167', 'FA_168', 'FA_169', 'FA_170', 'FA_171', 'FA_172', 'FA_173', 'FA_174', 'FA_175', 'FA_176', 'FA_177', 'FA_178', 'FA_179', 'FA_180')
LDSCOutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)

save(LDSCOutput, file = "LDSCoutput.RData")

correlation_matrix = cov2cor(LDSCOutput$S)
save(correlation_matrix, file = "genetic_correlation_matrix.RData")

#create Standardized S (genetic correlation matrix) using matrix algebra
D=sqrt(diag(diag(LDSCOutput$S)))
S_Stand=solve(D)%*%LDSCOutput$S%*%solve(D)

#obtain diagonals of the original V matrix and take their sqrt to get SE's
Dvcov1<-sqrt(diag(LDSCOutput$V))

#calculate the ratio of the rescaled and original S matrices
scaleO=as.vector(lowerTriangle((S_Stand/LDSCOutput$S),diag=T))

## Make sure that if ratio in NaN (division by zero) we put the zero back in: 
scaleO[is.nan(scaleO)] <- 0

#rescale the SEs by the same multiples that the S matrix was rescaled by
Dvcovl<-as.vector(Dvcov1*t(scaleO))

#obtain the sampling correlation matrix by standardizing the original V matrix
Vcor<-cov2cor(LDSCOutput$V)

#rescale the sampling correlation matrix by the appropriate diagonals
V_stand<-diag(Dvcovl)%*%Vcor%*%diag(Dvcovl)

#create empty SE matrix to store values
SE<-matrix(0, ncol(LDSCOutput$S), ncol(LDSCOutput$S))

#pull out SEs of genetic correlations (sqrt of diagonals of standardized V)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(V_stand))

save(SE, file = "rg_SE.RData")

heritability = as.data.frame(diag(LDSCOutput$S))

#create empty SE matrix to store values
SE_cov<-matrix(0, ncol(LDSCOutput$S), ncol(LDSCOutput$S))

#pull out SE for genetic covariances
SE_cov[lower.tri(SE_cov,diag=TRUE)] <-sqrt(diag(LDSCOutput$V))

heritability$se = diag(SE_cov)

setnames(heritability, 1, "heritability")

write.table(heritability, file = "SNPheritability.txt", row.names = F, col.names = T, quote = F)