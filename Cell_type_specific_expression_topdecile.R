#HMAGMA_celltypes

#We are first going to create multiple different datasets and save them for posterity. In all datasets, we will: 
#1. Calculate log-transformed of the gene expression data + 1 pseudocount, and then center the data (not scale it)
#2. Calculate the total mean and cell-type specific mean value
#3. Conduct a multivariate regression analyses with MAGMA/HMAGMA gene signficiance as a binary indepedent variable,using MAGMA

#Datasets used in this study are as follows:
# 1. Psychencode fetal
# 2. Psychencode adult
# 3. Nowakovski fetal
# 4. Midfetal from Geschwind lab


# Step 1
library(data.table); library(dplyr); library(tidyr)

##################Fetal PsychEncode#############################################################################################
#Step 1 organize the data
load("~/HMAGMA/Psychencode_singlecell/Sestan.fetalHuman.Psychencode.Rdata")

cpm2 = as.data.frame(cpm2)
cpm2$Geneid = row.names(cpm2)
cpm2 = separate(data = cpm2, col = Geneid, into = c("GENE", "HGNC_name"), sep = "\\|")
data_log = as.data.frame(log2(cpm2[,1:762] + 1))
data_scaled = as.data.frame(scale(data_log, center = T, scale=FALSE) + 1)
names = as.data.frame(cpm2[,763])
data_scaled = cbind(names, data_scaled)
setnames(data_scaled, 1, "GENE")

meta2$names = row.names(meta2)

mean_data = as.data.frame(data_scaled$GENE)

alpha = unique(meta2$ctype)

for (j in alpha){
  id_list = subset(meta2, ctype == j)
  data_scaled_2 = data_scaled[,colnames(data_scaled) %in% id_list$names]
  mean_data[[j]] = rowMeans(data_scaled_2)
}

mean_data$total_mean = rowMeans(data_scaled[,-1])
setnames(mean_data, "data_scaled$GENE", "GENE")

write.table(mean_data, file = "~/HMAGMA/Psychencode_fetal_means.txt", row.names = F, col.names = T, quote = F)


###Step 2: Obtain genes in the top decile by cell type


mean_data = fread("~/HMAGMA/Psychencode_fetal_means.txt")
alpha = names(mean_data)
alpha = alpha[-1]
alpha = alpha[-25]

specific_means = as.data.frame(mean_data$GENE)
specific_mean_genenames = matrix(0, ncol = 1, nrow = 6015)
specific_mean_genenames = data.frame(specific_mean_genenames)


for (j in alpha){
  specific_means[[j]] = mean_data[[j]]/mean_data$total_mean
  specific_means[[j]] = ntile(specific_means[[j]], 10)
  specific_means[[j]] = ifelse(specific_means[[j]] == 10, 1, 0)
  name = paste0(j)
  check1 = subset(specific_means, specific_means[[j]] == 1 )
  specific_mean_genenames = cbind(specific_mean_genenames, check1$`mean_data$GENE`)
}

specific_mean_genenames = specific_mean_genenames[,-1]
names(specific_mean_genenames) = alpha

specific_mean_genenames2 = t(specific_mean_genenames)
specific_mean_genenames2 = data.frame(specific_mean_genenames2)
row.names(specific_mean_genenames2) = names(specific_mean_genenames)
write.table(specific_mean_genenames2, file = "~/HMAGMA/Results/Psychencode_fetal_topdecile.txt", row.names = T, col.names = F, quote = F)



##################Adult PsychEncode#############################################################################################

load("~/HMAGMA/Psychencode_singlecell/Sestan.adultHumanNuclei.Psychencode.Rdata")

umi2 = as.data.frame(umi2)
umi2$Geneid = row.names(umi2)
umi2 = separate(data = umi2, col = Geneid, into = c("GENE", "HGNC_name"), sep = "\\|")
data_log = as.data.frame(log2(umi2[,1:17093] + 1))
data_scaled = as.data.frame(scale(data_log, center = T, scale=FALSE) + 1)
names = as.data.frame(umi2[,17094])
names = separate(data = names, col = `umi2[, 17094]`, into = c("GENE", "dot"), sep = "\\.")
names = as.data.frame(names$GENE)
setnames(names, 1, "GENE")
data_scaled = cbind(names, data_scaled)


meta2$names = row.names(meta2)

mean_data = as.data.frame(data_scaled$GENE)
alpha = unique(meta2$ctype)

for (j in alpha){
  id_list = subset(meta2, ctype == j)
  data_scaled_2 = data_scaled[,colnames(data_scaled) %in% id_list$names]
  mean_data[[j]] = rowMeans(data_scaled_2)
}

mean_data$total_mean = rowMeans(data_scaled[,-1])
setnames(mean_data, "data_scaled$GENE", "GENE")


write.table(mean_data, file = "~/HMAGMA/Psychencode_singlecell/celltype_adult_mean.txt", row.names = F, col.names = T, quote = F)


###Now run it: Step 2

specific_mean = fread("~/HMAGMA/celltype_adult_specificmean.txt")
alpha = names(specific_mean)
alpha = alpha[-1]
alpha = alpha[-1]

mean_data = specific_mean
rm(specific_mean)

specific_means = as.data.frame(mean_data$GENE)
specific_mean_genenames = matrix(0, ncol = 1, nrow = 2347)
specific_mean_genenames = data.frame(specific_mean_genenames)


for (j in alpha){
  specific_means[[j]] = mean_data[[j]]/mean_data$total_mean
  specific_means[[j]] = ntile(specific_means[[j]], 10)
  specific_means[[j]] = ifelse(specific_means[[j]] == 10, 1, 0)
  name = paste0(j)
  check1 = subset(specific_means, specific_means[[j]] == 1 )
  specific_mean_genenames = cbind(specific_mean_genenames, check1$`mean_data$GENE`)
}

specific_mean_genenames = specific_mean_genenames[,-1]
names(specific_mean_genenames) = alpha

specific_mean_genenames2 = t(specific_mean_genenames)
specific_mean_genenames2 = data.frame(specific_mean_genenames2)
row.names(specific_mean_genenames2) = names(specific_mean_genenames)
write.table(specific_mean_genenames2, file = "~/HMAGMA/Results/Psychencode_adult_topdecile.txt", row.names = T, col.names = F, quote = F)




###Kreigstein early fetal - 2021 Nat Neuro

data = fread("./Eze_earlyhumandev_2021/exprMatrix.tsv")

id = fread("./Eze_earlyhumandev_2021/meta.tsv")

##remove those pesky infinity rows, the data is already scaled
means_alpa = as.data.frame(colSums (data[,-1], na.rm = FALSE, dims = 1))
check1 = subset(means_alpa, means != "Inf")
check2 = as.data.frame(rownames(check1))

data = as.data.frame(data)
data1 = data[,colnames(data) %in% check2$`rownames(check1)`]

mean_data = as.data.frame(data$gene)
alpha = unique(id$Cell.Type)

for (j in alpha){
  id_list = subset(id,Cell.Type == j)
  data_scaled_2 = data1[,colnames(data1) %in% id_list$Cell]
  mean_data[[j]] = rowMeans(data_scaled_2, na.rm = TRUE)
}

mean_data$total_mean = rowMeans(data1[,-1], na.rm = TRUE)

file1 = fread("mart_export (3).txt")
merged = merge(mean_data, file1, by.x = "data$gene", by.y = "ID")
specific_mean = merged[,c(2:9)]

save(specific_mean, file = "./Eze_earlyhumandev_2021/Mean_data.RData")

mean_data2 = as.data.frame(data$gene)
alpha2 = unique(id$Cluster)

for (j in alpha2){
  id_list = subset(id,Cluster == j)
  data_scaled_2 = data1[,colnames(data1) %in% id_list$Cell]
  mean_data2[[j]] = rowMeans(data_scaled_2, na.rm = TRUE)
  
}

mean_data2$`data$gene` = data$gene
mean_data2$total_mean = rowMeans(data1[,-1], na.rm = TRUE)

merged = merge(mean_data2, file1, by.x = "data$gene", by.y = "ID")
specific_mean_cluster = merged[,c(2:63)]

save(specific_mean_cluster, file = "./Eze_earlyhumandev_2021/Mean_data_cluster.RData")



#######Run it now###############
Protein_coding = fread("~/HMAGMA/protein_coding_IDs.txt")
setnames(Protein_coding, 1, "GENE")

load("~/HMAGMA/Eze_earlyhumandev_2021/Mean_data.RData")

id = fread("~/HMAGMA/Eze_earlyhumandev_2021/meta.tsv")

alpha = unique(id$Cell.Type)


mean_data = specific_mean
rm(specific_mean)

specific_means = as.data.frame(mean_data$GENE)
specific_mean_genenames = matrix(0, ncol = 1, nrow = 1737)
specific_mean_genenames = data.frame(specific_mean_genenames)


for (j in alpha){
  specific_means[[j]] = mean_data[[j]]/mean_data$total_mean
  specific_means[[j]] = ntile(specific_means[[j]], 10)
  specific_means[[j]] = ifelse(specific_means[[j]] == 10, 1, 0)
  name = paste0(j)
  check1 = subset(specific_means, specific_means[[j]] == 1 )
  specific_mean_genenames = cbind(specific_mean_genenames, check1$`mean_data$GENE`)
}

specific_mean_genenames = specific_mean_genenames[,-1]
names(specific_mean_genenames) = alpha

specific_mean_genenames2 = t(specific_mean_genenames)
specific_mean_genenames2 = data.frame(specific_mean_genenames2)
row.names(specific_mean_genenames2) = names(specific_mean_genenames)
write.table(specific_mean_genenames2, file = "~/HMAGMA/Results/Eze_fetal_topdecile.txt", row.names = T, col.names = F, quote = F)




###Geschwind midfetal

library(Matrix)
load("~/HMAGMA/sc_dev_cortex_geschwind/raw_counts_mat.rdata")
raw.counts.mat <- as.matrix(raw_counts_mat)

id = read.csv("~/HMAGMA/sc_dev_cortex_geschwind/cell_metadata.csv")

##remove those pesky infinity rows, the data is already scaled

data_log = as.data.frame(log2(raw.counts.mat[,-1] + 1))
data_scaled = as.data.frame(scale(data_log, center = T, scale=FALSE) + 1)
names = as.data.frame(row.names(raw.counts.mat))
data_scaled = cbind(names, data_scaled[,-1])
setnames(data_scaled, 1, "GENE")

mean_data = as.data.frame(data_scaled$GENE)
alpha = unique(id$Cluster)

for (j in alpha){
  id_list = subset(id,Cluster == j)
  data_scaled_2 = data_scaled[,colnames(data_scaled) %in% id_list$Cell]
  mean_data[[j]] = rowMeans(data_scaled_2, na.rm = TRUE)
}

mean_data$total_mean = rowMeans(data_scaled[,-1], na.rm = TRUE)

setnames(mean_data, 1, "GENE")

specific_mean = mean_data

gene_names = fread("~/HMAGMA/sc_dev_cortex_geschwind/gene_names.txt")

specific_mean = merge(gene_names, specific_mean, by.x = "HGNC", by.y = "GENE")

specific_mean2 = specific_mean[!duplicated(specific_mean[,c("HGNC")]),]
specific_mean2 = specific_mean2[!duplicated(specific_mean2[,c("Ensembl")]),]

specific_mean = specific_mean2[,c(-1)]
setnames(specific_mean, 1, "GENE")


save(specific_mean, file = "~/HMAGMA/sc_dev_cortex_geschwind/Mean_data.RData")

table1 = NULL


