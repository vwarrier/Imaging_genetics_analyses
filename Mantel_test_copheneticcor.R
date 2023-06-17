##Mantel tests###
setwd("/mnt/beegfs/home4/arc/vw260/METAL")

distance = read.csv("geodesic_distance_HCP.csv", header = FALSE)

distance2 = as.matrix(distance[1:180, 1:180])

load("./Folding_meta/final_GWAS/munged/FI_genetic_correlation_matrix.RData")
mantel.rtest(dist(distance2, dist(correlation_matrix[-c(167), -c(167)]), nrepet = 999)
cor_cophenetic(hclust(dist(dist2)), hclust(dist(cor_2)))



