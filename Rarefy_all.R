library(vegan)
taxonomy<-read.table('ASV_taxonomy_lou.txt') #ASV's with name
ASV_counts<-read.table('ASVs_counts_fixed.tsv', header=T) #ASV's with counts 


ASV_table<-t(ASV_counts) 
###################################################
library(vegan)
rarecurve(ASV_table, step = 500, xlab = "Sample Size", ylab = "Species", label = F)

ASVXXX=ASV_table[rowSums(ASV_table) > 1000-1,]
rarefied=rrarefy(ASVXXX, sample= 1000)
sum(colSums(rarefied)==0) #5225 ASV's lost; 426 samples
ASVrarefiednozero =rarefied[,colSums(rarefied)>0]
length(ASVrarefiednozero[1,]) # 8804 ASV's 365 samples remain; not all of these are from ITEX
rarecurve(ASVrarefiednozero, step = 10, xlab = "Sample Size", ylab = "Species", label = F)


#write.table(ASVrarefiednozero, file="rarefied_fungi.txt")
