######################################################################################
############################ Testing random genes ####################################
######################################################################################
#Even though the probability value gives us the % our qvalues are by chance, it is a good idea to test 
#random genes to make sure we do not get a false positive. In order to that, we create the full gene list
#This is done as the bovine work, but with the human "mart"
install.packages("biomaRt")
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

Full_human_gene_list <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'strand', 
											 'ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 
											 'entrezgene_id', 'refseq_mrna', 'description'), mart = ensembl)		
											 
#Full_human_gene_list <- Full_human_gene_list[!duplicated(Full_human_gene_list$ensembl_gene_id), ]	
Full_human_gene_list <- Full_human_gene_list[!duplicated(Full_human_gene_list$entrezgene_id), ]
#Full_bovine_gene_list <- Full_gene_list[!duplicated(Full_gene_list$?), ]	

#clean chromosome names :
Full_human_gene_list$chromosome_name <- gsub("CHR_HSCHR", "", Full_human_gene_list$chromosome_name)

Full_human_gene_list$chromosome_name <- gsub("_.*","",Full_human_gene_list$chromosome_name)

write.csv(Full_human_gene_list, file ="Full_human_gene_list.csv")

library(dplyr)
#Now we create the random set of genes for integration testing. This is done with a small loop
#setwd("C:/Users/Thomas Hall/Dropbox/New Alv_mac analysis/Analysis/Combined analysis/Gene lists of integration/Random_genes")
for(i in 1:100){
    test <- sample_n(Full_human_gene_list, 250)
    name=paste(c("Random_gene_set_", i, ".csv"), collapse="")
    write.csv(test, file=name,row.names=FALSE)
	i = i + 1
}
