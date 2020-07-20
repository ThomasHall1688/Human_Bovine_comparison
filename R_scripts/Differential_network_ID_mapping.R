########################################################################################
###############################         DEN          ###################################
########################################################################################

#Id's need to be converted in order to import them into both innateDB and Cystoscape via Jactive modules 

library(org.Hs.eg.db)
columns(org.Hs.eg.db)
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bensembl <- useMart("ensembl", dataset = "btaurus_gene_ensembl")

#Read in data 
setwd("C:/Users/Thomas Hall/Dropbox/Human_bovine_comparison/Analysis/DEN analysis/Network input")

GeneCards_2.5_Humie <- read_csv("GeneCards_2.5(Humie).csv")
HAM_24hpi_DE_genes_FDR_0.05 <- read_csv("HAM_24hpi_DE_genes_FDR_0.05.csv")
MB_24hpi_DE_genes_FDR_0.05 <- read_csv("MB_24hpi_DE_genes_FDR_0.05.csv")
MDM_24hpi_DE_genes_FDR_0.05 <- read_csv("MDM_24hpi_DE_genes_FDR_0.05.csv")
TB_24hpi_DE_genes_FDR_0.05 <- read_csv("TB_24hpi_DE_genes_FDR_0.05.csv")

#Genecards output 
GeneCards_2.5_Humie$HS_ENTREZ<- (mapIds(org.Hs.eg.db, GeneCards_2.5_Humie$`Gene Symbol`,'ENTREZID', 'SYMBOL'))
GeneCards_2.5_Humie$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, GeneCards_2.5_Humie$`Gene Symbol`, 'ENSEMBL','SYMBOL'))

HAM_24hpi_DE_genes_FDR_0.05$HS_ENTREZ<- (mapIds(org.Hs.eg.db, HAM_24hpi_DE_genes_FDR_0.05$Gene_ID,'ENTREZID', 'SYMBOL'))
HAM_24hpi_DE_genes_FDR_0.05$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, HAM_24hpi_DE_genes_FDR_0.05$Gene_ID, 'ENSEMBL','SYMBOL'))

MDM_24hpi_DE_genes_FDR_0.05$HS_ENTREZ<- (mapIds(org.Hs.eg.db, MDM_24hpi_DE_genes_FDR_0.05$Gene_ID,'ENTREZID', 'SYMBOL'))
MDM_24hpi_DE_genes_FDR_0.05$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, MDM_24hpi_DE_genes_FDR_0.05$Gene_ID, 'ENSEMBL','SYMBOL'))

MB_24hpi_DE_genes_FDR_0.05$HS_ENTREZ<- (mapIds(org.Hs.eg.db, MB_24hpi_DE_genes_FDR_0.05$Gene_ID,'ENTREZID', 'SYMBOL'))
MB_24hpi_DE_genes_FDR_0.05$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, MB_24hpi_DE_genes_FDR_0.05$Gene_ID, 'ENSEMBL','SYMBOL'))

TB_24hpi_DE_genes_FDR_0.05$HS_ENTREZ<- (mapIds(org.Hs.eg.db, TB_24hpi_DE_genes_FDR_0.05$Gene_ID,'ENTREZID', 'SYMBOL'))
TB_24hpi_DE_genes_FDR_0.05$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, TB_24hpi_DE_genes_FDR_0.05$Gene_ID, 'ENSEMBL','SYMBOL'))


write.csv(GeneCards_2.5_Humie, file = "GeneCards_2.5_Humie_Annotated.csv", row.names = FALSE)
write.csv(HAM_24hpi_DE_genes_FDR_0.05, file = "HAM_24hpi_DE_genes_FDR_0.05_Annotated.csv", row.names = FALSE)
write.csv(MDM_24hpi_DE_genes_FDR_0.05, file = "MDM_24hpi_DE_genes_FDR_0.05_Annotated.csv", row.names = FALSE)
write.csv(MB_24hpi_DE_genes_FDR_0.05, file = "MB_24hpi_DE_genes_FDR_0.05_Annotated.csv", row.names = FALSE)
write.csv(TB_24hpi_DE_genes_FDR_0.05, file = "TB_24hpi_DE_genes_FDR_0.05_Annotated.csv", row.names = FALSE)


#As usual, entrex ids and everything else is missing from a lot of genes. Doing this manually will take forever, so lets try to fix this first with R as much as possible./ 
Missing_gene_hAM <- read_csv("Missing_gene_hAM.csv")

Missing_gene_hAM_annotation <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id','description' ),  
                              filter = 'external_gene_name',
                              values = Missing_gene_hAM$Gene_ID, 
                              mart = ensembl)


#testing for complete annotate
test <- read_csv("HAM_24hpi_DE_genes_FDR_0.05.csv")
test$HS_ENTREZ <- (mapIds(org.Hs.eg.db, test$Gene_ID, 'ENTREZID','SYMBOL'))
test$HS_SYMBOL<- (mapIds(org.Hs.eg.db, test$Gene_ID,'ALIAS', 'SYMBOL'))
test$HS_ENTREZ_2 <- (mapIds(org.Hs.eg.db, test$HS_SYMBOL, 'ENTREZID','ALIAS'))
test$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, test$Gene_ID, 'ENSEMBL','SYMBOL'))
test$HS_GENENAME <- (mapIds(org.Hs.eg.db, test$HS_ENTREZ, 'GENENAME','ENTREZID'))
test$HS_REFSEQ <- (mapIds(org.Hs.eg.db, test$Gene_ID, "REFSEQ",'SYMBOL'))
test$HS_ALIAS <- (mapIds(org.Hs.eg.db, test$Gene_ID, "ALIAS",'SYMBOL'))
test$HS_ENSEMBL_2 <- (mapIds(org.Hs.eg.db, test$HS_ALIAS, 'ENSEMBL','ALIAS'))
test$HS_ENSEMBL_3 <- (mapIds(org.Hs.eg.db, test$HS_ENTREZ, 'ENSEMBL','ENTREZID'))


test_annotation <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_id', 'external_gene_name', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene_id',
                              values = test$HS_ENTREZ, 
                              mart = ensembl)
							  
							  
							  
							  hgnc_symbol
							  
							
####################DEN annotation#####################

DEN_hMDM_24hpi$HS_SYMBOL<- (mapIds(org.Hs.eg.db, as.character(DEN_hMDM_24hpi$HS_ENTREZ), 'SYMBOL', 'ENTREZID'))
DEN_hMDM_24hpi$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, as.character(DEN_hMDM_24hpi$HS_ENTREZ), 'ENSEMBL','ENTREZID'))

DEN_hAM_24hpi$HS_SYMBOL<- (mapIds(org.Hs.eg.db, as.character(DEN_hAM_24hpi$HS_ENTREZ), 'SYMBOL', 'ENTREZID'))
DEN_hAM_24hpi$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, as.character(DEN_hAM_24hpi$HS_ENTREZ), 'ENSEMBL','ENTREZID'))

DEN_MB_24hpi$HS_SYMBOL<- (mapIds(org.Hs.eg.db, as.character(DEN_MB_24hpi$HS_ENTREZ), 'SYMBOL', 'ENTREZID'))
DEN_MB_24hpi$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, as.character(DEN_MB_24hpi$HS_ENTREZ), 'ENSEMBL','ENTREZID'))

DEN_TB_24hpi$HS_SYMBOL<- (mapIds(org.Hs.eg.db, as.character(DEN_TB_24hpi$HS_ENTREZ), 'SYMBOL', 'ENTREZID'))
DEN_TB_24hpi$HS_ENSEMBL <- (mapIds(org.Hs.eg.db, as.character(DEN_TB_24hpi$HS_ENTREZ), 'ENSEMBL','ENTREZID'))


DEN_hMDM_24hpi_annotation <- getBM(attributes = c('entrezgene_id','description', 'chromosome_name', 'start_position', 'end_position', 'strand'),  
                              filter = 'entrezgene_id',
                              values = DEN_hMDM_24hpi$HS_ENTREZ, 
                              mart = ensembl)

colnames(DEN_hMDM_24hpi_annotation)[1] <- "HS_ENTREZ"
DEN_hMDM_24hpi_annotation <- DEN_hMDM_24hpi_annotation[!duplicated(DEN_hMDM_24hpi_annotation$HS_ENTREZ), ]
DEN_hMDM_Annotated <- merge(DEN_hMDM_24hpi, DEN_hMDM_24hpi_annotation, by = "HS_ENTREZ", all.x = TRUE)



DEN_hAM_24hpi_annotation <- getBM(attributes = c('entrezgene_id','description', 'chromosome_name', 'start_position', 'end_position', 'strand'),  
                              filter = 'entrezgene_id',
                              values = DEN_hAM_24hpi$HS_ENTREZ, 
                              mart = ensembl)

colnames(DEN_hAM_24hpi_annotation)[1] <- "HS_ENTREZ"
DEN_hAM_24hpi_annotation <- DEN_hAM_24hpi_annotation[!duplicated(DEN_hAM_24hpi_annotation$HS_ENTREZ), ]
DEN_hAM_Annotated <- merge(DEN_hAM_24hpi, DEN_hAM_24hpi_annotation, by = "HS_ENTREZ", all.x = TRUE)



DEN_TB_24hpi_annotation <- getBM(attributes = c('external_gene_name', 'entrezgene_id','description', 'chromosome_name', 'start_position', 'end_position', 'strand'),  
                              filter = 'external_gene_name',
                              values = DEN_TB_24hpi$HS_SYMBOL, 
                              mart = bensembl)

colnames(DEN_TB_24hpi_annotation)[1] <- "HS_SYMBOL"
DEN_TB_24hpi_annotation <- DEN_TB_24hpi_annotation[!duplicated(DEN_TB_24hpi_annotation$HS_SYMBOL), ]
DEN_TB_Annotated <- merge(DEN_TB_24hpi, DEN_TB_24hpi_annotation, by = "HS_SYMBOL", all.x = TRUE)




DEN_MB_24hpi_annotation <- getBM(attributes = c('external_gene_name', 'entrezgene_id','description', 'chromosome_name', 'start_position', 'end_position', 'strand'),  
                              filter = 'external_gene_name',
                              values = DEN_MB_24hpi$HS_SYMBOL, 
                              mart = bensembl)

colnames(DEN_MB_24hpi_annotation)[1] <- "HS_SYMBOL"
DEN_MB_24hpi_annotation <- DEN_MB_24hpi_annotation[!duplicated(DEN_MB_24hpi_annotation$HS_SYMBOL), ]
DEN_MB_Annotated <- merge(DEN_MB_24hpi, DEN_MB_24hpi_annotation, by = "HS_SYMBOL", all.x = TRUE)

write.csv(DEN_MB_Annotated, file = "DEN_MB_Annotated.csv", row.names = FALSE)
write.csv(DEN_TB_Annotated, file = "DEN_TB_Annotated.csv", row.names = FALSE)
write.csv(DEN_hAM_Annotated, file = "DEN_hAM_Annotated.csv", row.names = FALSE)
write.csv(DEN_hMDM_Annotated, file = "DEN_hMDM_Annotated.csv", row.names = FALSE)
