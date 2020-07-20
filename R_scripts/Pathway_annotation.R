####################Pathway annotation#####################

library(org.Hs.eg.db)
MDM_TB_Pathway_genes_InnateDB <- read_csv("MDM_TB_Pathway_genes_InnateDB.csv")
hAM_TB_Pathway_genes_InnateDB <- read_csv("hAM_TB_Pathway_genes_InnateDB.csv")
bAM_TB_Pathway_genes_InnateDB <- read_csv("bAM_TB_Pathway_genes_InnateDB.csv")
bAM_MB_Pathway_genes_InnateDB <- read_csv("bAM_MB_Pathway_genes_InnateDB.csv")

PATH_bAM_MB_KEGG <- read_csv("PATH_bAM_MB_KEGG.csv")
PATH_bAM_TB_KEGG <- read_csv("PATH_bAM_TB_KEGG.csv")
PATH_hAM_TB_KEGG <- read_csv("PATH_hAM_TB_KEGG.csv")
PATH_hMDM_TB_KEGG <- read_csv("PATH_hMDM_TB_KEGG.csv")

library(org.Hs.eg.db)
hAM_TB_Pathway_genes_InnateDB$ENTREZID <- (mapIds(org.Hs.eg.db, as.character(hAM_TB_Pathway_genes_InnateDB$Ensembl), 'ENTREZID', 'ENSEMBL'))
MDM_TB_Pathway_genes_InnateDB$ENTREZID <- (mapIds(org.Hs.eg.db, as.character(MDM_TB_Pathway_genes_InnateDB$Ensembl), 'ENTREZID', 'ENSEMBL'))
PATH_hAM_TB_KEGG$ENTREZID <- (mapIds(org.Bt.eg.db, as.character(PATH_hAM_TB_KEGG$Ensembl), 'ENTREZID', 'ENSEMBL'))
PATH_hMDM_TB_KEGG$ENTREZID <- (mapIds(org.Bt.eg.db, as.character(PATH_hMDM_TB_KEGG$Ensembl), 'ENTREZID', 'ENSEMBL'))


library(org.Bt.eg.db)
bAM_MB_Pathway_genes_InnateDB$ENTREZID <- (mapIds(org.Bt.eg.db, as.character(bAM_MB_Pathway_genes_InnateDB$Ensembl), 'ENTREZID', 'ENSEMBL'))
bAM_TB_Pathway_genes_InnateDB$ENTREZID <- (mapIds(org.Bt.eg.db, as.character(bAM_TB_Pathway_genes_InnateDB$Ensembl), 'ENTREZID', 'ENSEMBL'))
PATH_bAM_MB_KEGG$ENTREZID <- (mapIds(org.Bt.eg.db, as.character(PATH_bAM_MB_KEGG$Ensembl), 'ENTREZID', 'ENSEMBL'))
PATH_bAM_TB_KEGG$ENTREZID <- (mapIds(org.Bt.eg.db, as.character(PATH_bAM_TB_KEGG$Ensembl), 'ENTREZID', 'ENSEMBL'))


write.csv(hAM_TB_Pathway_genes_InnateDB, file = "hAM_TB_Pathway_genes_InnateDB.csv", row.names = FALSE)
write.csv(MDM_TB_Pathway_genes_InnateDB, file = "MDM_TB_Pathway_genes_InnateDB.csv", row.names = FALSE)
write.csv(bAM_MB_Pathway_genes_InnateDB, file = "bAM_MB_Pathway_genes_InnateDB.csv", row.names = FALSE)
write.csv(bAM_TB_Pathway_genes_InnateDB, file = "bAM_TB_Pathway_genes_InnateDB.csv", row.names = FALSE)

write.csv(PATH_bAM_MB_KEGG, file = "PATH_bAM_MB_KEGG.csv", row.names = FALSE)
write.csv(PATH_bAM_TB_KEGG, file = "PATH_bAM_TB_KEGG.csv", row.names = FALSE)
write.csv(PATH_hAM_TB_KEGG, file = "PATH_hAM_TB_KEGG.csv", row.names = FALSE)
write.csv(PATH_hMDM_TB_KEGG, file = "PATH_hMDM_TB_KEGG.csv", row.names = FALSE)



MDM_TB_Pathway_genes_InnateDB <- read_csv("MDM_TB_Pathway_genes_InnateDB.csv")
hAM_TB_Pathway_genes_InnateDB <- read_csv("hAM_TB_Pathway_genes_InnateDB.csv")
bAM_TB_Pathway_genes_InnateDB <- read_csv("bAM_TB_Pathway_genes_InnateDB.csv")
bAM_MB_Pathway_genes_InnateDB <- read_csv("bAM_MB_Pathway_genes_InnateDB.csv")

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bensembl <- useMart("ensembl", dataset = "btaurus_gene_ensembl")


MDM_TB_Pathway_genes_InnateDB_annotation <- getBM(attributes = c('entrezgene_id', 'description', 'chromosome_name', 'start_position', 'end_position', 'strand'),  
                              filter = 'entrezgene_id',
                              values = MDM_TB_Pathway_genes_InnateDB$ENTREZID, 
                              mart = ensembl)

colnames(MDM_TB_Pathway_genes_InnateDB_annotation)[1] <- "ENTREZID"
MDM_TB_Pathway_genes_InnateDB_annotation <- MDM_TB_Pathway_genes_InnateDB_annotation[!duplicated(MDM_TB_Pathway_genes_InnateDB_annotation$ENTREZID), ]
MDM_TB_Pathway_genes_InnateDB_annotation <- merge(MDM_TB_Pathway_genes_InnateDB, MDM_TB_Pathway_genes_InnateDB_annotation, by = "ENTREZID", all.x = TRUE)


hAM_TB_Pathway_genes_InnateDB_annotation <- getBM(attributes = c('entrezgene_id', 'description', 'chromosome_name', 'start_position', 'end_position', 'strand'),  
                              filter = 'entrezgene_id',
                              values = hAM_TB_Pathway_genes_InnateDB$ENTREZID, 
                              mart = ensembl)

colnames(hAM_TB_Pathway_genes_InnateDB_annotation)[1] <- "ENTREZID"
hAM_TB_Pathway_genes_InnateDB_annotation <- hAM_TB_Pathway_genes_InnateDB_annotation[!duplicated(hAM_TB_Pathway_genes_InnateDB_annotation$ENTREZID), ]
hAM_TB_Pathway_genes_InnateDB_annotation <- merge(hAM_TB_Pathway_genes_InnateDB, hAM_TB_Pathway_genes_InnateDB_annotation, by = "ENTREZID", all.x = TRUE)


bAM_TB_Pathway_genes_InnateDB_annotation <- getBM(attributes = c('entrezgene_id', 'description', 'chromosome_name', 'start_position', 'end_position', 'strand'),  
                              filter = 'entrezgene_id',
                              values = bAM_TB_Pathway_genes_InnateDB$ENTREZID, 
                              mart = bensembl)

colnames(bAM_TB_Pathway_genes_InnateDB_annotation)[1] <- "ENTREZID"
bAM_TB_Pathway_genes_InnateDB_annotation <- bAM_TB_Pathway_genes_InnateDB_annotation[!duplicated(bAM_TB_Pathway_genes_InnateDB_annotation$ENTREZID), ]
bAM_TB_Pathway_genes_InnateDB_annotation <- merge(bAM_TB_Pathway_genes_InnateDB, bAM_TB_Pathway_genes_InnateDB_annotation, by = "ENTREZID", all.x = TRUE)


bAM_MB_Pathway_genes_InnateDB_annotation <- getBM(attributes = c('entrezgene_id', 'description', 'chromosome_name', 'start_position', 'end_position', 'strand'),  
                              filter = 'entrezgene_id',
                              values = bAM_MB_Pathway_genes_InnateDB$ENTREZID, 
                              mart = bensembl)

colnames(bAM_MB_Pathway_genes_InnateDB_annotation)[1] <- "ENTREZID"
bAM_MB_Pathway_genes_InnateDB_annotation <- bAM_MB_Pathway_genes_InnateDB_annotation[!duplicated(bAM_MB_Pathway_genes_InnateDB_annotation$ENTREZID), ]
bAM_MB_Pathway_genes_InnateDB_annotation <- merge(bAM_MB_Pathway_genes_InnateDB, bAM_MB_Pathway_genes_InnateDB_annotation, by = "ENTREZID", all.x = TRUE)

write.csv(hAM_TB_Pathway_genes_InnateDB_annotation, file = "hAM_TB_Pathway_genes_InnateDB.csv", row.names = FALSE)
write.csv(MDM_TB_Pathway_genes_InnateDB_annotation, file = "MDM_TB_Pathway_genes_InnateDB.csv", row.names = FALSE)
write.csv(bAM_MB_Pathway_genes_InnateDB_annotation, file = "bAM_MB_Pathway_genes_InnateDB.csv", row.names = FALSE)
write.csv(bAM_TB_Pathway_genes_InnateDB_annotation, file = "bAM_TB_Pathway_genes_InnateDB.csv", row.names = FALSE)
