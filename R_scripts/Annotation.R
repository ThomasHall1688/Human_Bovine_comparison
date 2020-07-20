########################################################################################
#############################         Annotation        ################################
########################################################################################

#Gene coordinates, as well as other information is needed for donwstream analysis and integration. 
#This pipeline is convulted, mainly because the annotation on most databases is still out of date. It will be refined in the future i.e. only the biomaRt step.


#first check available ID's 
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

#Then append ID's based on row names 
#Convert one list of genes to another id type. First '' is what you want, second '' is what you have i.e the first line converst Entrez to Symbol

HAM_24hpi_DE_top_500$ENTREZ<- (mapIds(org.Hs.eg.db, HAM_24hpi_DE_top_500$`Gene ID`,'ENTREZID', 'SYMBOL'))
HAM_24hpi_DE_top_500$ENSEMBL <- (mapIds(org.Hs.eg.db, HAM_24hpi_DE_top_500$ENTREZ, 'ENSEMBL','ENTREZID'))
HAM_24hpi_DE_top_500$GENENAME <- (mapIds(org.Hs.eg.db, HAM_24hpi_DE_top_500$ENTREZ, 'GENENAME','ENTREZID'))

MDM_24hpi_DE_top_500$ENTREZ<- (mapIds(org.Hs.eg.db, MDM_24hpi_DE_top_500$`Gene ID`,'ENTREZID', 'SYMBOL'))
MDM_24hpi_DE_top_500$ENSEMBL <- (mapIds(org.Hs.eg.db, MDM_24hpi_DE_top_500$ENTREZ, 'ENSEMBL','ENTREZID'))
MDM_24hpi_DE_top_500$GENENAME <- (mapIds(org.Hs.eg.db, MDM_24hpi_DE_top_500$ENTREZ, 'GENENAME','ENTREZID'))



library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

HAM_24hpi_DE_top_500_annotation <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = MB2_FDR_0.05_genes$`rownames(subset(MB2_res, padj < 0.05))`, 
                              mart = ensembl)

MDM_24hpi_DE_top_500_annotation<- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ), 
                              filter = 'ensembl_gene_id',
                              values = MDM_24hpi_DE_top_500$ENSEMBL, 
                              mart = ensembl)							  
							 

							 
#This didnt really work so I just did it on the website and it was fine. 
HAM_24hpi_DE_top_500_annotated <- merge(HAM_24hpi_DE_top_500, mart_export_HAM, by = "Gene_name", all.x = TRUE)	
MDM_24hpi_DE_top_500_annotated <- merge(MDM_24hpi_DE_top_500, mart_export_MDM, by = "Gene_name", all.x = TRUE)	 
				
write.csv(HAM_24hpi_DE_top_500_annotated  , file="HAM_24hpi_DE_top_500_annotated.csv", row.names  = FALSE)
write.csv(MDM_24hpi_DE_top_500_annotated  , file="MDM_24hpi_DE_top_500_annotated.csv", row.names  = FALSE)				
#Obviously, it missed some, which I did manually. These were (HAM)

MARCH1 (was a date)
ADCK3					
C10orf128
C15orf52
C19orf59
C5orf20
EMR3
FAM105A
FAM129B
FAM198B
FAM26F
FYB
IL8
LINC00515
LOC283070
LOC285972
LSMD1
MB21D1
METTL21B
MGC12916
MICALCL
MLLT4
MST4
O3FAR1
PQLC3
SGK110
TMEM194A
ZAK

#which convert to 

MARCH1
COQ8A
TMEM273
CCDC9B
MCEMP1
DCANP1
ADGRE3
OTULINL 
NIBAN2 
FAM198B-AS1 
CALHM6
FYB1
IL8
LINC00515
AIG1
IRF7
NAA38
CGAS
EEF1AKMT3 
MGC12916
MICAL2
AFDN-DT
MST4
FFAR4
SLC66A3
SBK3
NEMP1
MAP3K20


#MDM_24hpi_DE_top_500

AMICA1
C14orf159
C21orf7
C5orf4
C7orf29
C7orf41
CXorf57
EMR1
EMR3
FAIM3
FAM190A
FAM26F
FCGR1C
GPER
GPR126
IL8
KIAA0247
KIAA1147
KIAA1199
KIAA1462
LINC00341
LOC100289511
LOC100499227
LOC440896
LOC728431
MGC12916
PCNX
RARRES3
SEPP1

#Get FDR cutoff for the 4 datasets

HAM_24hpi_DE_genes_FDR_0.05  <- subset(HAM_24hpi_DE_genes, FDR <= 0.05)
MDM_24hpi_DE_genes_FDR_0.05  <- subset(MDM_24hpi_DE_genes, FDR <= 0.05)			
MB_24hpi_DE_genes_FDR_0.05 <- subset(MB_24hpi_DE_genes, FDR <= 0.05)
TB_24hpi_DE_genes_FDR_0.05  <- subset(TB_24hpi_DE_genes, FDR <= 0.05)

write.csv(HAM_24hpi_DE_genes_FDR_0.05, file = "HAM_24hpi_DE_genes_FDR_0.05.csv", row.names = FALSE)
write.csv(MDM_24hpi_DE_genes_FDR_0.05, file = "MDM_24hpi_DE_genes_FDR_0.05.csv", row.names = FALSE)
write.csv(MB_24hpi_DE_genes_FDR_0.05, file = "MB_24hpi_DE_genes_FDR_0.05.csv", row.names = FALSE)
write.csv(TB_24hpi_DE_genes_FDR_0.05, file = "TB_24hpi_DE_genes_FDR_0.05.csv", row.names = FALSE)


#LFC Annotation 
DE_hAM_TB_LFC$ENTREZ<- (mapIds(org.Hs.eg.db, DE_ham_TB_LFC$$`Gene ID`,'ENTREZID', 'SYMBOL'))
DE_hMDM_TB_LFC$ENTREZ<- (mapIds(org.Hs.eg.db, DE_hMDM_TB_LFC$$`Gene ID`,'ENTREZID', 'SYMBOL'))




DE_bAM_TB_LFC_annotation

DE_hAM_TB_LFC_annotation

DE_hMDM_TB_LFC_annotation



library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bensembl <- useMart("ensembl", dataset = "btaurus_gene_ensembl")

DE_bAM_TB_LFC_annotation <- getBM(attributes = c('entrezgene_id', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene_id',
                              values = DE_bAM_TB_LFC$row, 
                              mart = bensembl)


DE_hAM_TB_LFC_annotation <- getBM(attributes = c('entrezgene_id', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene_id',
                              values = DE_hAM_TB_LFC$ENTREZ, 
                              mart = ensembl)

DE_hMDM_TB_LFC_annotation <- getBM(attributes = c('entrezgene_id', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene_id',
                              values = DE_hMDM_TB_LFC$ENTREZ, 
                              mart = ensembl)							  

							  
colnames(DE_bAM_TB_LFC_annotation)[1] <- "row"
DE_bAM_TB_LFC_annotation <- DE_bAM_TB_LFC_annotation[!duplicated(DE_bAM_TB_LFC_annotation$row), ]
DE_bAM_TB_LFC_Annotated <- merge(DE_bAM_TB_LFC, DE_bAM_TB_LFC_annotation, by = "row", all.x = TRUE)		

colnames(DE_hAM_TB_LFC_annotation)[1] <- "ENTREZ"
DE_hAM_TB_LFC_annotation <- DE_hAM_TB_LFC_annotation[!duplicated(DE_hAM_TB_LFC_annotation$ENTREZ), ]
DE_hAM_TB_LFC_Annotated <- merge(DE_hAM_TB_LFC, DE_hAM_TB_LFC_annotation, by = "ENTREZ", all.x = TRUE)		

colnames(DE_hMDM_TB_LFC_annotation)[1] <- "ENTREZ"
DE_hMDM_TB_LFC_annotation <- DE_hMDM_TB_LFC_annotation[!duplicated(DE_hMDM_TB_LFC_annotation$ENTREZ), ]
DE_hMDM_TB_LFC_Annotated <- merge(DE_hMDM_TB_LFC, DE_hMDM_TB_LFC_annotation, by = "ENTREZ", all.x = TRUE)		

write.csv(DE_bAM_TB_LFC_Annotated, file = "DE_bAM_TB_LFC_Annotated.csv", row.names = FALSE)
write.csv(DE_hAM_TB_LFC_Annotated, file = "DE_hAM_TB_LFC_Annotated.csv", row.names = FALSE)
write.csv(DE_hMDM_TB_LFC_Annotated, file = "DE_hMDM_TB_LFC_Annotated.csv", row.names = FALSE)
