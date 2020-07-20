########################UPSETR##############################
#UpsetR takes a list of genes, with each column named after the DE list, with 1's (present gene) and 0's. I made a combined list of 
HAM_24hpi_DE_genes_FDR_0.05$HAM_24hpi <- rep(1,nrow(HAM_24hpi_DE_genes_FDR_0.05))
MDM_24hpi_DE_genes_FDR_0.05$MDM_24hpi <- rep(1,nrow(MDM_24hpi_DE_genes_FDR_0.05))
MB_24hpi_DE_genes_FDR_0.05$MB_24hpi <- rep(1,nrow(MB_24hpi_DE_genes_FDR_0.05))
TB_24hpi_DE_genes_FDR_0.05$TB_24hpi <- rep(1,nrow(TB_24hpi_DE_genes_FDR_0.05))

UpsetR_figure <- merge(Human_and_Bovine_gene_list, HAM_24hpi_DE_genes_FDR_0.05, by = "Gene_ID", all.x = TRUE)
UpsetR_figure$log2FoldChange <- NULL
UpsetR_figure$`p-value` <- NULL
UpsetR_figure$FDR <- NULL

UpsetR_figure <- merge(UpsetR_figure, MDM_24hpi_DE_genes_FDR_0.05, by = "Gene_ID", all.x = TRUE)
UpsetR_figure$log2FoldChange <- NULL
UpsetR_figure$`p-value` <- NULL
UpsetR_figure$FDR <- NULL

UpsetR_figure <- merge(UpsetR_figure, MB_24hpi_DE_genes_FDR_0.05, by = "Gene_ID", all.x = TRUE)
UpsetR_figure$log2FoldChange <- NULL
UpsetR_figure$`p-value` <- NULL
UpsetR_figure$FDR <- NULL

UpsetR_figure <- merge(UpsetR_figure, TB_24hpi_DE_genes_FDR_0.05, by = "Gene_ID", all.x = TRUE)
UpsetR_figure$log2FoldChange <- NULL
UpsetR_figure$`p-value` <- NULL
UpsetR_figure$FDR <- NULL

UpsetR_figure[is.na(UpsetR_figure)] <- 0

#simple
upset(UpsetR_figure, sets = c("hAM_TB_24hpi", "hMDM_TB_24hpi", "bAM_TB_24hpi", "bAM_MB_24hpi"), 
							sets.bar.color = c("#809fff", "#ffcc80", "#ff6666", "#8cd98c" ), order.by = "freq", empty.intersections = "on",
							matrix.color = "#3d3d5c", main.bar.color = "#52527a" )
							
#for getting avg logfc

UpsetR_figure_LFC <- merge(Human_and_Bovine_gene_list, HAM_24hpi_DE_genes_FDR_0.05, by = "Gene_ID", all.x = TRUE)
UpsetR_figure_LFC$`p-value` <- NULL


UpsetR_figure_LFC <- merge(UpsetR_figure_LFC, MDM_24hpi_DE_genes_FDR_0.05, by = "Gene_ID", all.x = TRUE)
UpsetR_figure_LFC$`p-value` <- NULL


UpsetR_figure_LFC <- merge(UpsetR_figure_LFC, MB_24hpi_DE_genes_FDR_0.05, by = "Gene_ID", all.x = TRUE)
UpsetR_figure_LFC$`p-value` <- NULL


UpsetR_figure_LFC <- merge(UpsetR_figure_LFC, TB_24hpi_DE_genes_FDR_0.05, by = "Gene_ID", all.x = TRUE)
UpsetR_figure_LFC$`p-value` <- NULL

colnames(UpsetR_figure_LFC)[2] <- "Log2FC_hAM_TB_24hpi"
colnames(UpsetR_figure_LFC)[5] <- "Log2FC_hMDM_TB_24hpi"
colnames(UpsetR_figure_LFC)[8] <- "Log2FC_bAM_MB_24hpi"
colnames(UpsetR_figure_LFC)[11] <- "Log2FC_bAM_TB_24hpi"

colnames(UpsetR_figure_LFC)[3] <- "FDR_hAM_TB_24hpi"
colnames(UpsetR_figure_LFC)[6] <- "FDR_hMDM_TB_24hpi"
colnames(UpsetR_figure_LFC)[9] <- "FDR_bAM_MB_24hpi"
colnames(UpsetR_figure_LFC)[12] <- "FDR_bAM_TB_24hpi"

UpsetR_figure_LFC <- subset(UpsetR_figure_LFC, select = c(1,2,5,8,11,3,6,9,12))

UpsetR_figure_LFC$Avg_Log2FC<- rowMeans(UpsetR_figure_LFC[,2:5], na.rm = TRUE)
UpsetR_figure_LFC$Avg_FDR<- rowMeans(UpsetR_figure_LFC[,6:9], na.rm = TRUE)

UpsetR_figure_avg <- merge(UpsetR_figure, UpsetR_figure_LFC, by = "Gene_ID", all.x = TRUE)
