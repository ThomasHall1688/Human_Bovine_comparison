########################################################################################
##############################         Sig SNPs        #################################
########################################################################################
#This is a loop to read in all new qval files, which have our significant SNPs, subset them absed on new qval <= 0.1, and creating new .csv's, so they do not override original files. 

library(dplyr)

files = list.files(pattern="*new_qval.txt")
for(k in files){
    Genes = read.delim(k) #read file from current directory
    Genes <- subset(Genes, qval <= 0.1)
	name=paste(c(k, "_qval_0.1.csv"))
    write.csv(Genes, file = paste(c(k, "_qval_0.1.csv"), collapse=""), row.names = F) 
}
csvfilenames = list.files(pattern="*qval_0.1.csv")
for (i in 1:length(csvfilenames)) assign(csvfilenames[i], read.csv(csvfilenames[i]))
dfList <- setNames(lapply(csvfilenames, get),csvfilenames)
Combined_snps <- bind_rows(dfList, .id = "Filename")
write.csv(Combined_snps, file = "Combined_snps.csv", row.names = F) 


#Gather positions and alleles from GWAS files 

HOFR <- read.csv(file ="/home/workspace/thall/Analysis/Human_Bovine_comparison_paper/GWAS_integration/GWAS_remapping_data/HOFR_ARS.csv")
HUMAN <- read.csv(file ="/home/workspace/thall/Analysis/Human_Bovine_comparison_paper/GWAS_integration/GWAS_remapping_data/Human_GWAS_complete.csv")
colnames(HOFR)[3] <- "SNP"


bAM_MB_SNPs <- read_csv("bAM_MB_SNPs.csv")
bAM_TB_SNPs <- read_csv("bAM_TB_SNPs.csv")
hAM_TB_SNPs <- read_csv("hAM_TB_SNPs.csv")
hMDM_TB_SNPs <- read_csv("hMDM_TB_SNPs.csv")

MB_test <- merge(bAM_MB_SNPs, HOFR, by = "SNP", all.x = TRUE )
TB_test <- merge(bAM_TB_SNPs, HOFR, by = "SNP", all.x = TRUE )
hAM_test <- merge(hAM_TB_SNPs, HUMAN, by = "SNP", all.x = TRUE )
hMDM_test <- merge(hMDM_TB_SNPs, HUMAN, by = "SNP", all.x = TRUE )

write.csv(MB_test, file = "bAM_MB_significant_SNPs.csv", row.names = FALSE)
write.csv(TB_test, file = "bAM_TB_significant_SNPs.csv", row.names = FALSE)
write.csv(hAM_test, file = "hAM_TB_significant_SNPs.csv", row.names = FALSE)
write.csv(hMDM_test, file = "hMDM_TB_significant_SNPs.csv", row.names = FALSE)
