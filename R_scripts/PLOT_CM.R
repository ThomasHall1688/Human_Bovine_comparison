CMplot(Bovine_orginal_CM,type="p",plot.type="c",chr.labels=paste(c(1:29),sep=""),r=0.4,cir.legend=TRUE,
        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,file="jpg", threshold = c(0.05, 0.1),
		amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","blue"),
        memo="Bovine_orginal_CM",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)
		
SNPs_bAM_MB <- bAM_MB_newQ$SNP   
CMplot(New_bAM_MB_CM,type="p",plot.type="c",chr.labels=paste(c(1:29),sep=""),r=0.4,cir.legend=TRUE,
        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,file="jpg", threshold = c(0.05, 0.1),
		amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","blue"),
        memo="bAM_MB_CM_HL",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10, highlight=SNPs_bAM_MB, highlight.col="green")
		
SNPs_bAM_TB <- bAM_TB_newQ$SNP   		
CMplot(New_bAM_TB_CM,type="p",plot.type="c",chr.labels=paste(c(1:29),sep=""),r=0.4,cir.legend=TRUE,
        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,file="jpg", threshold = c(0.05, 0.1),
		amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","blue"),
        memo="bAM_TB_CM_HL",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10, highlight=SNPs_bAM_TB, highlight.col="green")

CMplot(Human_orginal_CM,type="p",plot.type="c",chr.labels=paste(c(1:22,"X"),sep=""),r=0.4,cir.legend=TRUE,
        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,file="jpg", threshold = c(0.05, 0.1),
		amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","blue"),
        memo="Human_orginal_CM",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)	
		
SNPs_hMDM <- hMDM_newQ$SNP
CMplot(New_hMDM_TB_CM,type="p",plot.type="c",chr.labels=paste(c(1:22,"X"),sep=""),r=0.4,cir.legend=TRUE,
        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,file="jpg", threshold = c(0.05, 0.1),
		amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","blue"),
        memo="hMDM_TB_CM_HL",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10, highlight=SNPs_hMDM, highlight.col="green")	
		
SNPs_hAM <- hAM_newQ$SNP		
CMplot(New_hAM_TB_CM,type="p",plot.type="c",chr.labels=paste(c(1:22,"X"),sep=""),r=0.4,cir.legend=TRUE,
        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,file="jpg", threshold = c(0.05, 0.1),
		amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","blue"),
        memo="hAM_TB_CM_HL",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10, highlight=SNPs_hAM , highlight.col="green")	

		
		
New_hAM_TB_CM$Chromosome[is.na(New_hAM_TB_CM$Chromosome)] <- "X"
New_hMDM_TB_CM$Chromosome[is.na(New_hMDM_TB_CM$Chromosome)] <- "X"
Human_orginal_CM[1, 4] = 0.0006318753

New_hAM_TB_CM <- New_hAM_TB_CM[order(New_hAM_TB_CM$q_value),]
New_hAM_TB_CM[1, 4] = 0.0006318753


write.csv(New_bAM_MB_CM, file = "New_bAM_MB_CM.csv", row.names = FALSE)
write.csv(New_bAM_TB_CM, file = "New_bAM_TB_CM.csv", row.names = FALSE)
write.csv(New_hAM_TB_CM, file = "New_hAM_TB_CM.csv", row.names = FALSE)
write.csv(New_hMDM_TB_CM, file = "New_hMDM_TB_CM.csv", row.names = FALSE)
