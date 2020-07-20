##########################scatter and dot plots#####################################
#dot plots 
Common_DE_FDRgenes_all <- merge(Common_DE_FDRgenes_all, UpsetR_figure_LFC, by = "Gene_ID", all.x = TRUE)

Common_DE_FDRgenes_LFC <- subset(Common_DE_FDRgenes_all, select = c(1,2,3,4,5))
rownames(Common_DE_FDRgenes_LFC) <- Common_DE_FDRgenes_LFC[,1]
Common_DE_FDRgenes_LFC[,1] <- NULL

Common_DE_FDRgenes_LFC <- melt(Common_DE_FDRgenes_LFC)


#Dot plot#
p = ggplot(Common_DE_FDRgenes_LFC_melt, aes(x=Common_DE_FDRgenes_LFC_melt$`Sample at 24hpi`, y=Common_DE_FDRgenes_LFC_melt$`Log2 fold change`, fill = factor(colour))) +
    geom_dotplot(binaxis='y', stackdir='center', method = "dotdensity", stackratio = 0.7, dotsize = 2.5,
    binwidth = 0.2, show.legend = FALSE) +
    geom_label_repel(aes(label = Common_DE_FDRgenes_LFC_melt$Genes), box.padding   = 1.3, point.padding = 0.5, label.size = 0.5, label.padding = 0.35,
	show.legend = FALSE, na.rm = TRUE, colour = "#262626", label.r = 0.35) +
	geom_hline(yintercept=0, linetype="dashed", color = "#000000")
p + labs(x = "Group at 24hpi", y = "Log2 fold change")  +
    scale_fill_manual(values=setNames(Common_DE_FDRgenes_LFC_melt$Code, Common_DE_FDRgenes_LFC_melt$colour)) +
    scale_colour_manual(values=setNames(Common_DE_FDRgenes_LFC_melt$Code, Common_DE_FDRgenes_LFC_melt$colour)) +
	theme(axis.text.x = element_text(color="#262626", size=14),
          axis.text.y = element_text(color="#262626", size=14),
		  axis.title = element_text(color="#262626", size=18,face="bold"))

#	 geom_text(aes(label=Common_DE_FDRgenes_LFC_melt$Genes),hjust=0, vjust=0) 


#four group histo

#Create vectors with desirable variables for graphing
Time.vec<-c(rep("hAM TB",2),rep("hMDM TB",2),
            rep("bAM TB",2),rep("bAM MB",2))
Variable.vec<-c("Up","Down","Up","Down","Up","Down","Up","Down")
Variable.condition.vec<-c("hAM_up","hAM_down","hMDM_up","hMDM_down","bAM_up",
                          "bAM_down","bAM_MB_up","bAM_MB_down")
Code.vec <- c("#8cd98c", "#c6ecc6", "#ff6666", "#ffb3b3", "#809fff", "#ccd9ff", "#ffcc80", "#ffebcc")
						 
#Going to count how many genes are up and down so,
#set up empty vector with blank entries to store the gene counts
values.vec<-rep("x",8)

#here, we need the LF2C, so we use the LFC datasets



HAM_24hpi_DE_genes_FDR_0.05_LFC <- as.vector(HAM_24hpi_DE_genes_FDR_0.05["log2FoldChange"])
values.vec[1] <-sum(HAM_24hpi_DE_genes_FDR_0.05_LFC > 0)
values.vec[2] <-sum(HAM_24hpi_DE_genes_FDR_0.05_LFC < 0)

MDM_24hpi_DE_genes_FDR_0.05_LFC <- as.vector(MDM_24hpi_DE_genes_FDR_0.05["log2FoldChange"])
values.vec[3] <-sum(MDM_24hpi_DE_genes_FDR_0.05_LFC > 0)
values.vec[4] <-sum(MDM_24hpi_DE_genes_FDR_0.05_LFC < 0)

TB_24hpi_DE_genes_FDR_0.05_LFC <- as.vector(TB_24hpi_DE_genes_FDR_0.05["log2FoldChange"])
values.vec[5] <-sum(TB_24hpi_DE_genes_FDR_0.05_LFC > 0)
values.vec[6] <-sum(TB_24hpi_DE_genes_FDR_0.05_LFC < 0)

MB_24hpi_DE_genes_FDR_0.05_LFC <- as.vector(MB_24hpi_DE_genes_FDR_0.05["log2FoldChange"])
values.vec[7] <-sum(MB_24hpi_DE_genes_FDR_0.05_LFC > 0)
values.vec[8] <-sum(MB_24hpi_DE_genes_FDR_0.05_LFC < 0)


bar_data.raw<-data.frame(a=character(),b=character(),c=numeric(), d=character())
bar_data<-rbind(bar_data.raw, data.frame(a=Time.vec, b=Variable.vec, d=as.numeric(values.vec), e=Variable.condition.vec, f=Code.vec))
colnames(bar_data)<-c("Group","Variable","Genes","Variable.condition", "Codes")

#Make custom labels for legend of graph to include both italicised and plain text
label_5<-expression(paste(italic(""),"Increased expression"))
label_6<-expression(paste(italic(""),"Decreased expression"))

q<-ggplot(bar_data, aes(Group), ylim(-2500:2500)) + 
  geom_bar(data = subset(bar_data, Variable == "Up"), 
           aes(y = Genes, fill = Variable.condition), stat = "identity", position = "dodge",colour="black",size=0.4,width=0.7) +
#  scale_fill_manual(values=c("#8cd98c", "#c6ecc6", "#ff6666", "#ffb3b3", "#809fff", "#ccd9ff", "#ffcc80", "#ffebcc"), 
  scale_fill_manual(values=c("#ffebcc", "#ccd9ff", "#809fff", "#ffcc80", "#c6ecc6", "#8cd98c", "#ffb3b3", "#ff6666"),
                    name=" ",
                    breaks=c("delta_up", "delta_down"), #define the
                    #breaks so that you can relabel
                    labels=c(label_5, label_6)) +
  geom_bar(data = subset(bar_data, Variable == "Down"), #colours are bovis up, tb up, bovis down, tb down
           aes(y = -Genes, fill = Variable.condition), stat = "identity", position = "dodge",colour="black",size=0.4,width=0.7) + 
  geom_hline(yintercept = 0,colour = "black") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text.align = 0)  #aligning the legend labels to legend boxes

#Defining custom subsets for adding stars to the graph
x_up<-(data = subset(bar_data, Variable == "Up"))
x_down<-(data = subset(bar_data, Variable == "Down"))
x_delta_2_up<-subset(x_up[1,])
x_delta_2_down<-(x_down[1,])
x_delta_6_up<-(x_up[2,])
x_delta_6_down<-(x_down[2,])
x_delta_24_up<-(x_up[3,])
x_delta_24_down<-(x_down[3,])
x_delta_48_up<-(x_up[4,])
x_delta_48_down<-(x_down[4,])

q +
  geom_text(data = subset(bar_data, Variable == "Up"), 
            aes(Group, Genes, group=Variable.condition, label=Genes),
            position = position_dodge(width=0.9), vjust = -0.25, size=3.5) +
  geom_text(data = subset(bar_data, Variable == "Down"), 
            aes(Group, -Genes, group=Variable.condition, label=Genes),
            position = position_dodge(width=0.9), vjust = 1.25, size=3.5) +
  coord_cartesian(ylim = c(-2100, 2100)) +
  scale_x_discrete(name="Group at 24hpi",
                   labels=c("bAM MB","bAM TB","hAM TB","hMDM TB")) + #getting rid of the 0 in 02hr and 06hr
  scale_y_continuous("Number of differentially expressed genes") +
  theme(legend.text=element_text(size=9),legend.key.size=unit(0.4,"cm")) + #changing size of legend
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.title.y=element_text(size=11)) +
  theme(aspect.ratio = 1.3) +
  theme(legend.position="bottom", legend.box = "horizontal") + #horizontal legend at bottom of graph
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.2)) +
  annotate("text", x = 1, y = -2200, size=3.5, label = "FDR < 0.05")
