################################# SNP bar #####################################
bAM_MB <- read_excel("bAM_MB.xlsx")
View(bAM_MB)
library(readxl)
bAM_TB <- read_excel("bAM_TB.xlsx")
View(bAM_TB)
library(readxl)
hAM_TB <- read_excel("hAM_TB.xlsx")
View(hAM_TB)
library(readxl)
hMDM_TB <- read_excel("hMDM_TB.xlsx")
View(hMDM_TB)



my_title <- expression(paste("bAM infected with ", italic("Mycobacterium bovis"))) 
bAM_MB_plot <- ggplot(bAM_MB, aes(fill=Log2FC, y=SNPs, x=Gene)) + 
    geom_bar(position="dodge", stat="identity", width = 0.75) +
    coord_flip() +
    scale_fill_gradient(low = "#4d94ff", high = "#ff3333", na.value = "#81DA84") +
    ggtitle(my_title) +
    #geom_text(aes(label=Description), position=position_stack(vjust=0.75), colour="#404040") + 
    theme(
        axis.title.x = element_text(color="#404040", size=14, face="bold"),
        axis.title.y = element_text(color="#404040", size=14, face="bold"),
        axis.text.y = element_text(color="#404040",  size=12, face = "italic")) +
		ylim(0, 50)
plot(bAM_MB_plot)

	
my_title <- expression(paste("bAM infected with ", italic("Mycobacterium tuberculosis"))) 
bAM_TB_plot <- ggplot(bAM_TB, aes(fill=Log2FC, y=SNPs, x=Gene)) + 
    geom_bar(position="dodge", stat="identity", width = 0.75) +
    coord_flip() +
    scale_fill_gradient(low = "#4d94ff", high = "#ff3333", na.value = "#81DA84") +
    ggtitle(my_title) +
    #geom_text(aes(label=Description), colour="#404040", position=position_stack(vjust=0.5)) + 
    theme(
        axis.title.x = element_text(color="#404040", size=14, face="bold"),
        axis.title.y = element_text(color="#404040", size=14, face="bold"),
        axis.text.y = element_text(color="#404040",  size=12, face = "italic")) +
		ylim(0, 50)
plot(bAM_TB_plot)

my_title <- expression(paste("hAM infected with ", italic("Mycobacterium tuberculosis")))
hAM_TB_plot <- ggplot(hAM_TB, aes(fill=Log2FC, y=SNPs, x=Gene)) + 
    geom_bar(position="dodge", stat="identity", width = 0.75) +
    coord_flip() +
    scale_fill_gradient(low = "#4d94ff", high = "#ff3333", na.value = "#81DA84") +
    ggtitle(my_title) +
   #geom_text(aes(label=Description), colour="#404040", position=position_stack(vjust=0.5)) + 
    theme(
        axis.title.x = element_text(color="#404040", size=14, face="bold"),
        axis.title.y = element_text(color="#404040", size=14, face="bold"),
        axis.text.y = element_text(color="#404040",  size=12, face = "italic")) +
		ylim(0, 50)
plot(hAM_TB_plot)

my_title <- expression(paste("hMDM infected with ", italic("Mycobacterium tuberculosis"))) 
hMDM_TB_plot <- ggplot(hMDM_TB, aes(fill=Log2FC, y=SNPs, x=Gene)) + 
    geom_bar(position="dodge", stat="identity", width = 0.75) +
    coord_flip() +
    scale_fill_gradient(low = "#4d94ff", high = "#ff3333", na.value = "#81DA84") +
    ggtitle(my_title) +
    #geom_text(aes(label=Description), colour="#404040", position=position_stack(vjust=0.5), vjust="inward") + 
    theme(
        axis.title.x = element_text(color="#404040", size=14, face="bold"),
        axis.title.y = element_text(color="#404040", size=14, face="bold"),
        axis.text.y = element_text(color="#404040",  size=12, face = "italic")) +
		ylim(0, 50)
plot(hMDM_TB_plot)

#arrange them in one plot 
ggarranged <- ggarrange(bAM_MB_plot, bAM_TB_plot, hAM_TB_plot, hMDM_TB_plot, common.legend = TRUE)
ggarranged
