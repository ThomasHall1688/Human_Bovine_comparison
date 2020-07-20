######################################### Line plot ########################################
library(ggplot2)
library(ggthemes)
setwd("C:/Users/Thomas Hall/Dropbox/Human_bovine_comparison/Analysis/GWAS integration/GWASin Summaries/Line plot data")

rm(summary)
bAM_MB <- read.csv("bAM_MB.csv")
summary <- bAM_MB
summary <- summary[order(summary$SNPS_total),] 
summary$Distance <- as.character(summary$Distance)
colnames(summary)[9] <- "probability_0.1"
#summary$probability_0.1 <- summary$probability_0.1 * 0.01
summary$probability_0.1 <- -log10(summary$probability_0.1)
#summary$probability_0.1 <- abs(summary$probability_0.1)
my_title <- expression(paste("bAM infected with ", italic("Mycobacterium bovis"))) 

Line <- ggplot(data = summary,  aes(x= reorder(Distance, +SNPS_total), y=probability_0.1, group = Group, color = Group, shape = Group)) + 
			geom_line(size = 2.5, alpha = 0.7, show.legend = TRUE) +
			geom_point(size = 5, show.legend = TRUE) +
			ylim(-2, 3.5) +
			theme(axis.text=element_text(size=15, colour = "#595959"), axis.title=element_text(size=16, face = "bold", , colour = "#595959" )) +
			geom_hline(yintercept=-0.69897000433, linetype="dashed", color = "blue", alpha = 0.6) +
#			geom_hline(yintercept=1.30103, color = "red", alpha = 0.6) +
		    labs(y = "Probability(-log10)", x = "Distance from gene", title = my_title) 
 
plot(Line) +
theme_hc() + scale_colour_colorblind()
#theme_calc() + scale_color_calc()
#theme_hc() + scale_color_calc()
#theme_hc()+ scale_colour_hc()
#theme_economist() + scale_colour_economist() 

 

rm(summary)
bAM_MB <- read.csv("bAM_TB.csv")
summary <- bAM_MB
summary <- summary[order(summary$SNPS_total),] 
summary$Distance <- as.character(summary$Distance)
colnames(summary)[9] <- "probability_0.1"
#summary$probability_0.1 <- summary$probability_0.1 * 0.01
summary$probability_0.1 <- -log10(summary$probability_0.1)
#summary$probability_0.1 <- abs(summary$probability_0.1)
my_title <- expression(paste("bAM infected with ", italic("Mycobacterium tuberculosis"))) 

Line <- ggplot(data = summary,  aes(x= reorder(Distance, +SNPS_total), y=probability_0.1, group = Group, color = Group, shape = Group)) + 
			geom_line(size = 2.5, alpha = 0.7, show.legend = TRUE) +
			geom_point(size = 5, show.legend = TRUE) +
			ylim(-2, 3.5) +
			theme(axis.text=element_text(size=15, colour = "#595959"), axis.title=element_text(size=16, face = "bold", , colour = "#595959" )) +
			geom_hline(yintercept=-0.69897000433, linetype="dashed", color = "blue", alpha = 0.6) +
#			geom_hline(yintercept=1.30103, color = "red", alpha = 0.6) +
		    labs(y = "Probability(-log10)", x = "Distance from gene", title = my_title) 
 
plot(Line) +
theme_hc() + scale_colour_colorblind()
#theme_calc() + scale_color_calc()
#theme_hc() + scale_color_calc()
#theme_hc()+ scale_colour_hc()
#theme_economist() + scale_colour_economist() 


rm(summary)
bAM_MB <- read.csv("hAM_TB.csv")
summary <- bAM_MB
summary <- summary[order(summary$SNPS_total),] 
summary$Distance <- as.character(summary$Distance)
colnames(summary)[9] <- "probability_0.1"
#summary$probability_0.1 <- summary$probability_0.1 * 0.01
summary$probability_0.1 <- -log10(summary$probability_0.1)
#summary$probability_0.1 <- abs(summary$probability_0.1)
my_title <- expression(paste("hAM infected with ", italic("Mycobacterium tuberculosis"))) 

Line <- ggplot(data = summary,  aes(x= reorder(Distance, +SNPS_total), y=probability_0.1, group = Group, color = Group, shape = Group)) + 
			geom_line(size = 2.5, alpha = 0.7, show.legend = TRUE) +
			geom_point(size = 5, show.legend = TRUE) +
			ylim(-2, 3.5) +
			theme(axis.text=element_text(size=15, colour = "#595959"), axis.title=element_text(size=16, face = "bold", , colour = "#595959" )) +
			geom_hline(yintercept=-0.69897000433, linetype="dashed", color = "blue", alpha = 0.6) +
#			geom_hline(yintercept=1.30103, color = "red", alpha = 0.6) +
		    labs(y = "Probability(-log10)", x = "Distance from gene", title = my_title) 
 
plot(Line) +
theme_hc() + scale_colour_colorblind()
#theme_calc() + scale_color_calc()
#theme_hc() + scale_color_calc()
#theme_hc()+ scale_colour_hc()
#theme_economist() + scale_colour_economist() 


rm(summary)
bAM_MB <- read.csv("hMDM_TB.csv")
summary <- bAM_MB
summary <- summary[order(summary$SNPS_total),] 
summary$Distance <- as.character(summary$Distance)
colnames(summary)[9] <- "probability_0.1"
#summary$probability_0.1 <- summary$probability_0.1 * 0.01
summary$probability_0.1 <- -log10(summary$probability_0.1)
#summary$probability_0.1 <- abs(summary$probability_0.1)
my_title <- expression(paste("hMDM infected with ", italic("Mycobacterium tuberculosis"))) 

Line <- ggplot(data = summary,  aes(x= reorder(Distance, +SNPS_total), y=probability_0.1, group = Group, color = Group, shape = Group)) + 
			geom_line(size = 2.5, alpha = 0.7, show.legend = TRUE) +
			geom_point(size = 5, show.legend = TRUE) +
			ylim(-2, 3.5) +
			theme(axis.text=element_text(size=15, colour = "#595959"), axis.title=element_text(size=16, face = "bold", , colour = "#595959" )) +
			geom_hline(yintercept=-0.69897000433, linetype="dashed", color = "blue", alpha = 0.6) +
#			geom_hline(yintercept=1.30103, color = "red", alpha = 0.6) +
		    labs(y = "Probability(-log10)", x = "Distance from gene", title = my_title) 
 
plot(Line) +
theme_hc() + scale_colour_colorblind()
#theme_calc() + scale_color_calc()
#theme_hc() + scale_color_calc()
#theme_hc()+ scale_colour_hc()
#theme_economist() + scale_colour_economist() 
