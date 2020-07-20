##########################DE correlation plot#####################################
#For doing a Log2fc correlation plot of the common genes. 
library(RColorBrewer)
library(corrplot)

Correlation <- read_xlsx("Correlation.xlsx")
Correlation2 <- data.frame(Correlation[,-1], row.names=Correlation$`Gene ID`)
Correlation.result <- cor(Correlation2)
corrplot(Correlation.result, method="pie", type="upper", col=brewer.pal(n=8, name="RdYlBu"))
