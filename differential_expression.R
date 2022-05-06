
setwd("C:/Users/hncetin/OneDrive/Belgeler/Programming/MATLAB/MATLAB Drive/Thesis-Model/manuscript-scripts")
dir.create("differential_toptables")


## Collect gene expression data (requires connection)

library(Biobase)
library(GEOquery)
library(limma)

platform <- getGEO("GPL22543")	
platform <- cbind(platform@dataTable@table$ID, platform@dataTable@table$SystematicName)
colnames(platform) <- c("ID", "Gene")

geo = list()
expressions = list()
factors = list()


# Experiment: ethanol
# b2
geo$ethanolb2       <- getGEO("GSE78759", GSEMatrix=TRUE, getGPL=FALSE)
expressions$ethanolb2 <- exprs(geo$ethanolb2[[1]])
expressions$ethanolb2 <- expressions$ethanolb2[,c(1,2,3,4,5)]
colnames(expressions$ethanolb2) <- c("ETH-REF","ETH-REF","ETH-REF", "ETH-B2","ETH-B2")

factors$ethanolb2 <- as.factor(colnames(expressions$ethanolb2))
# b8
geo$ethanolb8       <- getGEO("GSE78759", GSEMatrix=TRUE, getGPL=FALSE)
expressions$ethanolb8 <- exprs(geo$ethanolb8[[1]])
expressions$ethanolb8 <- expressions$ethanolb8[,c(1,2,3,6,7,8)]
colnames(expressions$ethanolb8) <- c("ETH-REF2","ETH-REF2","ETH-REF2", "ETH-B8","ETH-B8","ETH-B8")
factors$ethanolb8 <- as.factor(colnames(expressions$ethanolb8))


# Experiment: caffeine
geo$caffeine      <- getGEO("GSE124452", GSEMatrix=TRUE, getGPL=FALSE) 
expressions$caffeine <- exprs(geo$caffeine[[1]])
colnames(expressions$caffeine) <- c("CAF-905","CAF-905","CAF-905", "CAF-REF","CAF-REF","CAF-REF")
factors$caffeine <- as.factor(colnames(expressions$caffeine))

# Experiment: coniferylaldehyde
geo$coniferylaldehyde          <- getGEO("GSE119240", GSEMatrix=TRUE, getGPL=FALSE) 
expressions$coniferylaldehyde <- exprs(geo$coniferylaldehyde[[1]])
colnames(expressions$coniferylaldehyde) <- c("CON-BH13","CON-BH13","CON-BH13","CON-REF","CON-REF","CON-REF")
factors$coniferylaldehyde <- as.factor(colnames(expressions$coniferylaldehyde))


# Experiment: iron
geo$iron          <- getGEO("GSE61317", GSEMatrix=TRUE, getGPL=FALSE) 
expressions$iron <- exprs(geo$iron[[1]])
expressions$iron <- expressions$iron[,c(1:6)]
colnames(expressions$iron) <- c("IRN-REF","IRN-REF","IRN-REF", "IRN-M8FE","IRN-M8FE","IRN-M8FE")
factors$iron <- as.factor(colnames(expressions$iron))



# Experiment: nickel
geo$nickel          <- getGEO("GSE50985", GSEMatrix=TRUE, getGPL=FALSE) 
expressions$nickel <- exprs(geo$nickel[[1]])
colnames(expressions$nickel) <- c("NIC-REF","NIC-REF","NIC-REF", "NIC-M9","NIC-M9","NIC-M9")
factors$nickel <- as.factor(colnames(expressions$nickel))


# Experiment: phenylethanol
geo$phenylethanol <- getGEO("GSE59353", GSEMatrix=TRUE, getGPL=FALSE)
expressions$phenylethanol <- exprs(geo$phenylethanol[[1]])
colnames(expressions$phenylethanol) <- c("PHE-REF","PHE-REF","PHE-REF", "PHE-C9","PHE-C9","PHE-C9")
factors$phenylethanol <- as.factor(colnames(expressions$phenylethanol))


# Experiment: silver
geo$silver <- getGEO("GSE143335", GSEMatrix=TRUE, getGPL=FALSE)
expressions$silver <- exprs(geo$silver[[1]]) 
colnames(expressions$silver) <- c("SIL-2E", "SIL-2E", "SIL-2E", "SIL-REF", "SIL-REF", "SIL-REF")
factors$silver <- as.factor(colnames(expressions$silver))
####################



######################################################################################
## Boxplot of all expressions


library(ggplot2)
library(tidyverse)
library(stringr)
library(reshape)

experiments = names(expressions)
df = list()
fac = list()
for (i in 1:length(experiments)) {
  df[[i]] <- data.frame(expressions[[i]])
} 
names(df) = experiments
summary(df)

# Cbind the df and melt for boxplot
allDF <- c(df$ethanolb2, df$ethanolb8, df$caffeine, df$coniferylaldehyde, df$iron,df$nickel,df$phenylethanol, df$silver)
meltDF <- melt(allDF)
names(meltDF) <- c("Gene Expressions", "Strain Name")

# Plot BW
ggplot(meltDF, aes(x=`Strain Name`, y=`Gene Expressions`, fill="black")) +
  # geom_jitter(color="#525252", size=0.01, alpha=0.05) +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(breaks = round(seq(-2, 2, by = 1), 1),  limits = c(-2,2)) + 
  theme(text = element_text(size=18),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "#ffffff"), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "bottom", legend.box="horizontal") 


######################################################################################
## Fit linear models and calculate differential expressions

colnames(expressions$ethanolb2) <- c("REF","REF","REF", "EV","EV")
colnames(expressions$ethanolb8) <- c("REF","REF","REF", "EV","EV","EV")
colnames(expressions$caffeine) <- c("EV","EV","EV", "REF","REF","REF")
colnames(expressions$coniferylaldehyde) <- c("EV","EV","EV", "REF","REF","REF")
colnames(expressions$iron) <- c("REF","REF","REF", "EV","EV","EV")
colnames(expressions$nickel) <- c("REF","REF","REF", "EV","EV","EV")
colnames(expressions$phenylethanol) <- c("REF","REF","REF", "EV","EV","EV")
colnames(expressions$silver) <- c("EV", "EV", "EV", "REF","REF","REF")

factors$ethanolb2 <- as.factor(colnames(expressions$ethanolb2))
factors$ethanolb8 <- as.factor(colnames(expressions$ethanolb8))
factors$caffeine <- as.factor(colnames(expressions$caffeine))
factors$coniferylaldehyde <- as.factor(colnames(expressions$coniferylaldehyde))
factors$iron <- as.factor(colnames(expressions$iron))
factors$nickel <- as.factor(colnames(expressions$nickel))
factors$phenylethanol <- as.factor(colnames(expressions$phenylethanol))
factors$silver <- as.factor(colnames(expressions$silver))

experiments = names(expressions)
for (i in 1:length(experiments)) { 
  
  exp = experiments[i]
  # Fit linear model for each gene given a series of arrays
  designmat <- model.matrix(~ factors[[i]] + 0, geo[[i]])
  colnames(designmat) <- levels(factors[[i]])
  fit <- lmFit(expressions[[i]], designmat)
  fit <- eBayes(fit) 
  
  # Make pair-wise comparison between the REF and Evolved 
  
  design_contrast <- makeContrasts(EV-REF, levels=designmat)
  fit2 <- contrasts.fit(fit, design_contrast)
  fit2 <- eBayes(fit2)
  
  # Collect Results
  results <- decideTests(fit2)
  toptable <- topTable(fit2, adjust="fdr", p.value=0.05, number = dim(expressions[[i]])[1])
  
  # Save Results
  write.table(toptable, file=paste0("differential_toptables/", exp, "_toptable.tsv"), row.names=T, col.names=T)
}

write.table(platform, file=paste0("differential_toptables/platform_for_genes.tsv"), row.names=F, col.names=T)



## PLOT

setwd("differential_toptables")
experiments = c('ethanolb2','ethanolb8', 'caffeine', 'coniferylaldehyde', 'iron', 'nickel', 'phenylethanol', 'silver')


ethanolb2 <- read.table(file = 'ethanolb2_toptable.tsv') 
ethanolb8 <- read.table(file = 'ethanolb8_toptable.tsv') 
caffeine <- read.table(file = 'caffeine_toptable.tsv') 
coniferylaldehyde <- read.table(file = 'coniferylaldehyde_toptable.tsv') 
iron <- read.table(file = 'iron_toptable.tsv') 
nickel <- read.table(file = 'nickel_toptable.tsv') 
phenylethanol <- read.table(file = 'phenylethanol_toptable.tsv') 
silver <- read.table(file = 'silver_toptable.tsv') 

exprs_neg = list(B2 = rownames(subset(ethanolb2, logFC < 0)),
                 B8 = rownames(subset(ethanolb8, logFC < 0)),
                 Caf9052 = rownames(subset(caffeine, logFC < 0)),
                 BH13 = rownames(subset(coniferylaldehyde, logFC < 0)),
                 M8FE = rownames(subset(iron, logFC < 0)),
                 M9 = rownames(subset(nickel, logFC < 0)),
                 C9 = rownames(subset(phenylethanol, logFC < 0)),
                 E2 = rownames(subset(silver, logFC < 0)))

exprs_pos = list(B2 = rownames(subset(ethanolb2, logFC > 0)),
                 B8 = rownames(subset(ethanolb8, logFC > 0)),
                 Caf9052 = rownames(subset(caffeine, logFC > 0)),
                 BH13 = rownames(subset(coniferylaldehyde, logFC > 0)),
                 M8FE = rownames(subset(iron, logFC > 0)),
                 M9 = rownames(subset(nickel, logFC > 0)),
                 C9 = rownames(subset(phenylethanol, logFC > 0)),
                 E2 = rownames(subset(silver, logFC > 0)))

library("UpSetR")
library("ggplot")
pltA <- upset(fromList(exprs_pos), nsets = 8, order.by = "degree", set_size.show=T,text.scale = 1.5, empty.intersections = NULL)
pltA

pltB <-upset(fromList(exprs_neg), nsets = 8, order.by = "degree", set_size.show=T, text.scale = 1.5, empty.intersections = NULL)
pltB 


