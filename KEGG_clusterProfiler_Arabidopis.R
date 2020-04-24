##########################################
##  KEGG pathways Enrichment Analysis   ##
##  with ClusterProfiler (Arabidopsis)  ##
##  Fidel Lozano-Elena                  ##
##  CRAG - 06/02/2020                   ##
##########################################


# Description of the script ----
# KEGG pathway Enrichment Analysis in Arabidopsis. This script uses Arabidopsis annotation and ClusterProfiler packages
# (Check https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf)
# Three difererent enrichment results representations are provided. 

# Body of the script ----

# Load dependencies
library(clusterProfiler)
library(org.At.tair.db)

# Define gene of interest as character vector
mygenes<-rownames(subset(DE, FDR<0.05, logFC>log(1.5)))

# Define universe (Optional, e.g. features identified in RNAseq)
myuniverse<-rownames(DE)

# Create ClusterProfiler object (KEGG enrichment pathway)
ego<-enrichKEGG(gene     = mygenes,
                universe = keys(org.At.tair.db),
                organism = "ath")

# Take a look
head(summary(ego))

# Save results in form of data.frame
myresult<-ego@result

# Representations
## Dotplot:
dotplot(ego,showCategory=30,title="",font.size=8)

## Barplot
barplot(keggresults,showCategory = 15,title="")

## EMAP-plot:
emapplot(ego,showCategory = 30)

# Deploy genes in a category
strsplit(ego@result["ath00020","geneID"],split = "/")[[1]]

#### END OF THE SCRIPT
