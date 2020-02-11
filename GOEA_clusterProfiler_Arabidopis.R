###################################
##  GO Enrichment Analysis with  ##
##  ClusterProfiler (Arabidopsis)##
##  Fidel Lozano-Elena           ##
##  CRAG - 06/02/2020            ##
###################################


# Description of the script ----
# Gene Ontology Enrichment Analysis in Arabidopsis. This script uses Arabidopsis annotation and ClusterProfiler packages
# (Check https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf)
# Three difererent GOEA results representations are provided. Also GO categories can be "pruned" at a certain level to get more generar or specific categories

# Body of the script ----

# Load dependencies
library(clusterProfiler)
library(org.At.tair.db)

# Define gene of interest as character vector
mygenes<-rownames(subset(DE, FDR<0.05, logFC>log(1.5)))

# Define universe (Optional, e.g. features identified in RNAseq)
myuniverse<-rownames(DE)

# Create ClusterProfiler object
ego<-enrichGO(gene          = mygenes,
              universe      = keys(org.At.tair.db), # Include all Arabidopsis genes
              OrgDb         = org.At.tair.db,
              keyType       = "TAIR", # Tell the function to use TAIR names
              ont           = "BP", 
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              qvalueCutoff  = 0.05)

# Take a look
head(summary(ego))

# Representations
## Dotplot:
dotplot(ego,showCategory=30,title="My title",font.size=8)

## GO graph:
plotGOgraph(ego, firstSigNodes = 10,useInfo = "all")

## EMAP-plot:
emapplot(ego,showCategory = 30)

# GO level filtering
ego_filtered<-gofilter(ego, level = 4) # Same functions than above can be applied to ego_filtered

# Access the data
ego@result

# Deploy genes in a category
strsplit(ego@result["GO:0000690","geneID"],split = "/")[[1]]

#### END OF THE SCRIPT