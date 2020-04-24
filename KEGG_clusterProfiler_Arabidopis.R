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
myresults<-ego@result

# Representations
## Dotplot:
dotplot(ego,showCategory=30,title="",font.size=8)

## Barplot
barplot(keggresults,showCategory = 15,title="")

## EMAP-plot:
emapplot(ego,showCategory = 30)

# Represent enrichment values in form of barplots
## Change margins, so GO description is readable
par(mar=c(5,18,2,2)+0.1)
# Only signficant categories
sigs<-myresults[which(myresults$p.adjust<0.15),]
{
  barplot(rev(-log(sigs$p.adjust)),
          horiz = TRUE,
          names.arg = rev(sigs$Description),
          cex.names = 1,
          col = heat.colors(n=max(sigs$Count))[rev(sigs$Count[order(sigs$Count,decreasing = TRUE)])],
          las = 1, 
          xlab = "-log(adj. p-value)", 
          main = "KEGG pathways - My genes")
  # Add guides
  abline(v = seq(0,-log(min(as.numeric(sigs$p.adjust))),2.5),
         lty = 3,
         col = "gray65")
  # Add line for default significance tHreshold
  abline(v = -log(0.05),
         lty = 2, 
         col = "darkgreen")
  # Add legend with color scale limits with the number of annotated genes per GO slim
  legend(x = "bottomright",
         legend = c(paste(max(sigs$Count),"genes",sep = " "),paste(min(sigs$Count),"genes",sep = " ")),
         title = "Nº annotated genes:",
         pch = 21,
         bty = "n",
         pt.bg = heat.colors(n=max(sigs$Count))[c(max(sigs$Count),1)],
         col = "black"
  )
}

# DEPLOY SELECTED CATEGORIES
strsplit(ego@result["ath00020","geneID"],split = "/")[[1]]

# Gene names, logFC and Description
DE[strsplit(ego@result["ath00020","geneID"],split = "/")[[1]],c("logFC","GeneName","Description")]

#### END OF THE SCRIPT
