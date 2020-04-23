####################################
##  GO Enrichment Analysis using  ##
##  Plant GO Slim annotations     ##
##  (Using ClusterProfiler)       ##
##  Fidel Lozano-Elena            ##
##  Barcelona - 20/04/2020        ##
####################################


# Description of the script ----
# Gene Ontology Enrichment Analysis using Plant GO Slim annotations:
# Slim annotations summarize general categories of functional annotations.
# This script uses Arabidopsis annotation database and the ClusterProfiler generic enrichment function for the enrichment analysis of (only) GO Slim terms
# Representation of the enrichment results is given in form of a barplot

# Body of the script ----

# Load dependencies
library(ontologyIndex)
library(clusterProfiler)
library(org.At.tair.db)

### 1 - LOAD AND PARSE ANNOTATION DATA
# Get Plant Slim annotations(Original obo file from current.geneontology.org/ontology/subsets/goslim_plant.obo)
go_slim<-get_ontology(file = "https://raw.githubusercontent.com/fle1/canolab_scripts/master/useful_data/goslim_plant.obo", propagate_relationships = "is_a", extract_tags = "minimal")

# Load all Arabidopsis GO annotations and mapping to genes
GO2TAIR <- as.list(org.At.tairGO2ALLTAIRS)

# Keep only GO categories of Plant Slim annotation with associated TAIR names
slim_annotations<-GO2TAIR[which(names(GO2TAIR)%in%go_slim$id)]

# Keep only term ID as names of the list with annotations (Save information on evidences mapping GO-Genes in a list for future)
evidences<-slim_annotations
for (i in 1:length(slim_annotations)) {
  slim_annotations[[i]]<-unname(slim_annotations[[i]])
}

# OPTIONAL: Delete the three source ontologies (they are not informative): MF:GO:0003674; CC:GO:0005575; BP:GO:0008150
slim_annotations<-slim_annotations[-which(names(slim_annotations)%in%c("GO:0003674","GO:0005575","GO:0008150"))]

# Generate data frames with teh mapping GOSlimID-Genes and GOSlimID-Description
GOSlim2Genes<-data.frame(TERM=names(unlist2(slim_annotations)),
                         GENE=unlist2(slim_annotations))

GOSlim2Names<-data.frame(TERM=names(slim_annotations),
                         NAME=go_slim$name[match(names(slim_annotations),go_slim$id)])
                         
                         
### 2 - LOAD USERS DATA AND DEFINE GENE OF INTEREST
# Define genes of interest as character vector                        
DE<-read.delim("Differential_exression_analysis_results.txt")              
mygenes<-rownames(subset(DE, FDR<0.05 & logFC < log(0.5)))

# Define universe (Optional, e.g. features identified in RNAseq; if not chosen the enrichment will be performed with the universe of genes annotated in GO Slims)
myuniverse<-rownames(DE)

### 3 - ENRICHMENT
# Create ClusterProfiler object with the generic enrichment function (It uses hypergeometric test)
ego<-enricher(gene = mygenes,
              TERM2GENE = GOSlim2Genes,
              TERM2NAME = GOSlim2Names,
              universe = myuniverse,
              minGSSize = 1,
              maxGSSize = 5000,
              pvalueCutoff = 1,
              qvalueCutoff = 1,
              pAdjustMethod = "BH")

# Convert results to a data frame
myresults<-as.data.frame(ego@result)

### 4 - REPRESENT ENRICHMENT RESULTS
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
          main = "GOSlimEA - My genes")
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

### 5 - DEPLOY GENES OF A PARTICULAR CATEGORY
# Deploy genes in a category
strsplit(myresults["GO:0003677","geneID"],split = "/")[[1]]
# Alternatively, use grep
strsplit(myresults[grep("DNA-binding",myresults$Description),"geneID"],split = "/")[[1]]

# Check genes in user's matrix
myset<-DE[strsplit(myresults[grep("DNA-binding",myresults$Description),"geneID"],split = "/")[[1]],]
View(myset)
