####################################
##  GO Enrichment Analysis using  ##
##  Plant GO Slim annotations     ##
##  (Using topGO)                 ##
##  Fidel Lozano-Elena            ##
##  Barcelona - 20/04/2020        ##
####################################


# Description of the script ----
# Gene Ontology Enrichment Analysis using Plant GO Slim annotations. Slim annotations summarize general categories of functional annotations.
# This script uses Arabidopsis annotation and topGO packages for the enrichment analysis
# Represntation of the enrichment results is given in form of barplots

# Body of the script ----

# Load dependencies
library(ontologyIndex)
library(topGO)
library(GO.db)
library(org.At.tair.db)

# Download Plant Slim annotations (current.geneontology.org/ontology/subsets/goslim_plant.obo)
download.file(url = "current.geneontology.org/ontology/subsets/goslim_plant.obo",
              destfile = "./goslim_plant.obo")
# Read Plant Slim Annotations
go_slim<-get_ontology("./goslim_plant.obo", propagate_relationships = "is_a", extract_tags = "minimal")

# Load all Arabidopsis GO annotations and mapping to genes
GO2TAIR <- as.list(org.At.tairGO2TAIR)

# Keep only GO categories of Plant Slim annotation with associated TAIR names
slim_annotations<-GO2TAIR[which(names(GO2TAIR)%in%go_slim$id)]

# Keep the information of the gene-GO annotation in other vector (And delete from the names of ID vectors)
evidences<-list()
for (i in 1:length(slim_annotations)) {
  evidences[[i]]<-names(slim_annotations[[i]])
  slim_annotations[[i]]<-unname(slim_annotations[[i]])
}
names(evidences)<-names(slim_annotations)

# Define genes of interest as a named integer vector (Names: GeneID, Integers: 1 for DEG, 0 for NonDEG)
mygenes<-factor(as.integer(rownames(DE) %in% rownames(subset(DE,FDR<0.05&logFC>log(1.5)))))
names(mygenes)<-rownames(DE)


setwd("/Users/fidel/Dropbox/IDRICA03_seq/Pairwise_brl3_CTRLvs.WT_CTRL/")
brl3CTRL<- read.delim("./All_genes.txt",sep="\t",header = TRUE)

mygenes<-factor(as.integer(rownames(brl3CTRL) %in% rownames(subset(brl3CTRL,FDR<0.05&logFC<log(0.5)))))
names(mygenes)<-rownames(brl3CTRL)

# Create a topGO class object. 
# Use Custom Annotation function: gene 2 GO
topGO<-new(Class = "topGOdata",
           description="My_genes",
           ontology="BP",
           allGenes=mygenes,
           nodeSize=1o,
           annotationFun= annFUN.GO2genes, GO2genes = slim_annotations)


# GO enrichment analysis. Classic algorithm with a Fisher test (See topGO manual for details and other algorithms )
myresults<-runTest(topGO,algorithm = "classic",statistic = "fisher",topNodes=0)

# Enriched GO Slim categories
which(score(myresults)<0.01)
length(score(myresults))
which(names(myresults)%in%names(slim_annotations)
# Distribution of p-values
hist(score(myresults), 50, xlab = "p-values")

# Save results in form of a table (human redeable)
GO<-GenTable(object = topGO, myresults)
View(GO)

printGenes(topGO,whichTerms = GO$GO.ID[2])
# Represent enrichment values in form of barplots
## Change margins, so GO description is readable
par(mar=c(5,24,2,2)+0.1)
{
  barplot(rev(-log(as.numeric(GO$result1))),
          horiz = TRUE,
          names.arg = rev(GO$Term),
          col=heat.colors(n=max(GO$Significant))[rev(GO$Significant[order(GO$Significant,decreasing = TRUE)])],
          las =1, 
          xlab = "-log(p-value)", 
          main="My genes"
  )
  # Add guides
  abline(v = seq(0,-log(min(as.numeric(GO$result1))),1),
         lty = 3,
         col = "gray65")
  # Add line for default significance tHreshold
  abline(v = -log(0.05),
         lty = 2, 
         col = "darkgreen")
  # Add legend with color scale limits with the number of annotated genes per GO slim
  legend(x = "bottomright",
         legend = c(paste(max(GO$Significant),"genes",sep = " "),paste(min(GO$Significant),"genes",sep = " ")),
         title = "Nº annotated genes:",
         pch = 21,
         bty = "n",
         pt.bg = heat.colors(n=max(GO$Significant))[c(max(GO$Significant),1)],
         col = "black"
  )
}

# Reset default
dev.off()

# Deploy particular categories (Generate a list)
genes_in_GO<-list()
evidences_my_genes<-list()

for (i in 1:length(GO$GO.ID)) {
  myDEG<-names(mygenes)[which(mygenes==1)]
  myGOterm<-GO$GO.ID[i]
  slim_annotations[[which(names(slim_annotations)==as.character(GO$GO.ID[i]))]
  evidences_my_genes[i]<-evidences[[i]][which(slim_annotations[[GO$GO.ID[i]]]%in%myDEG)]
  genes_in_GO[i]<-myDEG[which(myDEG%in%slim_annotations[[GO$GO.ID[i]]])]
}

"GO:0050794"%in%slim_annotations

genes_in_GO<-names(mygenes)[which(mygenes==1)]
genes_in_GO<-genes_in_GO[which(genes_in_GO%in%slim_annotations[[GO$GO.ID[1]]])]
View(brl3CTRL[genes_in_GO,])

#### END OF THE SCRIPT
