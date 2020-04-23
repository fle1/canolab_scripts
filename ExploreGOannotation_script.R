##==============================##
##  Explore GO annotations and  ##
##  annotated Arabidopsis genes ##
##  Fidel Lozano-Elena          ##
##  Barcelona - 22/04/2020      ##
##==============================##

# Description of the script ----
# Exploration of Gene Ontology annotation trees using GO.db tools
# Selected terms can be matched with Arabdidopsis GO annotations 

# Body of the script ---- Explore GO categories

# Load dependencies
library(org.At.tair.db)
library(GO.db)

### 1 - LOAD COMPLETE ONTOLOGY (BP) AND ARABIDOPSIS ANNOTATION (See GO.db documentation)

# GO Terms (terms object contains all information, GOterms only the names)
terms<-as.list(GOTERM)
GOterms<-lapply(terms,Term)

# Parent and children IDs (It is also possible to use ancestors and offspring)
children<- as.list(GOBPCHILDREN)
parents<-as.list(GOBPPARENTS)

# Load all Arabidopsis GO annotations and mapping to genes
GO2TAIR <- as.list(org.At.tairGO2ALLTAIRS)

### 3 - EXPLORE GO ANNOTATION AND SELECT CATEGORIES TO DEPLOY

# Find GO term by ID (e.g. Response to stress)
GOterms["GO:0006950"]
terms["GO:0006950"]

# Find GO term by pattern
GOterms[grep("response to stress",GOterms)]
GOterms[grep("photosynthesis",GOterms)]

# Find parents and childrens
parents["GO:0006950"]
GOterms[unlist(parents["GO:0006950"])]

children["GO:0006950"]
GOterms[unlist(children["GO:0006950"])]

# Create lists of genes annotated in a single term
DE<-read.delim("DEG.txt")
colnames(DE)
myset<-subset(DE,FDR<0.05&abs(logFC)>log(2))
myset<-myset[which(rownames(myset)%in%GO2TAIR[["GO:0006950"]]),c("ID","logFC","GeneName","Description")]

# Write list
write.table(x = myset, file = "Response_to_stress_DEG.txt",sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)

### 4 - SOURCE HEATMAP FUNCTION AND REPRESENT GENES OF INTEREST
source("https://raw.githubusercontent.com/fle1/canolab_scripts/master/DeployGO_heatmap_function.R")

png("Stress_genes.png",width = 18,height = 24,units = "cm",res=400)
go_heatmap(genes = rownames(subset(DE,FDR<0.05&abs(logFC)>log(2))),
           values = subset(DE,FDR<0.05&abs(logFC)>log(2))$logFC,
           go_parent = "GO:0006950",
           go_child = c("GO:0006970","GO:0009408","GO:0006979","GO:0009409","GO:0009414"),
           gene_names = paste(rownames(subset(DE,FDR<0.05&abs(logFC)>log(2))),
                              subset(DE,FDR<0.05&abs(logFC)>log(2))$GeneName,sep = " "),
           go_child_names = unlist(GOterms[c("GO:0006970","GO:0009408","GO:0006979","GO:0009409","GO:0009414")]),
           heatmap = TRUE,
           main = "Response to stress",
           cexRow = 0.4
)
dev.off()

### END OF THE SCRIPT
