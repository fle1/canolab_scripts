##==============================##
##  Function to deploy GO child ##
##  categories (from a parent)  ##
##  and plot them as heatmap    ##
##  Fidel Lozano-Elena          ##
##  Barcelona - 22/04/2020      ##
##==============================##

# Description of the script ----
# Selected child categories of a parent GO term can be deployed, so values for genes of interest 
# annotated in both terms (parent and children) are represented in form of heatmap
# If heatmap = FALSE, the output matrix is returned, which can be saved as variable
# gplots & RColoBrewer libraries are required


### 1 - DEFINE A FUNCTION TO DEPLOY GENES IN FORM OF HEATMAP
go_heatmap<-function(genes,values,go_parent,go_child,go_child_names=NULL,gene_names=NULL, heatmap=TRUE,main=NULL,cexRow=0.5) {
  
  ## Parse variables ####
  if (class(genes)!="character") {
    print("genes variables must be a character vector")
    stop()
  }
  if (class(values)!="numeric") {
    print("values must be a numeric vector")
    stop()
  }
  if (grep("GO:",go_parent)!=1) {
    print("go parent must be a single GO ID in the form of GO:XXXXXXX")
    stop()
  }
  if (class(go_child)!="character") {
    print("go_child argument must be a character vector of GO IDs (in form of GO:XXXXXXX) with length >= 2")
    stop()
  }
  if (length(go_child)<2) {
    print("go_child argument must be of length >= 2")
    stop()
  }
  if (length(genes)!=length(values)) {
    print("genes and values must have the same length")
    stop()
  }
  if (!is.null(go_child_names)) {
    if (length(go_child)!=length(go_child_names)) {
      print("The number of names for go_child categories do not match with the number of go_child categories provided")
      stop()
    }
  }
  if (!is.null(gene_names)) {
    if (length(gene_names)!=length(genes)) {
      print("The number of gene names provided must be exactly the same number of genes provided (and in the same order)")
      stop()
    }
  }
  # Check if Arabidopsis annoataion are loaded in the environment already ####
  if (!exists("GO2TAIR")) {
    
    print("Arabidopsis annotations are not loaded in the environment. They will be loaded now")
    
    # Check if org.At.tair.db is loaded
    if (!"org.At.tair.db" %in% (.packages())) {
      try(library("org.At.tair.db"))
    }
    
    GO2TAIR <<- as.list(org.At.tairGO2ALLTAIRS)
  }
  
  # Create the matrix of genes in selected categories with the given values (NA otherwise) ####
  myset<-data.frame(values=values,genes=genes,row.names = genes)
  myGO<-unname(GO2TAIR[[go_parent]])
  myset<-myset[which(rownames(myset)%in%myGO),]
  mymatrix<-data.frame(row.names = rownames(myset))
  for (i in 1:length(go_child)) {
    mymatrix<-cbind(mymatrix,rep(NA,nrow(mymatrix)))
    idx<-which(rownames(myset)%in%GO2TAIR[[go_child[i]]])
    mymatrix[idx,i]<-myset[idx,"values"]
  }
  # Add new names (if provided) ####
  if (!is.null(go_child_names)) {
    colnames(mymatrix)<-go_child_names
  } else {
    colnames(mymatrix)<-go_child
  }
  if (!is.null(gene_names)) {
    rownames(mymatrix)<-gene_names[match(rownames(mymatrix),genes)]
  }
  # Order the matrix and eliminate genes not annotated in any child go term ####
  mymatrix<-mymatrix[do.call(order,c(lapply(1:ncol(mymatrix), function(i) mymatrix[, i]),list(decreasing=TRUE))),]
  mymatrix<-mymatrix[-which(apply(mymatrix,1,sum,na.rm=TRUE)==0),]
  
  # Heatmap ####
  if (heatmap==TRUE) {
    
    # Check that gplot and RColorBrewer libraries are loaded
    if (!"gplots" %in% (.packages())) {
      try(library("gplots"))
    }
    if (!"RColorBrewer" %in% (.packages())) {
      try(library(RColorBrewer))
    }
    
    # My colors
    mycolors<-colorpanel(400,low="#08306b",high="#c71515",mid ="#ffffff")
    
    #Global layout of the heatmap
    lmat=rbind(c(3,3,3),c(2,1,0),c(4,4,0)) # Graph grid: 1 heatmap, 2:Col dendogram, 3:Row dendogram, 4:Color key
    
    # Heatmap proportions
    lwid=c(0.2,0.4,1) # Width (columns)
    lhei=c(0.15,1,0.15) # Heigth (Rows)
    
    # HEATMAP
    heatmap.2(x = data.matrix(mymatrix), # Data matrix
              margins = c(2,2),#margins for column and row names
              col=mycolors,
              na.color="white", #Color scale and NA in white
              scale="none", #No scaling of the color for  "row" or "column"
              dendrogram = "none", #No dendograms. "row" or "column"
              breaks <- seq(from=min(range(mymatrix)), to=max(range(mymatrix)), length.out=399), # Add hard limits to the range of the heatmap
              Rowv=FALSE,Colv=FALSE, #No reordering of the col/rows for dendogram (clustering)
              key=TRUE,keysize=0.2,key.xlab = "log(FC)", #Add Colorkey
              density.info="none", #No lines for density
              colsep=c(0:ncol(mymatrix)),rowsep = c(0,nrow(mymatrix)), # Separator lines
              sepwidth = c(0.025,0.1), sepcolor = "gray", # Width and color of separatior lines (row,col)
              trace="none",
              cexCol=0.75,cexRow=cexRow,srtCol = 45, #Size of the row/col labels and horizontal col labels. ºrotation collabel
              main = if (!is.null(main)) {main = main} else {main = go_parent}, #graph title
              lmat=lmat,lwid=lwid,lhei=lhei)# Heatmap layout
  } else {
    return(mymatrix)
  }
}
                                            
go_heatmap_dual<-function(genes,values,go_parent,go_child,go_child_names=NULL,gene_names=NULL, heatmap=TRUE,main=NULL,cexRow=0.5,asp=FALSE) {
  
  ## Parse variables ####
  if (class(genes)!="character") {
    print("genes variables must be a character vector")
    stop()
  }
  if (class(values)!="list") {
    print("values must be a list with two numeric vectors")
    stop()
  }
  if (grep("GO:",go_parent)!=1) {
    print("go parent must be a single GO ID in the form of GO:XXXXXXX")
    stop()
  }
  if (class(go_child)!="character") {
    print("go_child argument must be a character vector of GO IDs (in form of GO:XXXXXXX) with length >= 2")
    stop()
  }
  if (length(go_child)<2) {
    print("go_child argument must be of length >= 2")
    stop()
  }
  if (length(genes)!=length(values[[1]])) {
    print("genes and values must have the same length")
    stop()
  }
  if (length(values[[2]])!=length(values[[1]])) {
    print("values vectors must have the same length")
    stop()
  }
  if (!is.null(go_child_names)) {
    if (length(go_child)!=length(go_child_names)) {
      print("The number of names for go_child categories do not match with the number of go_child categories provided")
      stop()
    }
  }
  if (!is.null(gene_names)) {
    if (length(gene_names)!=length(genes)) {
      print("The number of gene names provided must be exactly the same number of genes provided (and in the same order)")
      stop()
    }
  }
  # Check if Arabidopsis annoataion are loaded in the environment already ####
  if (!exists("GO2TAIR")) {
    
    print("Arabidopsis annotations are not loaded in the environment. They will be loaded now")
    
    # Check if org.At.tair.db is loaded
    if (!"org.At.tair.db" %in% (.packages())) {
      try(library("org.At.tair.db"))
    }
    
    GO2TAIR <<- as.list(org.At.tairGO2ALLTAIRS)
  }
  
  # Create the matrix of genes in selected categories with the given values (NA otherwise) ####
  myset<-data.frame(values1=values[[1]],values2=values[[2]],genes=genes,row.names = genes)
  myGO<-unname(GO2TAIR[[go_parent]])
  myset<-myset[which(rownames(myset)%in%myGO),]
  mymatrix<-data.frame(row.names = rownames(myset))
  for (i in 1:length(go_child)) {
    mymatrix<-cbind(mymatrix,rep(NA,nrow(mymatrix)))
    idx<-which(rownames(myset)%in%GO2TAIR[[go_child[i]]])
    mymatrix[idx,2*i-1]<-myset[idx,"values1"]
    mymatrix[idx,2*i]<-myset[idx,"values2"]
  }
  # Add new names (if provided) ####
  if (!is.null(go_child_names)) {
    colnames(mymatrix)<-rep(go_child_names,each=2)
  } else {
    colnames(mymatrix)<-rep(go_child,each=2)
  }
  if (!is.null(gene_names)) {
    rownames(mymatrix)<-gene_names[match(rownames(mymatrix),genes)]
  }
  # Order the matrix and eliminate genes not annotated in any child go term ####
  mymatrix<-mymatrix[do.call(order,c(lapply(1:ncol(mymatrix), function(i) mymatrix[, i]),list(decreasing=TRUE))),]
  mymatrix<-mymatrix[-which(apply(mymatrix,1,sum,na.rm=TRUE)==0),]
  
  # Heatmap ####
  if (heatmap==TRUE) {
    
    # Check that gplot and RColorBrewer libraries are loaded
    if (!"gplots" %in% (.packages())) {
      try(library("gplots"))
    }
    if (!"RColorBrewer" %in% (.packages())) {
      try(library(RColorBrewer))
    }
    
    # My colors
    mycolors<-colorpanel(400,low="#08306b",high="#c71515",mid ="#ffffff")
    
    #Global layout of the heatmap
    lmat=rbind(c(3,3,3),c(2,1,0),c(4,4,0)) # Graph grid: 1 heatmap, 2:Col dendogram, 3:Row dendogram, 4:Color key
    
    # Heatmap proportions
    lwid=c(0.2,0.4,1) # Width (columns)
    lhei=c(0.15,1,0.15) # Heigth (Rows)
    
    # HEATMAP
    heatmap.2(x = data.matrix(mymatrix), # Data matrix
              margins = c(2,2),#margins for column and row names
              col=mycolors,
              na.color="white", #Color scale and NA in white
              scale="none", #No scaling of the color for  "row" or "column"
              dendrogram = "none", #No dendograms. "row" or "column"
              breaks <- seq(from=min(range(mymatrix)), to=max(range(mymatrix)), length.out=399), # Add hard limits to the range of the heatmap
              Rowv=FALSE,Colv=FALSE, #No reordering of the col/rows for dendogram (clustering)
              key=TRUE,keysize=0.2,key.xlab = "log(FC)", #Add Colorkey
              density.info="none", #No lines for density
              colsep=seq(0,ncol(mymatrix),2),rowsep = c(0,nrow(mymatrix)), # Separator lines
              sepwidth = c(0.025,0.1), sepcolor = "gray", # Width and color of separatior lines (row,col)
              trace="none",
              cexCol=0.75,cexRow=cexRow,srtCol = 45, #Size of the row/col labels and horizontal col labels. ºrotation collabel
              main = if (!is.null(main)) {main = main} else {main = go_parent}, #graph title
              lmat=lmat,lwid=lwid,lhei=lhei,asp=asp)# Heatmap layout
  } else {
    return(mymatrix)
  }
}
                                            
