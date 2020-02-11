####################################
##  Function to plot correlations ##
##  between pairs of DEG analyses ##
##  Fidel Lozano-Elena            ##
##  CRAG - 11/02/2020             ##
####################################

# Description of the script ----
# Function to automatically plot correlations after differential expression analysis.
# The function expect data.frames (1 and 2) and the name columns of numeric values to correlate (FC_column_1 and FC_column_2)
# Additionally if a p-value column is provided, it is possible to color significative genes in the union and intersection (and calculate correlations)
# If labels = TRUE, genes with extremes values and corelated will be tagged in the plot (It is necessary to provide column name for gene names)
# If return_intersection = TRUE, genes significant in both, matrix1 and matrix2 will be returned in form of a vector


# Function: ----
DEG_correlation<-function(
    
    ## Required arguments
    data_1,
    data_2,
    FC_column_1,
    FC_column_2,
    
    ## Default arguments
    corr.method = "pearson",
    colors = c("black","orange","darkred"),
    main = "",
    xlab = "Fold Change 1",
    ylab = "Fold Change 2",
    
    ## Optional arguments
    union_sigs = FALSE,
    intersection_sigs = FALSE,
    p.val_column_1,
    p.val_column_2,
    p.val_threshold = 0.05,
    labels = FALSE,
    gene_labels,
    return_intersection = FALSE
    
) {
    # Check that FC vectors has the same length and that are named and that their intersection is large
    if (is.null(rownames(data_1))|is.null(rownames(data_2))) {
        print("Data frames provided should have genes as row names")
        stop()
    }
    
    if (length(intersect(rownames(data_1),rownames(data_2))) < length(rownames(data_1))*0.5) {
        print("WARNING: data_1 and data_2 share less than 50% of the data_1 features... Take results with caution")
    }
    
    if (length(rownames(data_1)) != length(rownames(data_2))) {
        print("WARNING: FC_1 and FC_2 arguments have not the same length. Only genes present in both vectors will be considered for the plot")
        data_1<-data_1[intersect(rownames(data_1),rownames(data_2)),]
        data_2<-data_2[intersect(rownames(data_1),rownames(data_2)),]
    }
    
    # Reorder matrixes
    intersection<-intersect(rownames(data_1),rownames(data_2))
    data_1<-data_1[intersection,]
    data_2<-data_2[intersection,]
    
    # Plot
    plot(x=data_1[,FC_column_1],
         y=data_2[,FC_column_2],
         pch=16,cex=0.7,col=as.vector(colors)[1],
         main=main,
         xlab=xlab,
         ylab=ylab)
    
    # Guide
    abline(a = 0,b = 1,lty=2,col="tomato")
    
    # Correlation
    cors<-c(round(cor(x = data_1[,FC_column_1],
                      y = data_2[,FC_column_2],
                      use = "pairwise.complete.obs",
                      method = corr.method),2))
    
    # Draw poins of DEG genes (union AND/OR intersections)
    if (union_sigs == TRUE & intersection_sigs == FALSE ) {
        
        myset<-union(which(data_1[,p.val_column_1] < p.val_threshold),
                     which(data_2[,p.val_column_2] < p.val_threshold))
        points(x = data_1[myset,FC_column_1],
               y = data_2[myset,FC_column_2],
               pch=16,cex=0.5,col=colors[2])
        cors<-c(cors,round(cor(x = data_1[myset,FC_column_1],
                               y = data_2[myset,FC_column_2],
                               use = "pairwise.complete.obs",
                               method = corr.method),2))
        ## Legend
        legend("bottomright",legend = paste(c("Overall cor:","DEG union:"),cors,sep = " "),
               bty="n",text.col=colors[c(1,2)])
        
    } else if (union_sigs == FALSE & intersection_sigs == TRUE) {
        
        myset<-intersect(x = which(data_1[,p.val_column_1] < p.val_threshold),
                         y = which(data_2[,p.val_column_2] < p.val_threshold))
        points(x = data_1[myset,FC_column_1],
               y = data_2[myset,FC_column_2],
               pch=16,cex=0.5,col=colors[3])
        cors<-c(cors,round(cor(x = data_1[myset,FC_column_1],
                               y = data_2[myset,FC_column_2],
                               use = "pairwise.complete.obs",
                               method = corr.method),2))
        ## Legend
        legend("bottomright",legend = paste(c("Overall cor:","DEG intersect:"),cors,sep = " "),
               bty="n",text.col=colors[c(1,3)])
        
    } else if (union_sigs == TRUE & intersection_sigs == TRUE) {
        
        myset<-union(which(data_1[,p.val_column_1] < p.val_threshold),
                     which(data_2[,p.val_column_2] < p.val_threshold))
        points(x = data_1[myset,FC_column_1],
               y = data_2[myset,FC_column_2],
               pch=16,cex=0.5,col=colors[2])
        cors<-c(cors,round(cor(x = data_1[myset,FC_column_1],
                               y = data_2[myset,FC_column_2],
                               use = "pairwise.complete.obs",
                               method = corr.method),2))
        myset<-intersect(x = which(data_1[,p.val_column_1] < p.val_threshold),
                         y = which(data_2[,p.val_column_2] < p.val_threshold))
        points(x = data_1[myset,FC_column_1],
               y = data_2[myset,FC_column_2],
               pch=16,cex=0.5,col=colors[3])
        cors<-c(cors,round(cor(x = data_1[myset,FC_column_1],
                               y = data_2[myset,FC_column_2],
                               use = "pairwise.complete.obs",
                               method = corr.method),2))
        ## Legend
        legend("bottomright",legend = paste(c("Overall cor:","DEG union:","DEG intersect:"),cors,sep = " "),
               bty="n",text.col=colors)
        
    } else {
        
        ## Legend
        legend("bottomright",legend = paste("Overall cor:",cors,sep = " "),
               bty="n",text.col=colors[1])
    }
    
    # Text with extreme values (and correlated)
    if (labels) {
        # Highly correlated genes
        very_high<-which(data_1[,FC_column_1] > quantile(data_1[,FC_column_1],0.9995,na.rm = TRUE) & data_2[,FC_column_2] > quantile(data_2[,FC_column_2],0.9995,na.rm = TRUE))
        text(x = data_1[very_high,FC_column_1],
             y = data_2[very_high,FC_column_2],
             labels = data_1[very_high,gene_labels],
             cex = 0.7,pos=1)
        very_low<-which(data_1[,FC_column_1] < quantile(data_1[,FC_column_1],0.0005,na.rm = TRUE) & data_2[,FC_column_2] < quantile(data_2[,FC_column_2],0.0005,na.rm = TRUE))
        text(x = data_1[very_low,FC_column_1],
             y = data_2[very_low,FC_column_2],
             labels = data_1[very_low,gene_labels],
             cex = 0.7,pos=1)
    }
    
    # Return intersection genes??
    if (return_intersection) {
        myset<-intersect(x = which(data_1[,p.val_column_1] < p.val_threshold),
                         y = which(data_2[,p.val_column_2] < p.val_threshold))
        return(rownames(data_1)[myset])
    }
}

# END OF THE SCRIPT