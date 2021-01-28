###########################################
##  Plot to represent fold changes in    ##
##  metabolomic data respect a control   ##
##  Fidel Lozano-Elena                   ##
##  CRAG - 16/01/2020                    ##
###########################################

# Description of the script ----

# Script to automatically boxplot metabolomic data of two different "mutants" relativized to a common control
# Metabolites should be as row names of the matrixes and samples as colum names
# Matrixes should be ordered in the same way
# It performs t.tests (default behaviour) and only represent significant metabolites for either of the mutants (default p.val threshold = 0.05)
# metabolite_FC_dotplot performs the same basic function but representing the median with dots and the sd with bars

# Functions ----
metabolite_FC_boxplots<-function(
    
    # Required arguments:
    genotype1_matrix,
    genotype2_matrix,
    control_matrix,
    
    # Default arguments:
    normalization = "median",
    dotplot = TRUE,
    logarithmic = TRUE,
    horizontal = TRUE,
    color = c("#c71515","#08306b","darkgreen"),
    title = "",
    ylab = "Relative change",
    ylim = NULL,
    
    # Optional arguments:
    order = TRUE,
    t.test = TRUE,
    only_significant = FALSE,
    p_val_treshold = 0.05,
    separation_factor = 1,
    guides = TRUE
) {
    # Check that the matrixes have the same number of rows:
    if (nrow(genotype1_matrix)!=nrow(genotype2_matrix)) {
        print("ERROR: Matrix 1 has not the same number of rows than matrix 2")
        stop()
    } else if (nrow(genotype1_matrix)!=nrow(control_matrix)) {
        print("ERROR: Matrix 1 has not the same number of rows than the control matrix")
        stop()
    } else if (nrow(genotype2_matrix)!=nrow(control_matrix)) {
        print("ERROR: Matrix 2 has not the same number of rows than the control matrix")
        stop()
    }
    # Perform t.test
    if (t.test) {
        stat_1<-c()
        stat_2<-c()
        for (i in c(1:nrow(genotype1_matrix))) {
            stat_1<-c(stat_1,t.test(genotype1_matrix[i,],control_matrix[i,])$p.value)
            stat_2<-c(stat_2,t.test(genotype2_matrix[i,],control_matrix[i,])$p.value)
        }
        names(stat_1)<-rownames(genotype1_matrix)
        names(stat_2)<-rownames(genotype2_matrix)
        print(data.frame(p_vals_matrix1=stat_1,p_vals_matrix2=stat_2,row.names = names(stat_1)))
    }
    
    # Proceed only with significant metabolites:
    if (only_significant == TRUE) {
        union_met<-union(names(which(stat_1<p_val_treshold)),names(which(stat_2<p_val_treshold)))
        genotype1_matrix<-genotype1_matrix[union_met,]
        genotype2_matrix<-genotype2_matrix[union_met,]
        control_matrix<-control_matrix[union_met,]
    }
    
    # Normalization with the mean/median
    if (normalization == "mean") {
        val_1<-genotype1_matrix/apply(control_matrix,1,mean,na.rm=T)
        val_2<-genotype2_matrix/apply(control_matrix,1,mean,na.rm=T)
    } else {
        val_1<-genotype1_matrix/apply(control_matrix,1,median,na.rm=T)
        val_2<-genotype2_matrix/apply(control_matrix,1,median,na.rm=T)
    }
    
    # Logarithmic scale
    if (logarithmic == TRUE) {
        val_1<-log(val_1)
        val_2<-log(val_2)
    }
    
    # Reorder
    if (order == TRUE) {
        val_1<-val_1[order(apply(val_1,1,median),apply(val_2,1,median),decreasing = TRUE),]
        val_2<-val_2[rownames(val_1),]
    }
    
    # Determine plot positions
    n_metabolites<-nrow(genotype1_matrix)
    pos<-c(seq(1,n_metabolites*2))
    
    ## Modify the position according the separation factor
    # Avoid errors if there is 2 or less metabolites
    if (n_metabolites < 2) {
        pos<-c(seq(1,n_metabolites*2))
    } else {
        for (i in seq(1,n_metabolites*2)) {
            if (i%in%c(seq(3,n_metabolites*2,by = 2))) {
                pos[i:length(pos)]<-pos[i:length(pos)]+separation_factor
            }
        }
    }
    
    # Plot
    if (horizontal == TRUE) {
        ## Horizontal
        
        ### From top to bottom
        boxplot(t(val_2[c(nrow(val_2):1),]),at = pos[seq(1,length(pos)-1,by = 2)],cex=0.7,las=2,col = color[2], outline = FALSE, horizontal = horizontal, ylim = if (is.null(ylim)){ylim = c(min(min(val_1,na.rm = T),min(val_2,na.rm = T)),max(max(val_1,na.rm = T),max(val_2,na.rm = T)))} else{ylim = ylim},xlim=c(0,max(pos)+1),yaxt = "n",xlab = if (logarithmic==TRUE & ylab == "Relative change"){ylab="log(Relative change)"}else{ylab=ylab},main=title)
        boxplot(t(val_1[c(nrow(val_1):1),]),at = pos[seq(2,length(pos),by = 2)],cex=0.7,las=2,col = color[1],outline = FALSE, horizontal = horizontal, add = TRUE,yaxt = "n")
        ### Draw y-axis
        axis(side = 2, at = (pos[seq(1,length(pos)-1,by = 2)]+pos[seq(2,length(pos),by = 2)])/2,
             labels = rownames(val_1)[c(nrow(val_1):1)],las = 2)
        ### Guides
        if (guides == TRUE) {
            abline(h=(pos[seq(1,length(pos)-1,by = 2)]+pos[seq(2,length(pos),by = 2)])/2,lty=3,col="gray70")
        }
        ### Points
        if (dotplot == TRUE) {
            points(x = as.vector(t(val_2[c(nrow(val_2):1),])),
                   y = jitter(rep(pos[seq(1,length(pos)-1,by = 2)],each=ncol(val_1)),0.5),
                   pch = 21,cex = 0.5,bg="gray50")
            points(x = as.vector(t(val_1[c(nrow(val_1):1),])),
                   y = jitter(rep(pos[seq(2,length(pos),by = 2)],each=ncol(val_1)),0.5),
                   pch = 21,cex = 0.5,bg="gray50")
        }
        ### Asterisks? (to be implemented)
        
    } else {
        
        ## Vertical boxplot
        boxplot(t(val_1),at = pos[seq(1,length(pos)-1,by = 2)],cex=0.7,las=2,col = color[1],outline = FALSE, horizontal = horizontal,ylim=c(min(min(val_1,na.rm = T),min(val_2,na.rm = T)),max(max(val_1,na.rm = T),max(val_2,na.rm = T))),xlim=c(0,max(pos)+1),axt = "n",ylab = if (logarithmic==TRUE & ylab == "Relative change"){ylab="log(Relative change)"}else{ylab=ylab},main=title)
        boxplot(t(val_2),at = pos[seq(2,length(pos),by = 2)],cex=0.7,las=2,col = color[2], outline = FALSE, horizontal = horizontal, add = TRUE,xaxt = "n")
        
        ### Draw y-axis
        axis(side = 1, at = (pos[seq(1,length(pos)-1,by = 2)]+pos[seq(2,length(pos),by = 2)])/2,
             labels = rownames(val_1)[c(nrow(val_1):1)],las = 2)
        ### Guides
        if (guides == TRUE) {
            abline(v=(pos[seq(1,length(pos)-1,by = 2)]+pos[seq(2,length(pos),by = 2)])/2,lty=3,col="gray70")
        }
        ### Points
        if (dotplot == TRUE) {
            points(y = as.vector(t(val_1)),
                   x = jitter(rep(pos[seq(1,length(pos)-1,by = 2)],each=ncol(val_1)),0.5),
                   pch = 21,cex = 0.5,bg="gray50")
            points(y = as.vector(t(val_2)),
                   x = jitter(rep(pos[seq(2,length(pos),by = 2)],each=ncol(val_1)),0.5),
                   pch = 21,cex = 0.5,bg="gray50")
        }
    }
    
    # WT base line
    if (horizontal == TRUE) {
        abline(v=(if (logarithmic) v=0 else v=1),col=color[3],lty=2,lwd=1.5)
    } else {
        abline(h=(if (logarithmic) h=0 else h=1),col=color[3],lty=2,lwd=1.5)
    }
    
}
