###################################
##  Plot to explore fold-change  ##
##  thresholds (RNAseq)          ##
##  Fidel Lozano-Elena           ##
##  CRAG - 16/01/2020            ##
###################################


# Description of the script ----

# This script defines a function to automatically explore the number of differentially regulated genes (DEG),
# given a hard p-value (numeric) and a range of relative transcript change (vector, default relative changes = c(1,1.2,1.4,1.5,1.6,1.8,2)).
# Additional arguments might be given to change colors, title, draw number of genes at a specific FC treshold and to save plot (see below)

# Function: ----
plotFCthreshold<-function(
    ## Required arguments:
    logFC_vector,
    p_val_vector,
    
    ## Default arguments:
    p_val_threshold = 0.05,
    FC_thresholds = c(1,1.2,1.4,1.5,1.6,1.8,2),
    
    ## Ploting arguments:
    title = "",
    colors = c("gray55","#c71515","#08306b"),
    selected_threshold = FALSE,
    
    ## Save file:
    saveplot = FALSE,
    plotwidth_cm = 14,
    plotheigth_cm = 14,
    filename = "myfile"
) {
    ## Compute thresholds for total, up and down regulated genes
    logFC <- data.frame(logFC=logFC_vector,pval=p_val_vector)
    NoDEG <- c()
    NoGup <- c()
    NoGdown <- c()
    for (i in 1:length(FC_thresholds)) {
        NoDEG <- c(NoDEG,
                   length(which(logFC$pval<p_val_threshold& abs(logFC$logFC) > log(FC_thresholds[i])))
        )
        NoGup <- c(NoGup,
                   length(which(logFC$pval<p_val_threshold& logFC$logFC > log(FC_thresholds[i])))
        )
        NoGdown <- c(NoGdown,
                     length(which(logFC$pval<p_val_threshold& logFC$logFC < log(1/FC_thresholds[i])))
        )
    }
    
    ## Create labels
    if (1 %in% FC_thresholds) {
        labels = c("No threshold",paste(as.character((FC_thresholds[-1]-1)*100),"%",sep = ""))
    } else {
        labels = paste(as.character((FC_thresholds-1)*100),"%",sep = "")
    }
    
    ## Save plot?
    if (saveplot == TRUE) {
        png("NumberDEG.png",width = 14,height = 14,units = "cm",res = 300)
    }
    
    ## Plot:
    ### Total deregulated
    plot(
        x = seq(1,length(FC_thresholds),1),
        y = NoDEG,
        type = "b",
        bg = colors[1],
        pch = 21,
        xlab = "Relative transcript deregulation",
        ylab = "Number of genes", 
        main = title,
        xaxt = "n",
        yaxt = "n",
        ylim = c(0,max(NoDEG)),
        cex.lab=0.8,
        cex.axis=0.8
    )
    
    ### Upregulated
    lines(
        x = seq(1,length(FC_thresholds),1),
        y = NoGup,
        type = "b",
        bg = colors[2],
        pch = 21
    )
    
    ### Downregulated
    lines(
        x = seq(1,length(FC_thresholds),1),
        y = NoGdown,
        type = "b",
        bg = colors[3],
        pch = 21
    )
    
    ### Add custom axes
    axis(side = 2, at = round(seq(0,max(NoDEG),length.out = 6),-2), labels = round(seq(0,max(NoDEG),length.out = 6),-2), las=2, cex.lab=0.8, cex.axis=0.8) 
    axis(side = 1, at = seq(1,length(FC_thresholds),1), labels = labels, las=1, cex.lab=0.8, cex.axis=0.8) 
    
    ### Add guides
    abline(h=round(seq(0,max(NoDEG),length.out = 6),-2), col="gray65",lty=3,lwd=0.5)
    
    ### Add legend
    legend(
        "topright",
        pch = 21,
        pt.bg = colors,
        bty = "n",
        legend = c("Deregulated (total)", "Upregulated", "Downregulated")
    )
    ### Selected threshold
    if (selected_threshold != FALSE) {
        # Sep line
        abline(v= which(FC_thresholds==selected_threshold), lty=3,col="gray25")
        # Up
        text(x = which(FC_thresholds==selected_threshold),
             y = NoGup[which(FC_thresholds==selected_threshold)],
             pos = if (NoGup[which(FC_thresholds==selected_threshold)] > NoGdown[which(FC_thresholds==selected_threshold)]) {
                 3
             } else {
                 1
             },
             labels = as.character(NoGup[which(FC_thresholds==selected_threshold)]),
             col = colors[2],
             cex = 0.8
        )
        # Down
        text(x = which(FC_thresholds==selected_threshold),
             y = NoGdown[which(FC_thresholds==selected_threshold)],
             pos = if (NoGup[which(FC_thresholds==selected_threshold)] > NoGdown[which(FC_thresholds==selected_threshold)]) {
                 1
             } else {
                 3
             },
             labels = as.character(NoGdown[which(FC_thresholds==selected_threshold)]),
             col = colors[3],
             cex = 0.8
        )
    }
    
    ## Save plot?
    if (saveplot == TRUE) {
        dev.off()
    }
    
}
