###########################
##  Function to plot PCA ##
##  Fidel Lozano Elena   ##
##  Barcelona 2017       ##
###########################

# Calculate PCA in R:
## PCA <- princomp(scale(mymatrix,scale=F), cor = TRUE)

# Define ploting function

myPCA.plot <- function (PCA, main, col = rep(1,7)) {
  P <- PCA$loadings
  TOTvar <- cumsum(PCA$sdev^2/sum(PCA$sdev^2))
  varPC1 <- round(TOTvar[1] * 100)
  varPC2 <- round((TOTvar[2] - TOTvar[1])*100) 
  plot (P[,1], P[,2], 
        col = "grey", pch=4,
        main = main, 
        xlim = c(range(P[,1])[1]* 1.05,range(P[,1])[2]*0.95),
        ylim = range(P[,2]) * 1.25,
        xlab = paste("PC1: ", varPC1, "% expl.var.", sep = ""),
        ylab = paste("PC2: ", varPC2, "% expl.var.", sep = "")
  ) ## plot sample values for the two first PCs
  text (P[,1], P[,2], rownames(P), offset = 0,
        col = col, cex = 0.9
  )
}
