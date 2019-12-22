##################################
##  Script for ploting heatmaps ##
##  Nature Communications Style ##
##  (Fàbregas et al., 2018)     ##
##  Fidel Lozano Elena, PhD     ##
##  2018, CRAG, Barcelona       ##
##################################


# 1 - Load requiered libraries

library(RColorBrewer)
library(gplots)

# 2 - Create color palette (as in Fàbregas et al., 2018)
mycolors<-colorpanel(400,low="#08306b",high="#c71515",mid ="#ffffff")

# TESTING PURPOSES - Simulate data
mymatrix<-data.matrix(cbind(sample(rnorm(20,1,1)),sample(rnorm(10,1,2)),sample(rnorm(20,2,2))))

# 3 - Global layout of the heatmap
lmat=rbind(c(3,3,3),c(2,1,0),c(4,4,0)) # Graph grid: 1 heatmap, 2:Col dendogram, 3:Row dendogram, 4:Color key

# 4 - Relative size of each element of the graph
lwid=c(0.2,0.4,1) # Width (columns)
lhei=c(0.25,1,0.25) # Heigth (Rows)

# Alternatively set the size in cm
lwid=c(lcm(1.5),lcm(2.5),lcm(8))
lwid=c(lcm(2.5),lcm(10),lcm(2.5))

# 5- Save the plot?
#This code is for saving high resolution graphs, it is necessary adjsut the size (figure with margins has to fit), and then dev.off() to finally write the file
png("myheatmap.png",width = 12, height = 18, units = "cm",res = 400)

# 6 - HEATMAP
heatmap.2(x = mymatrix, # Data matrix
          margins = c(2,2),#margins for column and row names
          col=mycolors, na.color="white", #Color scale and NA in white
          scale="none", #No scaling of the color for  "row" or "column"
          dendrogram = "none", #No dendograms. "row" or "column"
          breaks <- seq(from=min(range(mymatrix)), to=max(range(mymatrix)), length.out=399), # Add hard limits to the range of the heatmap
          Rowv=FALSE,Colv=FALSE, #No reordering of the col/rows for dendogram (clustering)
          key=TRUE,keysize=0.2,key.xlab = "log2(FC)", #Add Colorkey
          density.info="none", #No lines for density
          colsep=c(0:ncol(mymatrix)),rowsep=c(0,6,18,nrow(mymatrix)), # Separator lines
          sepwidth = c(0.05,0.1), sepcolor = "gray", # Width and color of separatior lines (row,col)
          trace="none",labCol = "",
          cexCol=0.9,cexRow=0.7,srtCol = 45, #Size of the row/col labels and horizontal col labels. ºrotation collabel
          main= "My title", #graph title
          lmat=lmat,lwid=lwid,lhei=lhei) # Heatmap layout

# 7 - Shut off the device (Necessary for saving the plot)
dev.off()
