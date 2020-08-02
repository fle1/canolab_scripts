##############################
##  Handy functions to be   ##
##  used with boxplots      ##
##  Fidel Lozano, PhD       ##
##  Barcelona, 13/11/2019   ##
##############################

# 1- Define functions ####

# Automatic spacing per factor (or factor combinations) for box/violin plots
automatic_spacing<-function(factor1,factor2=NULL,factor3=NULL,sep_factor=1) {
  
  # Number of levels per factors
  n_levels1<-length(levels(factor1))
  n_levels2<-length(levels(factor2))
  n_levels3<-length(levels(factor3))
  
  # Parse optional arguments
  if (is.null(factor3)&is.null(factor2)) {
    total_levels<-n_levels1
    min_factor<-n_levels1
  } else if (is.null(factor3)) {
    total_levels<-n_levels1*n_levels2
    min_factor<-min(c(n_levels1,n_levels2))
  } else {
    total_levels<-n_levels1*n_levels2*n_levels3
    min_factor<-min(c(n_levels1,n_levels2,n_levels3))
  }
  
  # Positional vector
  pos_vector<-seq(1,total_levels)
  
  # Modify the vector according the grouping
  for (i in seq(1,total_levels)) {
    if (i%in%c(seq(n_levels1,total_levels,by = n_levels1)+1)) {
      pos_vector[i:length(pos_vector)]<-pos_vector[i:length(pos_vector)]+sep_factor
    }
  }
  return(pos_vector)
}

# Get an average position in the plot for each factor level combinations
automatic_spacing_mean_pos<-function(factor1,factor2=NULL,factor3=NULL,sep_factor=1) {
  
  # Number of levels per factors
  n_levels1<-length(levels(factor1))
  n_levels2<-length(levels(factor2))
  n_levels3<-length(levels(factor3))
  
  # Parse optional arguments
  if (is.null(factor3)&is.null(factor2)) {
    total_levels<-n_levels1
    min_factor<-n_levels1
  } else if (is.null(factor3)) {
    total_levels<-n_levels1*n_levels2
    min_factor<-min(c(n_levels1,n_levels2))
  } else {
    total_levels<-n_levels1*n_levels2*n_levels3
    min_factor<-min(c(n_levels1,n_levels2,n_levels3))
  }
  
  # Positional vector
  pos_vector<-seq(1,total_levels)
  
  # Modify the vector according the grouping
  for (i in seq(1,total_levels)) {
    if (i%in%c(seq(n_levels1,total_levels,by = n_levels1)+1)) {
      pos_vector[i:length(pos_vector)]<-pos_vector[i:length(pos_vector)]+sep_factor
    }
  }
  
  # Calculate mean positions
  mean_pos<-unlist(lapply(split(x = pos_vector,f = ceiling(seq_along(pos_vector)/min_factor)),mean),use.names = F)
  
  return(mean_pos)
  
}

# Multifactorial (max.2) color palette
autocolors<-function(factor1,factor2=NULL,factor3=NULL,intercalate=FALSE) {
  
  # MY COLOR PALETTES (6 levels in factor1 and 4 levels i factor2):
  mycolors<-c("lightblue","lightgoldenrod","lightsalmon","olivedrab","wheat","slateblue")
  mycolors1<-c("lightblue1","lightgoldenrod1","lightsalmon1","olivedrab1","wheat1","slateblue1")
  mycolors2<-c("lightblue2","lightgoldenrod2","lightsalmon2","olivedrab2","wheat2","slateblue2")
  mycolors3<-c("lightblue3","lightgoldenrod3","lightsalmon3","olivedrab3","wheat3","slateblue3")
  mycolors4<-c("lightblue4","lightgoldenrod4","lightsalmon4","olivedrab4","wheat4","slateblue4")
  
  # Parse optional arguments
  if (is.null(factor2)&is.null(factor3)) {
    n_levels_f1<-length(levels(factor1))
    colors<-mycolors[1:n_levels_f1]
  } else if (is.null(factor3)) {
    n_levels_f1<-length(levels(factor1))
    n_levels_f2<-length(levels(factor2))
    multicolor<-cbind(mycolors[1:n_levels_f1],
                      mycolors1[1:n_levels_f1],
                      mycolors2[1:n_levels_f1],
                      mycolors3[1:n_levels_f1],
                      mycolors4[1:n_levels_f1])
    colors<-c(multicolor[,1:n_levels_f2])
  } else {
    n_levels_f1<-length(levels(factor1))
    n_levels_f2<-length(levels(factor2))
    n_levels_f3<-length(levels(factor3))
    multicolor<-cbind(mycolors[1:n_levels_f1],
                      mycolors1[1:n_levels_f1],
                      mycolors2[1:n_levels_f1],
                      mycolors3[1:n_levels_f1],
                      mycolors4[1:n_levels_f1])
    colors<-c(multicolor[,1:n_levels_f2])
    
    if (intercalate==TRUE) {
      colors<-rep(colors,
                  each=1,
                  times=n_levels_f3)
    } else {
      colors<-rep(colors,
                  each=n_levels_f3,
                  times=1)
    }
    
  }
  
  return(colors)
}

# Get coordinates por ploting data points over the plot. With jittering
autopoints<-function(x,factors,pos,proportional=FALSE,jitter_factor=1) {
  
  # Check is numeric the x and a list the factors
  if (!is.numeric(x)) {
    print("x must be numeric")
    stop()
  } else if (class(factors)!="list" & class(factors)!="factor") {
    print("factors must be a factor or a list of factors")
    stop()
  }
  
  # Get values per factor levels
  values<-split(x = x,f = factors)
  
  # Collect the coordinates in a matrix
  myjitter<-c()
  for(i in 1:length(values)) {
    level_proportion<-length(values[[i]])/length(x)
    myjitter<-c(myjitter,jitter(x = rep(pos[i],length(values[[i]])), factor = jitter_factor, amount = if(proportional) {amount = as.vector(t(level_proportion))/2}))
  }
  
  ## Return the coordinates
  return(cbind(myjitter,unlist(values)))
}


# Modification of Tukey Post-Hoc test included in Agricolae package
## Modify the order of the grouping letter, so they match the factor level's order
myHSD.test<-function(aov,factor,alpha=0.05) {
  
  # Check if multcompView is loaded
  if (!"multcompView" %in% (.packages())) {
    try(library("multcompView"))
  }
  
  # Check aov argument is oav class
  if (!"aov"%in%class(aov)) {
    print("aov argument is not of 'aov' class")
    stop()
  }
  
  # Levels order:
  order<-levels(aov$model[,2])
  # Perform the Tukey test
  T<-TukeyHSD(x = aov,which = factor,conf.level = 1-alpha)
  # Obtain the p-values of Tukey in form of a matrix
  p_matrix<-vec2mat(T[[1]][,4],sep = "-")
  # Order the p values matrix according to specified order
  p_matrix<-p_matrix[order,order]
  # Obtain the letter vector
  letras<-multcompLetters(p_matrix,threshold = 0.05)
  # Obtain mean values from lm
  means<-tapply(aov$model[,1],aov$model[,2],mean)
  # Bind means and letters
  table<-data.frame(Means=means,Letters=letras$Letters[match(names(letras$Letters),names(means))])
  return(table)
  
}


# Modification of myHSD function (see above) for oav object created with sufficient statistics

## Modify the order of the grouping letter, so they match the factor level's order
myHSD.test.Sufficient<-function(aov,factor,alpha=0.05) {
    
    # Check if multcompView is loaded
    if (!"multcompView" %in% (.packages())) {
        try(library("multcompView"))
    }
    
    # Check aov argument is oav class
    if (!"aov"%in%class(aov)) {
        print("aov argument is not of 'aov' class")
        stop()
    }
    
    # Levels order:
    order<-aov$model[,2]
    # Perform the Tukey test
    T<-TukeyHSD(x = aov,which = factor,conf.level = 1-alpha)
    # Obtain the p-values of Tukey in form of a matrix
    p_matrix<-vec2mat(T[[1]][,4],sep = "-")
    # Order the p values matrix according to specified order
    p_matrix<-p_matrix[order,order]
    # Obtain the letter vector
    letras<-multcompLetters(p_matrix,threshold = 0.05)
    # Obtain mean values from lm
    means<-aov$model[,1]
    names(means)<-aov$model[,2]
    # Bind means and letters
    table<-data.frame(Means=means,Letters=letras$Letters[match(names(letras$Letters),names(means))])
    return(table)
    
}
