mypfamenrichment<-function(gene_universe,mysubset,description,annotation="Arabidopsis",write_tables=FALSE) {
    
    # Is Pfam annotation file loaded? If not load and parse ####
    if (!exists("mypfam")) {
        print("Warning: Pfam annotation file not detected. It will be loaded for GitHub repo")
        Pfam.A.hmm <- read.delim(file = "https://raw.githubusercontent.com/fle1/canolab_scripts/master/useful_data/Pfam-A.hmm.dat" ,header=FALSE, comment.char="/", stringsAsFactors=FALSE)
        
        # Only gather compulsory fields (info from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/userman.txt)
        mypfam<-data.frame(Accession=Pfam.A.hmm$V1[grep(pattern = "#=GF AC",x = Pfam.A.hmm$V1)],
                           ID=Pfam.A.hmm$V1[grep(pattern = "#=GF ID",x = Pfam.A.hmm$V1)],
                           Description=Pfam.A.hmm$V1[grep(pattern = "#=GF DE",x = Pfam.A.hmm$V1)],
                           Type=Pfam.A.hmm$V1[grep(pattern = "#=GF TP",x = Pfam.A.hmm$V1)])
        # Remove comment fields
        mypfam<-apply(mypfam,2,function(x){x<-substr(x,8,nchar(x))})
        mypfam<-data.frame(mypfam)
        # Remove the PFxxxxx.X in the Accession field (Does this mean something?)
        mypfam$Accession<-as.character(mypfam$Accession)
        # head(nchar(mypfam$Accession)) # Very strange!!! This is because the apply function?
        mypfam$Accession<-substr(mypfam$Accession,4,10)
        
        mypfam<<-mypfam
    }
    
    # Is Arabidopsis Annotation loaded
    if (annotation=="Arabidopsis") {
        print("Default annotation is Araport11, it will be used")
        if (!exists("Araport11_annotation")) {
            print("No Arabidopsis annotation file is detected in the environment...")
            user<-readline(prompt = "Do you want to load Araport11 annotation? (y/n):   ")
            user<-tolower(user)
            if (user=="y") {
                # Original file from phytozome
                Araport11_annotation<<-read.delim(file = "https://raw.githubusercontent.com/fle1/canolab_scripts/master/useful_data/Araport11.annotation_info.txt",stringsAsFactors = TRUE)
                annotation<-data.frame(Araport11_annotation)
            } else {
                stop()
            }
        } else {
            print("Araport11 annotation detected in the environment...")
            annotation<-Araport11_annotation
        }
        
    }
    
    # Reduce the number of genes provided to only genes with Pfam annotation
    my_new_subset<-annotation[grep("PF",annotation$Pfam),c("locusName","Pfam")]
    my_new_subset<-my_new_subset[which(my_new_subset$locusName%in%mysubset),"locusName"]
    my_new_subset<-unique(as.character(my_new_subset))
    
    # Create the gene universe based in the number of genes identified in the RNAseq with Pfam ID annotations
    my_gene_universe<-annotation[grep("PF",annotation$Pfam),c("locusName","Pfam")]
    my_gene_universe<-my_gene_universe[which(my_gene_universe$locusName%in%gene_universe),]
    
    # Reduce redundancy, put together same genes (isoforms)
    if (mean(nchar(as.character(my_gene_universe$locusName))) > 9) {
        print("isoforms annotation for Arabidopsis will be collapsed")
        my_gene_universe<-aggregate(Pfam~locusName,my_gene_universe,paste,collapse=",")
    } else {
        
    }
    
    # Give some info:
    print(paste("Your subset had",
                length(mysubset),
                "entries but only",
                length(unique(my_new_subset)),
                "have Pfam annotations.",
                "Only annotated genes will be used in the enrichment",
                sep = " ")
    )
    
    print(paste("Your gene universe has",
                length(gene_universe),
                "entries but only",
                length(unique(my_gene_universe$locusName)),
                "have Pfam annotations.",
                "Only annotated genes will be used in the enrichment",
                sep = " ")
    )
    
    
    # Unlist mygene universe data.frame (test only genes present in my subset)
    # Make no sense to test all and make the function very slow
    mydomains<-unique(unlist(
        strsplit(
            as.character(
                my_gene_universe[match(my_new_subset,my_gene_universe$locusName),"Pfam"]),
            split = ",")
        )
    )
    if (length(which(is.na(mydomains)))>0) {
        mydomains<-mydomains[-which(is.na(mydomains))]  #I expect only one NA, coming from non-matched
    }
    mydomains<-as.character(mydomains)
    
    
    # Loop for each Pfam ID (Core of the enrichment function)
    pvector<-c()
    for (i in 1:length(mydomains)){
        #L1 count my Pfam ID in DE genes
        L1<-my_gene_universe[grep(pattern = mydomains[i],my_gene_universe$Pfam),1]
        L1<-length(which(L1%in%my_new_subset))
        #L2 Count my Pfam ID in the NonDE genes
        L2<-my_gene_universe[grep(pattern = mydomains[i],my_gene_universe$Pfam),1]
        L2<-length(which(L2%in%setdiff(gene_universe,my_new_subset)))
        #L3 Count the Non myPfam ID in the DE genes
        L3<-my_gene_universe[-grep(pattern = mydomains[i],my_gene_universe$Pfam),1]
        L3<-length(which(L3%in%my_new_subset))
        #L4 Count the Non myPfam ID in the NonDE genes
        L4<-my_gene_universe[-grep(pattern = mydomains[i],my_gene_universe$Pfam),1]
        L4<-length(which(L4%in%setdiff(gene_universe,my_new_subset)))
        # Create the contingency table
        contig<-cbind(c(L1,L3),c(L2,L4))
        colnames(contig)<-c("DE","NonDE")
        rownames(contig)<-c("MyPfam","NoMyPfam")
        #Test Fisher/Hypergeometric
        pvector[i]<-fisher.test(contig,alternative = "greater")$p.value
        names(pvector)[i]<-mydomains[i]
    }
    
    # Generate a data.frame containing all info. Plot p-values
    # Order and adjust the p-vector
    pvector<-p.adjust(pvector,method = "BH")
    pvector<-pvector[order(pvector,decreasing = FALSE)]
    # Plot p-values distribution
    plot(-log(pvector),type="l",main = "p-value distribution",ylab = "-log(p-value)")
    abline(h=-log(0.05),lty=2,col="tomato")
    # Find enriched Pfam in annotation
    finalmatrix<-cbind(mypfam[match(names(pvector),mypfam$Accession),],
                       pvalues_BH=pvector)
    print(finalmatrix[1:10,])
    # Find the genes in the subset Write the tables (Optional)
    mygenesannotated<-annotation[match(my_new_subset,annotation$locusName),c("locusName","Pfam")]
    finalmatrix<-data.frame(finalmatrix,Genes_in_subset=rep(NA,length(pvector)))
    
    for (i in 1:length(pvector)){
        finalmatrix$Genes_in_subset[i]<-paste(
            as.vector(
                mygenesannotated[
                    grep(
                        names(
                            pvector)[i],
                        mygenesannotated$Pfam),
                    "locusName"]
            ),
            collapse = ",")
    }
    
    return(finalmatrix)
    
    # Write tables (Optional)
    if (write_tables==TRUE) {
        write.table(finalmatrix,
                    paste("Pfam_enrichment",paste(description,"txt",sep="."),sep = "_"),
                    quote = FALSE,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep = "\t")
        print(paste("Tables were written on:",getwd(),sep=" "))
    }
    
    else print("Analysis DONE")
    
}
