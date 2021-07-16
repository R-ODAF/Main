######################################################
### DEseq2 for normalization and DE identification ###
######################################################

library(DESeq2)
require("DESeq2")
require("edgeR")

###################################################################################
###################################################################################
# PARAMETERS TO SET MANUALLY                            		 

# Set file locations				
sampledir <- "/Path/to/data/TestData/"
outputdir <- paste(sampledir, "Output/", sep="")
if(!dir.exists(outputdir)) {dir.create(outputdir)}
# Names of files to load
SampleDataFile <- "sampleData.csv" #This comma delimited file contains the merged RSEM.genes.results files
SampleKeyFile <- "samplekeyTEST.csv" #This comma delimited file contains at least 2 columns: NAME (sample names identical to the column names of sampleData) and Compound (needs to identify to which group the sample belongs -> ExperimentalGroup & ControlGroup)

# Specify which groups need to be compared 
Samp4compare<- c("PIR") # Experimental group, needs to correspond with the SampleKeyFile. Several condition can be tested c("test1","test2), with the corresponding control
Cont4compare<- c("Vehicle") # Control group, needs to correspond with the SampleKeyFile
DESIGN<- "Compound"	#Column name samplekeyTEST.csv which defines the groups to be compared

# Set analysis ID. This ID will be used as prefix for the output files
analysisID <-"R-ODAF_test_PIR"
# Specify used platform/technology for data generation:
Platform <- "RNA-Seq" # Specify "RNA-Seq" or "TempO-seq"


###################################################################################
###################################################################################

#Set parameters according to platform
if (Platform=="RNA-Seq"){
  minCoverage <- 5000000
  MinCount<- 1
  pAdjValue<- 0.01 
} else if (Platform=="TempO-seq") {
  minCoverage <- 1000000
  MinCount<- 0.5
  pAdjValue<- 0.05 
} else { print("Platform/technology not recognized") }


# Load input files 
setwd(sampledir)
sampleData <- read.delim(SampleDataFile, sep=",", stringsAsFactors=FALSE, header=TRUE,  quote="\"", row.names=1)
DESeqDesign <- read.delim(SampleKeyFile, stringsAsFactors=FALSE, sep=",", header=TRUE,  quote="\"", row.names="NAME")

NORM_TYPE<-paste0(analysisID, "_DESeq2_", Platform)
print(NORM_TYPE)

# First data clean-up: replace NA & remove samples with total readcount < threshold
sampleData[ is.na(sampleData) ] <- 0 
sampleData<- sampleData[,(colSums(sampleData)> minCoverage)]
DESeqDesign <- DESeqDesign[rownames(DESeqDesign) %in% colnames(sampleData),]

# Generate a PCA plot to exclude the post-processing outliers replicates
# vst <- vst(dds)
# plotPCA(vst,intgroup=c(DESIGN))
# exclude outlier from the dataset 

##########
# DESeq2 #
##########

for (x in 1:length(Samp4compare)){	## for all comparisons to be done	
  condition1<- Cont4compare[x]	    		
  condition2<- Samp4compare[x]  
  
  DE_Design <- matrix(data=NA, ncol=2)
  DE_Design <- DESeqDesign [c(grep(condition1,DESeqDesign[,DESIGN[x]]), grep(condition2,DESeqDesign[,DESIGN[x]])),]
  samples <- sampleData[, rownames(DE_Design) ]
  
  ###########
  print(paste(condition2, " vs ", condition1, ":", NORM_TYPE))		
  
  colnames(samples)<-NULL
  dds <- DESeqDataSetFromMatrix(countData = round(samples), colData = as.data.frame(DE_Design), design = as.formula(paste0("~", DESIGN[x])))
  
  print("Wait... (dds step executing)")				
  dds <- DESeq(dds, quiet=TRUE)
  
  print(paste0("Filtering genes: 75% of at least 1 group need to be above ", MinCount, " CPM"))
  
  
  #Filter low readcounts (genes not meeting the condition : at least one condition with 75% of the samples above 1 CPM)
  
  SampPerGroup<-table(DE_Design[,DESIGN])
  Counts<-counts(dds, normalized=TRUE)
  CPMdds<-cpm(counts(dds, normalized=TRUE))
  
  Filter <- matrix(data=NA, ncol=3, nrow= nrow(Counts))
  rownames(Filter) <- rownames(Counts)
  colnames(Filter) <- c("Low","quantile","spike")
  
  # Apply the "Relevance" condition
  
  for (gene in 1:nrow(dds)) {
    
    CountsPass<-NULL
    for (group in 1:length(SampPerGroup)) { 
      sampleCols<-grep(dimnames(SampPerGroup)[[1]][group],DE_Design[,DESIGN])
      Check<-sum(CPMdds[gene,sampleCols] >= MinCount)>= 0.75*SampPerGroup[group]
      CountsPass<-c(CountsPass, Check)
    }
    
    if ( sum(CountsPass) > 0 ) {Filter[gene,1] <- 1 }	else { Filter[gene,1] <- 0 }
    
  }	
  
  compte <- Counts[Filter[,1] == 1,]
  Filter <- Filter[rownames(Filter) %in% rownames(compte),]
  
  print(paste("Relevance filtering removed ",nrow(dds)- nrow(Filter)," genes from the ",nrow(dds)," assessed. ", nrow(Filter)," genes remaining",sep=""))
  
  
  print("Obtaining the DESeq2 results")
  
  # compute the DEGs on the genes passing the Relevance condition
  
  res <- results(dds[rownames(compte),], alpha = 0.01, cooksCutoff =F, independentFiltering=F, contrast=c(DESIGN[x], condition2, condition1), pAdjustMethod= 'fdr')
  
  setwd(outputdir)
  FileName<-paste(NORM_TYPE, condition2,"vs",condition1, "FDR", pAdjValue, sep="_")
  
  #Save output tables		
  norm_data <<- counts(dds[rownames(compte)],normalized=TRUE) 
  DEsamples <<- subset(res,res$padj < pAdjValue)	
  
  DECounts <- compte[rownames(compte) %in% rownames(DEsamples),]
  Filter <- Filter[rownames(Filter) %in% rownames(DECounts),]
  
  print("Check median against third quantile" )
  print("AND")
  print("Check the presence of a spike" )
  
 if (NROW(DECounts) > 0) {
  for (gene in 1:NROW(DECounts)) {
    
    # Check the median against third quantile
    quantilePass <-NULL
    sampleColsg1 <- grep(dimnames(SampPerGroup)[[1]][1],DE_Design[,DESIGN])
    sampleColsg2 <- grep(dimnames(SampPerGroup)[[1]][2],DE_Design[,DESIGN])
    
    Check <- median(DECounts[gene,sampleColsg1]) > quantile(DECounts[gene,sampleColsg2], 0.75)[[1]]
    quantilePass <-c(quantilePass, Check)
    Check <- median(DECounts[gene,sampleColsg2]) > quantile(DECounts[gene,sampleColsg1], 0.75)[[1]]
    quantilePass <-c(quantilePass, Check)
    
    if ( sum(quantilePass) > 0 ) {Filter[gene,2] <- 1 }	else { Filter[gene,2] <- 0 }
    
    # Check for spike 
    spikePass <- NULL
    for (group in 1:length(SampPerGroup)) { 
      sampleCols<-grep(dimnames(SampPerGroup)[[1]][group],DE_Design[,DESIGN])
      if (max(DECounts[gene,sampleCols]) ==0) {Check <- FALSE} else {
        Check <- (max(DECounts[gene,sampleCols])/sum(DECounts[gene,sampleCols])) >= 1.4*(SampPerGroup[group])^(-0.66)
        spikePass<-c(spikePass, Check)
      }
    }
    if ( sum(spikePass) > 1 ) {Filter[gene,3] <- 0 }	else { Filter[gene,3] <- 1 }
    
  }		
  
  # extract the final list of DEGs
  
  DECounts_real <- DEsamples[rowSums(Filter) == 3 ,]
  DECounts_no_quant <- DEsamples[Filter[,2] == 0 ,]
  DECounts_spike <- DEsamples[Filter[,3] == 0 ,]
}  
  print(paste("A total of ",NROW(DECounts_real), " DEGs were selected, after ",NROW(DECounts_no_quant)," genes(s) removed by the quantile rule and ", NROW(DECounts_spike)," gene(s) with a spike",sep=""))
  
  # save the normalized counts and the list of DEGs
  write.table(norm_data,file=paste0(FileName, "_Norm_Data.txt"), sep="\t", quote=FALSE)
  write.table(DECounts_real,file=paste0(FileName,"_DEG_table.txt"), sep="\t", quote=FALSE)
  # save the filtered DEGs by either the quantile or the Spike rules
  write.table(DECounts_no_quant,file=paste0(FileName, "_quantile_filt.txt"), sep="\t", quote=FALSE)
  write.table(DECounts_spike,file=paste0(FileName,"_quantile_filt.txt"), sep="\t", quote=FALSE)
  
  print("DESeq2 Done")
  
}


