# R-ODAF
Omics Data Analysis Framework for Regulatory application

  
Dear future R-ODAF user,

The CEFIC project has developed Omics Data Analysis Frameworks for Regulatory application (R-ODAFs) in order to supply robust and reproducible data analysis pipelines for transcriptomics analysis. While developed for regulatory purposes, the pipelines can also be applied for research purposes.

This GitHub contains scripts/files nessesary for the sequencing R-ODAF (developped for both RNA-Seq and adaptable for TempO-seq technology).
You can either: 
- run the docker container (/R-ODAF_docker). The read_me file in the directory will guide you in the different step to follow
- independently run the "R-ODAF_sequencing_DataPreprocessing.sh" followed by "R_ODAF_DEGs.R" from the "Scripts" directory. For this, you will need to install all required package used in the R-ODAF(fastp, multiqc, STAR, RSEM, DESeq2)


The scripts/files necessary for the Microarray R-ODAF can be found in the ArrayAnalysis repository: https://github.com/arrayanalysis/

Thank you for using the R-ODAF!

Kind Regards,
The CEFIC C4 team
