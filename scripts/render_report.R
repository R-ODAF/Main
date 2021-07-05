#!/usr/bin/R
# Custom parameters for the report
library(tidyverse)
params <- list()
params$projectdir <- here::here() # If you open this as an R project, it should find your root directory
params$project_name <- "GT21"
params$species <- "human"     # one of human, mouse, rat, hamster
params$design <- "group"    # single experimental group of interest
params$nuisance <- NA # "batch"
params$intgroup <- c("group") # Can be multiple columns of interest (covariates, etc)
params$flag <- TRUE         # runs all analysis by default when set to TRUE
params$platform <- "TempO-Seq" # Choose from RNA-Seq or TempO-Seq
params$group_facet <- "chemical" # If you have many different experimental groups, you may subset the report by specifying a column in the metadata to filter groups, and then setting the group of interest in group_filter
params$group_filter <- NA  # If this is set, it will ONLY run the analysis on these groups
params$exclude_samples <- NA  # Optionally, a vector of sample names to exclude from the analysis
params$exclude_groups <-  NA # Optionally, a vector of groups to exclude from the analysis. By default this is assumed to be in the column specified by params$design.
params$include_only_column <- NA # Restrict analysis to group(s) in the column listed here based on params$include_only_group.
params$include_only_group <- NA # Restrict analysis to this/these group(s) within the column listed in params$include_only_column
params$use_cached_RData <- FALSE # If possible, load the saved RData for dds object and gene IDs
params$cpus <- 41 # Set to a lower number (e.g., 2 to 4) if you aren't working in a server environment
params$run_pathway_analysis <- TRUE

skip_extra <- "DMSO" # Remove DMSO controls as a facet

# Input file - Rmd
inputFile <- file.path(params$projectdir, "Rmd", "DESeq2_report.rnaseq.Rmd")

# Identify where metadata can be found
SampleKeyFile <- file.path(params$projectdir,
                           "metadata/metadata.QC_applied.txt")

# Read in metadata
DESeqDesign <- read.delim(SampleKeyFile,
                          stringsAsFactors=FALSE,
                          sep="\t",
                          header=TRUE,
                          quote="\"",
                          row.names=1) # Column must have unique IDs!!
DESeqDesign$original_names <- rownames(DESeqDesign)

# Run DESeq2 and make reports
if(is.na(params$group_facet)){
  message("Writing a single report for whole experiment.")
  # Output file - HTML
  filename <- paste0(params$platform, "_",
                     params$project_name, "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'),
                     ".html")
  outFile <- file.path(params$projectdir,
                       "reports",
                       filename)
  rmarkdown::render(input = inputFile,
                    encoding = encoding,
                    output_file = outFile,
                    params = params,
                    envir = new.env())
} else if(any(!is.na(params$group_filter))) {
  message(paste0("The group(s) of interest is (are) ",
                 paste(params$group_filter, collapse=" and "),".\n",
                 "Writing a single report for that (those) groups."))
  # Output file - HTML
  filename <- paste0(params$platform, "_",
                     params$project_name, "_",
                     paste(params$group_filter, collapse="_"), "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'),
                     ".html")
  outFile <- file.path(params$projectdir,
                       "reports",
                       filename)
  rmarkdown::render(input = inputFile,
                    encoding = encoding,
                    output_file = outFile,
                    params = params,
                    envir = new.env())
} else {
  # Remove params$exclude_groups
  facets <- DESeqDesign %>%
    filter(!(!!sym(params$group_facet)) %in%
             c(params$exclude_groups, skip_extra)) %>%
    pull(params$group_facet) %>% 
    unique()
  message(paste0("Making multiple reports based on ",
                 params$group_facet ,"..."))
  for(i in facets[c(1,12)]){
    message(paste0("Building report for ", i, "..."))
    params$group_filter <- i
    filename <- paste0(params$platform, "_",
                       params$project_name, "_",
                       i, "_",
                       format(Sys.time(),'%d-%m-%Y.%H.%M'),
                       ".html")
    outFile <- file.path(params$projectdir,
                         "reports",
                         filename)
    rmarkdown::render(input = inputFile,
                      encoding = encoding,
                      output_file = outFile,
                      params = params,
                      envir = new.env())
  }
}
