# This is a script containing manual methods that might need required to
# tidy up the data (count table from TempO-Seq experiments)

library(tidyverse)
library(data.table)

rootdir <- here::here()

################################################################################
# Merge multiple separate files - this code is not necessarily complete
################################################################################

# csv_files <- fs::dir_ls("data/processed", regexp = "\\.csv$", recurse=1)

# results <- csv_files %>%
#   map_dfr(read_csv, col_names=T, .id="source") %>%
#   tidyr::separate(source, into = c("delete","Day","flow_cell","delete2"), sep = "\\.") %>%
#   dplyr::select(-c(delete,delete2,flow_cell))

# View(as.data.frame(colnames(results)))

# 23816 by 866

# sampleDataTable <- data.table(results %>% select(-flow_cell,-Day))
# sampleDataTable[is.na(sampleDataTable)] <- 0
# sampleData <- sampleDataTable[, lapply(.SD,sum), by=c("X1")]

# write.table(sampleData,
#             file = file.path("data/processed","genes.data.extra_columns.tsv"),
#             sep=",",
#             quote=F,
#             row.names=F)
# 
# sampleData <- read.delim("data/processed/genes.data.tsv",
#                          sep=",",
#                          stringsAsFactors=FALSE,
#                          header=TRUE, 
#                          quote="\"",
#                          row.names=1,
#                          check.names=FALSE)

################################################################################
# Merge multiple separate files - easier method
################################################################################

# csv_files <- fs::dir_ls("data/processed", regexp = "\\.csv$", recurse=1)
# 
# d1fc1 <- read.csv(csv_files[1], row.names = 1)
# d1fc2 <- read.csv(csv_files[2], row.names = 1)
# d10fc1 <- read.csv(csv_files[3], row.names = 1)
# d10fc2 <- read.csv(csv_files[4], row.names = 1)
# d14fc1 <- read.csv(csv_files[5], row.names = 1)
# d14fc2 <- read.csv(csv_files[6], row.names = 1)
# d4fc1 <- read.csv(csv_files[7], row.names = 1)
# d4fc2 <- read.csv(csv_files[8], row.names = 1)
# 
# d1 <- as.data.frame(d1fc1 + d1fc2)
# d4 <- as.data.frame(d4fc1 + d4fc2)
# d10 <- as.data.frame(d10fc1 + d10fc2)
# d14 <- as.data.frame(d14fc1 + d14fc2)
# 
# results <- cbind(d1,d4,d10,d14)
# 


################################################################################
# Merge a single table across different columns by regex
################################################################################

df <- read.table(file = file.path(rootdir, "data/output/count_table.csv"),
                 sep=",",
                 header=T)

names(df) = gsub(pattern = "_trimmed", replacement = "", x = names(df))

df <- df %>% pivot_longer(cols = -c("X"),
                          names_to = c('.value','flowcell'),
                          names_pattern='(.*)_(.*)')

df[is.na(df)] <- 0
sampleDataTable <- data.table(df %>% dplyr::select(-c(flowcell)))
results <- sampleDataTable[, lapply(.SD,sum), by="X"]

write.table(results,
            file = file.path(rootdir,"data/processed","count_table.csv"),
            sep=",",
            quote=F,
            row.names=F)
