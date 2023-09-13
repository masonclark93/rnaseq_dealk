
library(tidyverse)


# Read the GMT file into a data frame
gmt_file <- "/Users/masonclark/Dropbox/Mac/Desktop/rnaseq_dealk/Data/RNA_seq/sequence_data/my_analysis/new_r_analysis_071322/normalized_data_analyses/GSEA_BROAD/Hzea_GOALL.gmt"
gmt_data <- read_delim(gmt_file, delim = "\t", col_names = FALSE)

# Function to aggregate genes in each gene set
aggregate_genes <- function(x) {
  paste(unique(unlist(x)), collapse = "\t")
}

# Group by the first column (gene set) and aggregate genes
aggregated_data <- gmt_data %>%
  group_by(X1) %>%
  summarize(across(starts_with("X"), aggregate_genes)) %>%
  ungroup()



#Write to file
write.table(aggregated_data, "Hzea_GOALL_FINAL_072723.gmt", col.names = FALSE)
