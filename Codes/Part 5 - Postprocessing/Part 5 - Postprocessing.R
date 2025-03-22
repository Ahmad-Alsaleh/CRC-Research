# libraries ----
rm(list = ls())

library(readxl)
library(dplyr)

# Read the features data
features <- read_excel(
  "Output Datasets/Final Selected OTUs.xlsx",
  col_types = c("text", "numeric", "skip", "skip", "skip")
) %>% as.data.frame()

# Get BAM score summary
bam_summary <- features$BAM.Score %>%
  summary()
bam_summary <- data.frame(
  Statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."),
  Value = as.numeric(bam_summary)
)
bam_summary

# Get OTU data summary
otu_data <- read_excel("Output Datasets/Analysis Dataset.xlsx") %>%
  as.data.frame() %>%
  select(features$OTU)

# Create a summary for each OTU
otu_summary <- sapply(otu_data, summary) %>%
  as.data.frame()
otu_summary

# You can also save as CSV files if needed
write.csv(bam_summary, "Output Datasets/final_features_bam_score_summary.csv", row.names = FALSE)
write.csv(otu_summary, "Output Datasets/final_features_otu_summary.csv")
