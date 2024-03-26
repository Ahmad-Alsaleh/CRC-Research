"
This code is used to load the datasets and combine them together
One excel files is exported: dataset.xlsx which contains the
actual data of patients, i.e.: their demographics
(Age, Sex, and Diagnoses) and OTU values

Spring 2023
Ahmad Alsaleh

"

# libraries ----
rm(list = ls())
library(readxl)
library(xlsx)
library(dplyr)

# importing datasets ----

# first dataset - Metadata ---

# importing dataset
metadata <- read_excel(
  "Original Datasets/Bataineh Metadata CRC 2020 - CRC vs. Normal.xlsx",
  col_types = c("text", "numeric", "text", "text", "skip", "skip")
) %>% as.data.frame()

metadata$Sex <- as.factor(metadata$Sex)
metadata$Diagnosis <- as.factor(metadata$Diagnosis)

# removing duplicates
metadata <- metadata[!duplicated(metadata[, 1]), ]

# second dataset - OTU table ---

# importing dataset
otu_table <- read.delim2("Original Datasets/OTU_table.txt", header = T, row.names = 1, quote = "")

# transposing OTU table
otu_table <- as.data.frame(t(otu_table))
rownames(otu_table) <- otu_table[, 1]
otu_table <- otu_table[, -1]

# removing duplicates
otu_table <- otu_table[!duplicated(rownames(otu_table)), ]

# removing taxonomy (can be obtained from the .biom file)
otu_table <- otu_table[-nrow(otu_table), ]

# converting OTU features to numeric
row_names <- rownames(otu_table)
otu_table <- as.data.frame(apply(otu_table, 2, as.numeric))
rownames(otu_table) <- row_names
rm(row_names)

# combining the two datasets ----

# replacing '-' and '/' by '.' in  `Accession No` for consistency b/w
# between the two datasets
metadata$`Accession No` <- sub("[-/]", ".", metadata$`Accession No`)

# sorting the datasets for easier mapping
otu_table <- otu_table[order(rownames(otu_table)), ]
metadata <- metadata[order(metadata$`Accession No`), ]
rownames(metadata) <- seq(1, nrow(metadata))

# mapping the intersection of the two datasets
# (ignoring observations not in both datasets)
common_accession_no <- intersect(rownames(otu_table), metadata$`Accession No`)
combined_data <- cbind(
  metadata[metadata$`Accession No` %in% common_accession_no, ],
  otu_table[common_accession_no, ]
)

# sorting columns (OTUs)
otus <- colnames(combined_data)[-(1:4)]
otus <- as.numeric(sub("Otu", "", otus))
otus <- c(-3:0, otus)
combined_data <- combined_data[, order(otus)]

# changing row names
combined_data <- textshape::column_to_rownames(combined_data, loc = 1)

# exporting datasets to excel files ----
write.xlsx(combined_data, "Output Datasets/dataset.xlsx")
