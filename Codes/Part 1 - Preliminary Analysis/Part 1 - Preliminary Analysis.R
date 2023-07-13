"
This code analyzes microbial data by performing data preprocessing,
statistical tests, and visualization. It reads an Excel file containing
microbial data, and calculates summary statistics. It conducts univariate
tests and filters data based on the Wilcox test. It saves the processed data
to a new Excel file and analyzes taxonomy data using specialized libraries.
The script aims to provide insights into the relationships between variables
and identify significant differences between groups in microbial data.

"

# libraries ----
set.seed(42)
rm(list = ls())
source("Codes/Part 1 - Preliminary Analysis/utility.R")
library(biomformat)
library(qiimer)
library(stringi)
library(readxl)
library(taRifx)
library(summarytools)
library(ggplot2)
library(xlsx)
library(vegan)
library(dplyr)
library(phyloseq)
library(stringr)

# loading data set ----
data = as.data.frame(read_excel("Output Datasets/dataset.xlsx"))
data = textshape::column_to_rownames(data)

data$Sex = as.factor(data$Sex)
data$Diagnosis = as.factor(data$Diagnosis)

data = japply(data, which(sapply(data, is.numeric)), as.integer)

meta.data = data[1:3]
otu.data = data[-(1:3)]

# summary statistics ----
summarytools::view(stby(meta.data$Age, meta.data$Diagnosis, descr,transpose = T,
												stats = c("mean", "sd")))
summarytools::view(stby(meta.data$Age, meta.data$Sex, descr, transpose = T,
												stats = c("mean", "sd")))
summarytools::view(dfSummary(meta.data, plain.ascii = F, style = "grid",
														 graph.magnif = 0.75, valid.col = F, tmp.img.dir = "/tmp"))

# univariate tests on metadata ----

# Age vs. Diagnosis
p_value = wilcox.test(meta.data$Age ~ meta.data$Diagnosis, exact = F)$p.value
p_value

if (p_value <= 0.05) {
	print("Age and Diagnosis are different")
} else {
	print("Age and Diagnosis are not different")
}

# Gender vs. Diagnosis
p_value = fisher.test(meta.data$Sex, meta.data$Diagnosis)$p.value
p_value
if (p_value <= 0.05) {
	print("Sex and Diagnosis are different")
} else {
	print("Sex and Diagnosis are not different")
}

# Age vs. Gender
p_value = wilcox.test(meta.data$Age ~ meta.data$Sex, exact = F)$p.value
p_value
if (p_value <= 0.05) {
	print("Age and Gender are different")
} else {
	print("Age and Gender are not different")
}

rm(p_value)

# converting counts to proportions ----
otu.data = otu.data / rowSums(otu.data) * 100

'
# removing near zero variance variables ----
nzv.vars = caret::nearZeroVar(otu.data)
otu.data[, -nzv.vars]
rm(nzv.vars)
'
# finding number of OTU's per group ----
otu.sum = rowsum(otu.data, meta.data$Diagnosis)

crc.otus = which(otu.sum["CRC", ] != 0)
normal.otus = which(otu.sum["Normal", ] != 0)

length(crc.otus)
length(normal.otus)
length(intersect(crc.otus, normal.otus)) / ncol(otu.data)

rm(crc.otus, normal.otus, otu.sum)

# filtering data using Wilcox Test ----
p_values = sapply(otu.data, function(col)
	wilcox.test(col ~ meta.data$Diagnosis, exact = F)$p.value)

p_values = sort(p_values, decreasing = T)

# using a scree plot to choose a cut-off (at the elbow)
scree.plot = ggplot(data.frame(Index = 1:length(p_values), p_value = p_values),
										aes(x = Index, y = p_value)) + geom_line() + labs(x = "OTU Index")
scree.plot

cutoff = 955
scree.plot + geom_vline(xintercept = cutoff, linetype = "dashed", color = "red") +
	geom_text(x = cutoff, y = 0.1, hjust = -0.1, label = paste0("Index = ", cutoff), color = "red")

selected.features = names(p_values[(cutoff + 1):length(p_values)])
otu.data = otu.data[, selected.features]

rm(scree.plot, selected.features, cutoff, p_values)

# adding response to the last column
otu.data = data.frame(otu.data, Diagnosis = meta.data$Diagnosis)

# shuffling observations
otu.data = otu.data[sample(nrow(otu.data)), ]

# exporting analysis dataset ----
# this is the cleaned data set that to be used for further analysis.
write.xlsx(otu.data, "Output Datasets/Analysis Dataset.xlsx")

# Ecological Assessment Methods ----

# reading taxonomy data
biom.data = read_biom("Original Datasets/BIOM Files/OTU_table.biom")
taxonomy.data = biom_taxonomy(biom.data)
tax = as.data.frame(t(stri_list2matrix(taxonomy.data)))
rownames(tax) = names(taxonomy.data)
colnames(tax) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# adding prefix to each taxa level
is.na(tax) <- tax == "Unclassified"
tax$Kingdom = str_c('K_', tax$Kingdom)
tax$Phylum = str_c('P_', tax$Phylum)
tax$Class = str_c('C_', tax$Class)
tax$Order = str_c('O_', tax$Order)
tax$Family = str_c('F_', tax$Family)
tax$Genus = str_c('G_', tax$Genus)
tax$Species = str_c('S_', tax$Species)

# handling unclassified (NA) values
taxonomy.data = handle.unclassified.values(tax)

rm(tax, biom.data)

write.xlsx(taxonomy.data, "Output Datasets/taxonomy.xlsx")

# creating phyloseq object
OTU = t(data[, -(1:3)]) %>% as.matrix %>% otu_table(taxa_are_rows = T)
TAX = taxonomy.data %>% as.matrix %>% tax_table

phy.seq = phyloseq(OTU, TAX)
sample_data(phy.seq) = meta.data

# pruning unnecessary columns
phy.seq = prune_taxa(taxa_sums(phy.seq) > 0, phy.seq)

# creating LEfSe dataset for Galaxy platform
lefse.dataset = create.lefse.dataset(phy.seq)

# exporting LEfSe dataset
lines <- c(paste(c('ID', names(lefse.dataset)), collapse = '\t'),
					 sapply(seq_len(nrow(lefse.dataset)), function(i)
					 	paste(c(row.names(lefse.dataset)[i],
					 					lefse.dataset[i, ]), collapse = '\t')))

writeLines(lines, con = "Output Datasets/LEfSe_dataset.txt")

# alpha diversity ----
plot = plot_richness(phy.seq, x = "Diagnosis", color = "Sex",
										 measures = c("Chao1", "Shannon", "simpson", "ace"))
plot = plot + geom_boxplot(data = plot$data, aes(x = Diagnosis, y = value,
																								 color = "Sex"), alpha = 0.1)
plot

richness = estimate_richness(phy.seq, split = T, measures =
														 	c("Chao1", "Shannon", "simpson", "ace"))
richness = cbind(richness, Diagnosis = sample_data(phy.seq)$Diagnosis)
richness$se.chao1 = NULL
richness$se.ACE = NULL

richness.p.values = sapply(richness[, -ncol(richness)], function(x)
	wilcox.test(x ~ richness$Diagnosis, exact = F)$p.value)
p.adjust(richness.p.values, method = "fdr")

# beta diversity ----
phy.seq = transform_sample_counts(phy.seq, function(x) x / sum(x) * 100)
bray.dist = distance(phy.seq, method = "bray")
anosim(bray.dist, sample_data(phy.seq)$Diagnosis)

# creating abundance plots of most abundant OTUs by phylum and genus levels ----
phy.seq.merged = merge_samples(phy.seq, group = "Diagnosis")
sample_data(phy.seq.merged)$Diagnosis = levels(sample_data(phy.seq)$Diagnosis)
phy.seq.merged = transform_sample_counts(phy.seq.merged, function(x) x / sum(x) * 100)
top.20.otus = names(sort(taxa_sums(phy.seq), decreasing = T)[1:20])
phy.seq.merged.top.20 = prune_taxa(top.20.otus, phy.seq.merged)

# Phylum-level
plot = plot_bar(phy.seq.merged.top.20, "Diagnosis", fill = "Phylum") + coord_flip() +
	ylab("Relative Abundance (%)") + labs(fill = "Abundant Phylum-level OTU")
plot

# comparing abundance values
abundances = data.frame(
	Diagnosis = as.factor(plot$data$Diagnosis),
	Bacteria = as.factor(plot$data$Phylum),
	Abundance = plot$data$Abundance
)
abundances = abundances %>%
	group_by(Diagnosis, Bacteria) %>%
	summarize(Abundance = sum(Abundance)) %>%
	ungroup %>% as.data.frame

abundances.crc = abundances[abundances$Diagnosis == "CRC", ]
abundances.normal = abundances[abundances$Diagnosis == "Normal", ]
abundances.table = data.frame(
	Abundances.CRC = abundances.crc$Abundance,
	Abundances.Normal = abundances.normal$Abundance,
	row.names = abundances.crc$Bacteria
)
abundances.table

# sorting by abundances
abundances.crc = abundances.crc[order(-abundances.crc$Abundance), ]
abundances.normal = abundances.normal[order(-abundances.normal$Abundance), ]
cbind(abundances.crc, abundances.normal)

# grouping by bacteria
abundances.by.bacteria = abundances %>% group_by(Bacteria) %>% summarize(Abundance = sum(Abundance)) %>% as.data.frame
abundances.by.bacteria = abundances.by.bacteria[order(-abundances.by.bacteria$Abundance), ]
abundances.by.bacteria

# grouping by Diagnosis
abundances.by.diagnosis = abundances %>% group_by(Diagnosis) %>% summarize(Abundance = sum(Abundance)) %>% as.data.frame
abundances.by.diagnosis = abundances.by.diagnosis[order(-abundances.by.diagnosis$Abundance), ]
abundances.by.diagnosis

# Genus-level
plot = plot_bar(phy.seq.merged.top.20, "Diagnosis", fill = "Genus") + coord_flip() +
	ylab("Relative Abundance (%)") + labs(fill = "Abundant Genus-level OTU")
plot

abundances = data.frame(
	Diagnosis = as.factor(plot$data$Diagnosis),
	Bacteria = as.factor(plot$data$Genus),
	Abundance = plot$data$Abundance
)
abundances = abundances %>%
	group_by(Diagnosis, Bacteria) %>%
	summarize(Abundance = sum(Abundance)) %>%
	ungroup %>% as.data.frame

abundances.crc = abundances[abundances$Diagnosis == "CRC", ]
abundances.normal = abundances[abundances$Diagnosis == "Normal", ]
abundances.table = data.frame(
	Abundances.CRC = abundances.crc$Abundance,
	Abundances.Normal = abundances.normal$Abundance,
	row.names = abundances.crc$Bacteria
)
abundances.table

# sorting by abundances
abundances.crc = abundances.crc[order(-abundances.crc$Abundance), ]
abundances.normal = abundances.normal[order(-abundances.normal$Abundance), ]
cbind(abundances.crc, abundances.normal)

# grouping by bacteria
abundances.by.bacteria = abundances %>% group_by(Bacteria) %>% summarize(Abundance = sum(Abundance)) %>% as.data.frame
abundances.by.bacteria = abundances.by.bacteria[order(-abundances.by.bacteria$Abundance), ]
abundances.by.bacteria

# grouping by Diagnosis
abundances.by.diagnosis = abundances %>% group_by(Diagnosis) %>% summarize(Abundance = sum(Abundance)) %>% as.data.frame
abundances.by.diagnosis = abundances.by.diagnosis[order(-abundances.by.diagnosis$Abundance), ]
abundances.by.diagnosis

# finding difference between in abundance for same bacteria
abundances.difference = abundances.table$Abundances.CRC - abundances.table$Abundances.Normal
abundances.difference = data.frame(Bacteria = rownames(abundances.table),
																	 Abundances.Difference = abundances.difference)
abundances.difference = abundances.difference[order(-abs(abundances.difference$Abundances.Difference)), ]
abundances.difference
