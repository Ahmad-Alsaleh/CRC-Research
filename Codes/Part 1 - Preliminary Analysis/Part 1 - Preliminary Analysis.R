"
This code analyzes microbial data by performing data preprocessing,
statistical tests, and visualization. It reads an Excel file containing
microbial data, and calculates summary statistics. It conducts univariate
tests and filters data based on the Wilcox test. It saves the processed data
to a new Excel file and analyzes taxonomy data using specialized libraries.
The script aims to provide insights into the relationships between variables
and identify significant differences between groups in microbial data.

"

set.seed(42)

# loading data set ----
library(readxl)
data = as.data.frame(read_excel("Output Datasets/dataset.xlsx"))
data = textshape::column_to_rownames(data)

data$Sex = as.factor(data$Sex)
data$Diagnosis = as.factor(data$Diagnosis)

library(taRifx)
data = japply(data, which(sapply(data, is.numeric)), as.integer)

meta.data = data[1:3]
otu.data = data[-(1:3)]

# summary statistics ----
library(summarytools)
summarytools::view(stby(meta.data$Age, meta.data$Diagnosis, descr, transpose = T, stats = c("mean", "sd")))
summarytools::view(stby(meta.data$Age, meta.data$Sex, descr, transpose = T, stats = c("mean", "sd")))
summarytools::view(dfSummary(meta.data, plain.ascii = F, style = "grid", graph.magnif = 0.75,
							 valid.col = F, tmp.img.dir = "/tmp"))

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
# TODO: I think no need for this (check the paper). if removed, make it in a separate commit (to go back easily)
otu.sum = rowsum(otu.data, meta.data$Diagnosis)

crc.otus = which(otu.sum["CRC", ] != 0)
normal.otus = which(otu.sum["Normal", ] != 0)

length(crc.otus)
length(normal.otus)
length(intersect(crc.otus, normal.otus))

rm(crc.otus, normal.otus, otu.sum)

# filtering data using Wilcox Test ----
p_values = sapply(otu.data, function(col)
	wilcox.test(col ~ meta.data$Diagnosis, exact = F)$p.value)

p_values = sort(p_values, decreasing = T)

# using a scree plot to choose a cut-off (at the elbow)
library(ggplot2)
scree.plot = ggplot(data.frame(index = 1:length(p_values), p_value = p_values),
			 aes(x = index, y = p_value)) +	geom_line()
scree.plot

# TODO: try to run the whole analysis when cutoff = 440 (export the graph with x axis way longer than y)
cutoff = 955

scree.plot + geom_vline(xintercept = cutoff, linetype = "dashed", color = "red") +
	geom_text(x = cutoff - 150, y = 0.1, label = paste("index =", cutoff), color = "red")

selected.features = names(p_values[(cutoff+1):length(p_values)])
otu.data = otu.data[, selected.features]

rm(scree.plot, selected.features, cutoff, p_values)

# adding response to the last column
otu.data = data.frame(otu.data, Diagnosis = meta.data$Diagnosis)

# shuffling the data set
otu.data = otu.data[sample(nrow(otu.data)), ]

# exporting analysis dataset ----
# this is the cleaned data set that to be used for further analysis. 
library(xlsx)
write.xlsx(otu.data, "Output Datasets/Analysis Dataset.xlsx")

# Ecological Assessment Methods ----

# reading taxonomy data
library(biomformat)
library(qiimer)
library(stringi)

biom.data = read_biom("Original Datasets/BIOM Files/OTU_table.biom")
taxonomy.data = biom_taxonomy(biom.data)
res = as.data.frame(t(stri_list2matrix(taxonomy.data)))
rownames(res) = names(taxonomy.data)
colnames(res) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# TODO: change this to match the name of the previous taxa level concatenated
# with `_Unclassified`
res[is.na(res)] = "Empty Cell"

taxonomy.data = res
rm(res, biom.data)

write.xlsx(taxonomy.data, "Output Datasets/taxonomy.xlsx")

# creating phyloseq object
library(phyloseq)
otu_test = as.matrix(t(data[, -(1:3)]))
tax_test = as.matrix(taxonomy.data)

OTU = otu_table(otu_test, taxa_are_rows = T)
TAX = tax_table(tax_test)

phy.seq = phyloseq(OTU, TAX)
sample_data(phy.seq) = meta.data

phy.seq = prune_taxa(taxa_sums(phy.seq) > 0, phy.seq)

# TODO: (Reem)
# creating LEfSe dataset for Galaxy
create.lefse.dataset = function() {
	# Code here ...
}
lefse.dataset = create.lefse.dataset()

# exporting LEfSe dataset
write.table(lefse.dataset, "/LEfSe_dataset.txt")

# alpha diversity ----
plot = plot_richness(phy.seq, x = "Diagnosis", color = "Sex", measures = c("Chao1", "Shannon","simpson","ace"))
plot = plot + geom_boxplot(data = plot$data, aes(x = Diagnosis, y = value, color = "Sex"), alpha = 0.1)
plot

richness = estimate_richness(phy.seq, split = T, measures = c("Chao1", "Shannon", "simpson", "ace"))
richness = cbind(richness, Diagnosis = sample_data(phy.seq)$Diagnosis)
richness$se.chao1 = NULL
richness$se.ACE = NULL

richness.p.values = sapply(richness[, -ncol(richness)], function(x)
	wilcox.test(x ~ richness$Diagnosis, exact = F)$p.value)
p.adjust(richness.p.values, method = "fdr")

# beta diversity ----
library(vegan)
phy.seq = transform_sample_counts(phy.seq, function(x) x / sum(x) * 100)##
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
	ylab("Relative Abundance (%)") + labs(fill="Abundant Phylum-level OTU")
plot

# comparing abundance values
library(dplyr)
abundances = data.frame(Diagnosis = as.factor(plot$data$Diagnosis), Bacteria = as.factor(plot$data$Phylum), Abundance = plot$data$Abundance)
abundances = abundances %>%
	group_by(Diagnosis, Bacteria) %>%
	summarize(Abundance = sum(Abundance)) %>%
	ungroup %>% as.data.frame

abundances.crc = abundances[abundances$Diagnosis == "CRC", ]
abundances.normal = abundances[abundances$Diagnosis == "Normal", ]
abundances.table = data.frame(Abundances.CRC = abundances.crc$Abundance, Abundances.Normal = abundances.normal$Abundance, row.names = abundances.crc$Bacteria)
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
	ylab("Relative Abundance (%)") + labs(fill="Abundant Genus-level OTU")
plot

abundances = data.frame(Diagnosis = as.factor(plot$data$Diagnosis), Bacteria = as.factor(plot$data$Genus), Abundance = plot$data$Abundance)
abundances = abundances %>%
	group_by(Diagnosis, Bacteria) %>%
	summarize(Abundance = sum(Abundance)) %>%
	ungroup %>% as.data.frame

abundances.crc = abundances[abundances$Diagnosis == "CRC", ]
abundances.normal = abundances[abundances$Diagnosis == "Normal", ]
abundances.table = data.frame(Abundances.CRC = abundances.crc$Abundance, Abundances.Normal = abundances.normal$Abundance, row.names = abundances.crc$Bacteria)
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
abundances.difference = data.frame(Bacteria = rownames(abundances.table), Abundances.Difference = abundances.difference)
abundances.difference = abundances.difference[order(-abs(abundances.difference$Abundances.Difference)), ]
abundances.difference
