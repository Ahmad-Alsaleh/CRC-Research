create.lefse.dataset <- function(physeq_pruned) {
  physeq_pruned <- transform_sample_counts(physeq_pruned, function(x) x / sum(x))

  taxonomy_data <- as.data.frame(tax_table(physeq_pruned))
  colnames(taxonomy_data) <- paste0("taxonomy", seq(1:7))

  ps_kingdom <- tax_glom(physeq_pruned, "Kingdom")
  taxa_names(ps_kingdom)
  taxa_names(ps_kingdom) <- tax_table(ps_kingdom)[, 1]
  taxa_names(ps_kingdom)
  otu_kingdom <- as.data.frame(otu_table(ps_kingdom))

  ps_phylum <- tax_glom(physeq_pruned, "Phylum")
  taxa_names(ps_phylum)
  taxa_names(ps_phylum) <- tax_table(ps_phylum)[, 2]
  taxa_names(ps_phylum)
  otu_phylum <- as.data.frame(otu_table(ps_phylum))

  match_phylum <- taxonomy_data$taxonomy1[match(rownames(otu_phylum), taxonomy_data$taxonomy2)]
  new_names_otu_phylum <- paste0(match_phylum, "|", rownames(otu_phylum))
  rownames(otu_phylum) <- new_names_otu_phylum

  ps_class <- tax_glom(physeq_pruned, "Class")
  taxa_names(ps_class)
  taxa_names(ps_class) <- tax_table(ps_class)[, 3]
  taxa_names(ps_class)
  otu_class <- as.data.frame(otu_table(ps_class))

  match_class <- taxonomy_data$taxonomy2[match(rownames(otu_class), taxonomy_data$taxonomy3)]
  match_phylum <- taxonomy_data$taxonomy1[match(match_class, taxonomy_data$taxonomy2)]
  new_names_otu_class <- paste0(
    match_phylum, "|",
    match_class, "|", rownames(otu_class)
  )
  rownames(otu_class) <- new_names_otu_class

  ps_order <- tax_glom(physeq_pruned, "Order")
  taxa_names(ps_order)
  taxa_names(ps_order) <- tax_table(ps_order)[, 4]
  taxa_names(ps_order)
  otu_order <- as.data.frame(otu_table(ps_order))

  match_order <- taxonomy_data$taxonomy3[match(rownames(otu_order), taxonomy_data$taxonomy4)]
  match_class <- taxonomy_data$taxonomy2[match(match_order, taxonomy_data$taxonomy3)]
  match_phylum <- taxonomy_data$taxonomy1[match(match_class, taxonomy_data$taxonomy2)]
  new_names_otu_order <- paste0(
    match_phylum, "|", match_class, "|", match_order, "|",
    rownames(otu_order)
  )

  rownames(otu_order) <- new_names_otu_order

  ps_family <- tax_glom(physeq_pruned, "Family")
  taxa_names(ps_family)
  taxa_names(ps_family) <- tax_table(ps_family)[, 5]
  taxa_names(ps_family)
  otu_family <- as.data.frame(otu_table(ps_family))

  match_family <- taxonomy_data$taxonomy4[match(rownames(otu_family), taxonomy_data$taxonomy5)]
  match_order <- taxonomy_data$taxonomy3[match(match_family, taxonomy_data$taxonomy4)]
  match_class <- taxonomy_data$taxonomy2[match(match_order, taxonomy_data$taxonomy3)]
  match_phylum <- taxonomy_data$taxonomy1[match(match_class, taxonomy_data$taxonomy2)]
  new_names_otu_family <- paste0(
    match_phylum, "|", match_class, "|", match_order,
    "|", match_family, "|", rownames(otu_family)
  )

  rownames(otu_family) <- new_names_otu_family

  ps_genus <- tax_glom(physeq_pruned, "Genus")
  taxa_names(ps_genus)
  taxa_names(ps_genus) <- tax_table(ps_genus)[, 6]
  taxa_names(ps_genus)
  otu_genus <- as.data.frame(otu_table(ps_genus))

  match_genus <- taxonomy_data$taxonomy5[match(rownames(otu_genus), taxonomy_data$taxonomy6)]
  match_family <- taxonomy_data$taxonomy4[match(match_genus, taxonomy_data$taxonomy5)]
  match_order <- taxonomy_data$taxonomy3[match(match_family, taxonomy_data$taxonomy4)]
  match_class <- taxonomy_data$taxonomy2[match(match_order, taxonomy_data$taxonomy3)]
  match_phylum <- taxonomy_data$taxonomy1[match(match_class, taxonomy_data$taxonomy2)]
  new_names_otu_genus <- paste0(
    match_phylum, "|", match_class, "|", match_order,
    "|", match_family, "|", match_genus, "|", rownames(otu_genus)
  )
  rownames(otu_genus) <- new_names_otu_genus

  Diagnosis <- as.character(sample_data(physeq_pruned)$Diagnosis)
  galaxy <- rbind(
    Diagnosis = Diagnosis, otu_kingdom, otu_phylum, otu_class, otu_order,
    otu_family, otu_genus, make.row.names = T
  )

  return(galaxy)
}

# replaces NA with [PREVIOUS]_Unclassified
handle.unclassified.values <- function(tax) {
  for (i in 1:7) {
    tax[, i] <- as.character(tax[, i])
  }

  tax[is.na(tax)] <- ""
  for (i in 1:nrow(tax)) {
    if (tax[i, 2] == "") {
      kingdom <- paste0(tax[i, 1], "_Unclassified")
      tax[i, 2:7] <- kingdom
    } else if (tax[i, 3] == "") {
      phylum <- paste0(tax[i, 2], "_Unclassified")
      tax[i, 3:7] <- phylum
    } else if (tax[i, 4] == "") {
      class <- paste0(tax[i, 3], "_Unclassified")
      tax[i, 4:7] <- class
    } else if (tax[i, 5] == "") {
      order <- paste0(tax[i, 4], "_Unclassified")
      tax[i, 5:7] <- order
    } else if (tax[i, 6] == "") {
      family <- paste0(tax[i, 5], "_Unclassified")
      tax[i, 6:7] <- family
    } else if (tax[i, 7] == "") {
      tax$Species[i] <-
        paste0(tax$Genus[i], "_Unclassified")
    }
  }

  return(tax)
}

# adds prefix to each taxa level
add.prefix <- function(tax) {
  is.na(tax) <- tax == "Unclassified"
  tax$Kingdom <- str_c("K_", tax$Kingdom)
  tax$Phylum <- str_c("P_", tax$Phylum)
  tax$Class <- str_c("C_", tax$Class)
  tax$Order <- str_c("O_", tax$Order)
  tax$Family <- str_c("F_", tax$Family)
  tax$Genus <- str_c("G_", tax$Genus)
  tax$Species <- str_c("S_", tax$Species)

  return(tax)
}
