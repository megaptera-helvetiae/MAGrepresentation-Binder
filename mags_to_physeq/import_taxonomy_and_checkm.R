

#This function reads in a
taxonomy_filename <- "Tara_Oceans_Med/TOBG-MED-READCOUNTMATCH.bac120.tsv"
checkm_filename <- "Tara_Oceans_Med/TOBG-MED_qa.txt"

import_gtdbtk_taxonomy_and_checkm <- function(taxonomy_filename, checkm_filename) {
  
  library(phyloseq)
  library(tidyverse)
  
tax_raw <- read.csv(taxonomy_filename, sep = "\t", header = FALSE, as.is=TRUE)
split_first_col <- matrix(unlist(strsplit(tax_raw[,1], split = ".", fixed = TRUE)), ncol=5, byrow = TRUE)

magnames <- split_first_col[,1]
taxonomy <- matrix(unlist(strsplit(tax_raw[,2], split = ";", fixed = TRUE)), ncol=7, byrow = TRUE)

rownames(taxonomy) <-magnames

checkm_dat <- read.csv(checkm_filename, sep = "\t", as.is = TRUE)
rownames(checkm_dat) <- checkm_dat$Bin.Id

common_samples <- intersect(rownames(checkm_dat),rownames(taxonomy))
checkm_subsetted <- checkm_dat[common_samples,]
taxonomy_subsetted <- taxonomy[common_samples,]

combined <- cbind(taxonomy_subsetted, checkm_subsetted)

colnames(combined)[1:7] <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
combined_mat <- as.matrix(combined)
combined_tbl <- as_tibble(combined_mat)
sub_prefixes <- function(x){gsub(x, pattern = "[a-z]__", replacement = "")}
combined_tbl_no_prefixes <- dplyr::mutate_at(.tbl = combined_tbl, .vars =colnames(combined)[1:7], .funs = sub_prefixes )

#label_unassigned <- function(x){gsub(x, pattern = "", replacement = "Unassigned")}
#combined_tbl_no_prefixes <- dplyr::mutate_at(.tbl = combined_tbl_no_prefixes, .vars =colnames(combined)[1:7], .funs = label_unassigned )


combined_matrix_2 <- as.matrix(combined_tbl_no_prefixes,ncol = ncol(combined_tbl_no_prefixes))
rownames(combined_matrix_2) <- rownames(combined_mat)
taxtab <- tax_table(combined_matrix_2)

return(taxtab)
}

# Import the tax data
tax <- import_gtdbtk_taxonomy_and_checkm(taxonomy_filename = "Tara_Oceans_Med/TOBG-MED-READCOUNTMATCH.bac120.tsv", checkm_filename = "TOBG-MED_qa.txt")

# Import Readcounts
otu <- read.delim("Tara_Oceans_Med/TOBG-MED-TOTAL.readcounts") %>%
  dplyr::select(-Length) %>%
  column_to_rownames(var = "X") %>%
  otu_table(taxa_are_rows = TRUE)



tara_med <- phyloseq(tax, otu)


