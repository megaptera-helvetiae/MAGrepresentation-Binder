# This file contains a few functions including:
#     1.import_gtdbtk_taxonomy_and_checkm 
#     2. import_metadata 
#     3. import_readcounts
#     4. import_tree


#### 1.
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

#### 2.
import_metadata<- function(province_filename, sizeFraction_filename) {
  # Read in the provice info
  prov_dat <- read.table(province_filename, sep = "\t",  as.is=TRUE) 
  
  # Read in the size fraction info
  sizeFrac_dat <- read.table(sizeFraction_filename, sep = "\t",  as.is=TRUE)
  
  # Join together
  metadat <- left_join(prov_dat, sizeFrac_dat, by = "V1") %>%
    distinct() %>%
    column_to_rownames(var = "V1") %>%
    rename(province = "V2.x", size_fraction = "V2.y")
  metadat_df <- sample_data(metadat)
  return(metadat_df)
}


#### 3.
import_readcounts <- function(readcounts_filename){
  otu_physeq <- read.delim(readcounts_filename,  as.is=TRUE) %>%
    dplyr::select(-Length) %>%
    column_to_rownames(var = "X") %>%
    otu_table(taxa_are_rows = TRUE)
  # Return the otu object
  return(otu_physeq)
}    

#### 4.
import_tree <- function(tree_filename){
  tree <- read_tree(tree_filename)
  taxa_names(tree) <- gsub(taxa_names(tree), pattern = "-", replacement = "_")
  # Return the tree object
  return(tree)
}


