#***************************
# MBINF 6210
#
# Yumna Khan - 1094825
#
# 2025-11-18
#***************************


#Goal
# To determine whether sequence similarity patterns (based on protein amino acid composition using k-mer features) naturally separate luciferases into their known biochemical types, or whether they cluster more strongly by taxonomy.



# ----- SECTION 0: Dependencies -------

# Uncomment to install if needed
#install.packages("kmer")
#install.packages("pheatmap")
#install.packages("umap")
#install.packages("Rtsne")
#install.packages("patchwork")
#install.packages("mclust")


# Libraries
library(tidyverse)
library(viridis)
library(readxl)
library(rentrez)
library(seqinr)
library(cluster)
library(kmer)
library(stats)
library(pheatmap)
library(umap)
library(patchwork)
library(mclust)

# Set preferred conflicts here:
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::slice)
conflicted::conflicts_prefer(dplyr::rename)


# *********************************************************************
# -------------------- SECTION 1: Data Input --------------------------
# *********************************************************************

#Read Excel file
raw_data = read_excel("../data/Table_1_Leaving the Dark Side_ Insights Into the Evolution of Luciferases.XLSX", sheet = 1)

# Show structure of file
class(raw_data)
summary(raw_data)
dim(raw_data)
names(raw_data)
lapply(raw_data, head)


# *********************************************************************
# ---------------- SECTION 2: Initial Data Exploration ----------------
# *********************************************************************

# View key variables
table(raw_data$`Luciferase/Photoprotein-type`)
table(raw_data$`NCBI Accession`)
table(raw_data$`Taxon`)

# View NA's
sum(is.na(raw_data$Taxon))
sum(is.na(raw_data$`Luciferase/Photoprotein-type`))

# View unique Taxon, Luciferase, and Accessions
(unique(raw_data$`Taxon`))
(unique(raw_data$`Luciferase/Photoprotein-type`))
(unique(raw_data$`NCBI Accession`))

# Notice there were no NA's. The data doesn't contain NA's or empty cells specifically, but rather an arbitrary assignment like "Unknown" or "/".

# Taxon was used instead of species because species labels were too sparse and many appeared only once, while phylum labels would be too general. Thus, taxon groups had enough replication to support meaningful comparisons


# To view any biases within the data, plot each category (luciferse type and taxonomy)

# Create new dataframe of only luciferase to sort counts
top_luc <- raw_data %>%
  count(`Luciferase/Photoprotein-type`, sort = TRUE)

# Plot Luciferase type - provide x and y data for sorting
ggplot(top_luc, aes(x = reorder(`Luciferase/Photoprotein-type`, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Luciferase Types",
    y = "Count",
    title = "Luciferase Types Vs Count"
  )

# Create new dataframe of only taxon to sort counts
top_tax <- raw_data %>%
  count(Taxon, sort = TRUE)

# Plot Taxonomy
ggplot(top_tax, aes(x = reorder(Taxon, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Taxonomy Names",
    y = "Count",
    title = "Taxonomy vs Count"
  )


# *********************************************************************
# ----------------- SECTION 3: Cleaning/Preprocessing -----------------
# *********************************************************************


### ------------- Preliminary Filtering & Splitting -------------

# Create new dataframe to discard unused columns, rename, filter, and trim excessive spaces
luc_data = raw_data %>% 
  select("Luciferase/Photoprotein-type", "Taxon", "NCBI Accession") %>%
  rename(
    Luciferase_type= "Luciferase/Photoprotein-type", 
    Accession = `NCBI Accession`) %>% 
  filter(!(Luciferase_type %in% c("Unknown", NA)))

luc_data

# Split accessions accordingly (i.e. if there are 2 accessions in a cell separated by "/" or ",", make it another row) and filters empty accession columns
luc_data_split <- luc_data %>%
  mutate(Accession_list = strsplit(Accession, "[/,]")) %>% 
  mutate(Accession_list = lapply(Accession_list, function(x) str_trim(x))) %>% 
  unnest(cols = c(Accession_list)) %>%
  select(-Accession) %>% 
  rename(Accession = Accession_list) %>% 
  filter(Accession != "")
luc_data_split  


# Sanity check: verify row counts, missing/empty values, and duplicates in accessions
list(
  total_rows = nrow(luc_data_split),
  empty_accessions = sum(luc_data_split$Accession == ""),
  missing_accessions = sum(is.na(luc_data_split$Accession)),
  unique_accessions = length(unique(luc_data_split$Accession)),
  duplicate_accessions = luc_data_split %>% 
    count(Accession) %>% 
    filter(n > 1)
)

# Sanity check: assign summary statistics for later reference
total_rows_before <- nrow(luc_data_split)
unique_accessions_before <- length(unique(luc_data_split$Accession))



### ---------------- Sequence Extraction ----------------

# Manual changes to NCBI Accessions
# These manual corrections were applied because the parsed accessions contained broken IDs in the original dataset when attempting to use entrez_fetch

luc_data_split$Accession[25] = "O77206"
luc_data_split$Accession[56] = "1SL9_A"
luc_data_split$Accession[57] = "1QV1_A" 
luc_data_split$Accession[67] = "3KPX_A"

luc_data_split$Accession[92] = "AAB26932"
luc_data_split$Accession[93] = "AAN73267"
luc_data_split$Accession[21] = "AAV35380"
luc_data_split$Accession[95] = "ABA03040"
luc_data_split$Accession[64] = "AEP19818"
luc_data_split$Accession[96] = "BAU71688"
luc_data_split$Accession[94] = "BAW94790"
luc_data_split$Accession[145] = "JAV81274"
luc_data_split$Accession[49] = "P02592"
luc_data_split$Accession[54] = "4MRY_A"
luc_data_split$Accession[63] = "4NQG_A"
luc_data_split$Accession[62] = "4NQG_B"


# Sanity check to confirm the manual corrections were applied correctly and do not introduce duplicates
unique(luc_data_split$Accession)


# Classify the accession into 3 categories (the most common NCBI databases were protein and nucleotides, remaining are unknown)
classify_accession <- function(acc) {
  
  # Try protein database
  is_protein <- tryCatch(
    entrez_search(db = "protein", term = acc)$count, # Search database and count how many hits were found
    error = function(e) 0 # If search fails, return 0
  )
  if (is_protein > 0) # If hit found, return protein
    return("protein")
  
  # Try nuccore
  is_nuccore <- tryCatch(
    entrez_search(db = "nuccore", term = acc)$count,
    error = function(e) 0
  )
  if (is_nuccore > 0) 
    return("nuccore")
  
  # If neither works
  return("unknown")
}

# Create new column for accession type 
luc_data_split$Accession_type <- sapply(luc_data_split$Accession, classify_accession)

# Remove the rows that aren't protein
luc_data_clean = luc_data_split %>% 
  filter(Accession_type == "protein")

# View counts for each accession type
table(luc_data_split$Accession_type)



# Function to retrieve protein sequence using NCBI accession id
fn_sequence = function(accession_id){
  
  # Fetch raw FASTA format from the protein database
  fasta_raw <- entrez_fetch(db = "protein", id = accession_id, rettype = "fasta")
  
  # Check if the fetch returned a valid FASTA (must start with '>')
  if (!startsWith(fasta_raw, ">")) {
    return(NA) # Return NA if sequence not found or invalid
  }
  
  # Split the FASTA into lines
  lines <- unlist(strsplit(fasta_raw, "\n")) 
  
  # Remove header (starts with '>')
  seq_lines <- lines[!grepl("^>", lines)] 
  
  # Collapse into long sequence
  sequence <- paste(seq_lines, collapse = "")
  
  # Remove unwanted characters from sequence
  sequence_clean = str_remove_all(sequence,"[XBZ\\-*]")
  
  # Return the cleaned protein sequence
  return(sequence_clean)
}


# In a new column, add sequence
luc_data_clean$Sequence <- sapply(luc_data_clean$Accession, fn_sequence)

# Compute the sequence length
luc_data_clean$Length = str_count(luc_data_clean$Sequence)


# Sanity check: view length of sequences, counts per luciferase type and taxon, count if unwanted characters still present
summary(luc_data_clean$Length)
table(luc_data_clean$Luciferase_type)
table(luc_data_clean$Taxon)
table(str_count(luc_data_clean$Sequence, "[XBZ\\-*]"))

# Sanity check: assign summary statistics
total_rows_after <- nrow(luc_data_clean)
unique_accessions_after <- length(unique(luc_data_clean$Accession))

# Sanity check: create a function to compare before and after unique accession change, NAs, and total rows to confirm no data was lost
sanity_summary <- function(df, name="Dataset") {
  cat(paste0("--- ", name, " ---\n"))
  cat("Total rows: ", nrow(df), "\n")
  cat("Unique accessions: ", length(unique(df$Accession)), "\n")
  cat("Sequences with NA: ", sum(is.na(df$Sequence)), "\n\n")
}

# Sanitiy check: apply function 
sanity_summary(luc_data_split, "Before Filtering") # Will give a warning message for unknown 'Sequence' since it was not created by then
sanity_summary(luc_data_clean, "After Filtering")



### ---------------- Visualization of Cleaned Data ----------------

# Histogram of sequence lengths to detect unusually short or long sequences
# Shows the distribution of protein lengths- 50 bins  chosen to balance detail and readability
ggplot(luc_data_clean, aes(x = Length)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(
    x = "Sequence Length",
    y = "Number of Proteins"
  ) +
  labs(title = "Histogram of Protein Sequence Lengths")


# Bar plot of taxonomic composition to check for over or under represented taxons
# Counts sequences per taxon to assess dataset after sequence filtering
ggplot(luc_data_clean, aes(x = Taxon)) +
  geom_bar() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Taxon",
    y = "Count"
  ) +
  labs(title="Taxonomic Composition")


# Note: Removing outliers could remove meaningful biological variation, especially since rare or unusually long sequences may still contain important functional or evolutionary information. Thus, no outliers were removed in this study.




# *********************************************************************
# ------------------ SECTION 4: Feature Engineering -------------------
# *********************************************************************

### ---------------- Amino Acid Heatmap ----------------

# Preprocessing for Heatmap - amino acid composition patterns across luciferase types/taxonomy can reveal if certain amino acids are enriched or depleted in functional classes, providing exploratory insight before clustering

# Define standard amino acids
aa_letters <- c("A","C","D","E","F","G","H","I","K","L",
                "M","N","P","Q","R","S","T","V","W","Y")

# Compute amino acid counts for each sequence
aa_counts <- lapply(luc_data_clean$Sequence, function(seq) {
  seq_vec <- strsplit(seq, "")[[1]] # Split sequence into individual letters
  seq_vec <- seq_vec[seq_vec %in% aa_letters]  # Drop non standard AAs
  counts <- table(seq_vec) # Frequency of all AAs
  counts_full <- setNames(rep(0, length(aa_letters)), aa_letters) # Create a vector of 0s for all amino acids so all sequences have the same vector length
  counts_full[names(counts)] <- counts # Fill in the actual counts
  return(counts_full)
})

# Convert to a data frame and add grouping variables
aa_df <- as.data.frame(do.call(rbind, aa_counts))
aa_df$Luciferase_type <- luc_data_clean$Luciferase_type
aa_df$Taxon <- luc_data_clean$Taxon

# Function to compute mean AA composition per group and convert to matrix
aa_matrix <- function(df, group_col) {
  aa_mean <- df %>%
    group_by(across(all_of(group_col))) %>%
    summarise(across(all_of(aa_letters), mean), .groups = "drop") # get mean of all AA counts per group
  mat <- as.matrix(aa_mean[, aa_letters])
  rownames(mat) <- aa_mean[[group_col]]
  mat # return matrix
}

aa_matrix_luc <- aa_matrix(aa_df, "Luciferase_type")
aa_matrix_tax <- aa_matrix(aa_df, "Taxon")

# Plot heatmap of luciferase
pheatmap(aa_matrix_luc,
         main = "Average Amino Acid Composition by Luciferase Type",
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE)

# Plot heatmap of luciferase
pheatmap(aa_matrix_tax,
         main = "Average Amino Acid Composition by Taxonomy",
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE)



### ---------------- Kmer Construction ----------------

# Convert sequences into character vectors of single letters
seq_list <- lapply(luc_data_clean$Sequence, function(s) strsplit(s, "")[[1]])

# Function for kmer length of each sequence
get_kmers <- function(seq, k) {
  n <- nchar(seq) # Get seq length
  if (n < k) return(NULL) # return null if seq < k
  sapply(1:(n - k + 1), function(i) substr(seq, i, i + k - 1)) # Slide a window of length k along the sequence and extract each substring (kmer)
}

# Build k-mer table for k=2
# Kmer was changed between 1 - 5 before kmer = 2 was chosen as the optimal k. The mean silhouette width was used as a reference (closer to 1 = better clustering). 
k <- 2
kmer_list <- lapply(luc_data_clean$Sequence, get_kmers, k = k)

# Combine into a matrix
all_kmers <- unique(unlist(kmer_list))
kmer_matrix <- t(sapply(kmer_list, function(kv) {
  tab <- table(kv) # Count kmers
  s <- rep(0, length(all_kmers)) # Initialize vector of 0s
  names(s) <- all_kmers
  if (!is.null(tab)) s[names(tab)] <- tab # Fill counts for present kmers
  s # Return s
}))

# Check matrix dimensions and view first 10 kmer columns
dim(kmer_matrix)
head(kmer_matrix[,1:10])



### ---------------- Normalization/Scaling ----------------

# Identify kmers that have zero variance across all sequences (same count everywhere)
# argument 2 means apply the function over the columns (kmers)
zero_var_cols <- apply(kmer_matrix, 2, function(x) var(x) == 0)
sum(zero_var_cols)

# Remove zero variance kmers; they don’t provide any information for clustering
kmer_matrix_clean <- kmer_matrix[, !zero_var_cols]

# Scale features for clustering/analysis
kmer_scaled <- scale(kmer_matrix_clean)
head(kmer_scaled[, 1:10])

# Check if scaling was done correctly
colMeans(kmer_scaled)  # should be ~0
apply(kmer_scaled, 2, sd)  # should be ~1
apply(kmer_scaled, 2, var)  # should be 1

# Quick plot if there are any outliers from scaling
boxplot(kmer_scaled[, 1:10], 
        main="First 10 Scaled K-mers",
        xlab = "K-mer",
        ylab = "Scaled Value"
        )


# Compute pairwise distance matrix using Manhattan distance
# Manhattan distance was chosen as it emphasizes larger differences without squaring them like Euclidean, giving better clustering results
dist_matrix <- dist(kmer_scaled, method="manhattan")



# *********************************************************************
# ------------------ SECTION 5: Clustering -----------------
# *********************************************************************


### ---------------- PCA ----------------

# Perform PCA on the k-mer distance matrix to reduce high dimensional features and overall cluster patterns

# PCA
# center = TRUE - subtracts the mean of each feature before PCA
# scale. = TRUE - scales features to unit variance so all kmers contribute equally
pca_results <- prcomp(dist_matrix, center = TRUE, scale. = TRUE)

# Create a data frame for plotting
pca_df <- data.frame(
  PC1 = pca_results$x[,1],
  PC2 = pca_results$x[,2],
  Luciferase_type = luc_data_clean$Luciferase_type,
  Taxon = luc_data_clean$Taxon
)

# Calculate the % variance explained by the first two PCs. Tells us how much of the total variation in the data is captured by these components
var_explained <- summary(pca_results)$importance[2,1:2] * 100

# PCA coloured by Luciferase type
ggplot(pca_df, aes(x = PC1, y = PC2, color = Luciferase_type)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, option = "D") +
  theme_bw() +
  labs(
    title = "PCA of Kmer Features Coloured by Luciferase Type",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  )

# PCA coloured by Taxon
ggplot(pca_df, aes(x = PC1, y = PC2, color = Taxon)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, option = "C") +
  theme_bw() +
  labs(title = "PCA of Kmer Features Coloured by Taxon",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  )


### ---------------- UMAP ----------------

# Perform UMAP on the distance matrix
# UMAP is a nonlinear dimensionality reduction method that preserves local and global structure better than PCA, especially for clustering and complex high dimensional data like kmer frequencies

# Create a data frame with UMAP coordinates and metadata for plotting
u <- umap(as.matrix(dist_matrix))
umap_df <- data.frame(
  UMAP1 = u$layout[,1],
  UMAP2 = u$layout[,2],
  Luciferase_type = luc_data_clean$Luciferase_type,
  Taxon = luc_data_clean$Taxon
)

# UMAP coloured by Luciferase type
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Luciferase_type)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, option = "H") +
  theme_bw() +
  labs(title = "UMAP of Kmer Features Coloured by Luciferase Type",
       color = "Luciferase Type")

# UMAP coloured by Taxon
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Taxon)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, option = "H") +
  theme_bw() +
  labs(title = "UMAP of Kmer Features Coloured by Taxon",
       color = "Taxonomy")



### ---------------- Hierarchical Clustering ----------------

# Hierarchical clustering
# Used ARI index to determine which method works best (tried different methods such as complete, ward.D2, etc) - (closer to 1 = better clustering)
# Used Ward.D2 linkage because it minimizes within cluster variance and produces compact, interpretable clusters
hc <- hclust(dist_matrix, method = "ward.D2")

# Plot dendrogram for Luciferase
plot(hc, labels = luc_data_clean$Luciferase_type,
     cex = 0.5,    # shrink text
     las = 2,      # rotate label text vertically
     main = "Dendrogram by Luciferase Type",
     xlab = "Luciferase Type", sub = "")

# Plot dendrogram for Taxon
plot(hc, labels = luc_data_clean$Taxon,
     cex = 0.5,
     las = 2,
     main = "Dendrogram by Taxonomy",
     xlab = "Taxonomy", sub = "")



### ---------------- K means Clustering ----------------

set.seed(123) # Ensures k-means clustering gives reproducible results

# Unique length of each category
kL = length(unique(luc_data_clean$Luciferase_type))
kT = length(unique(luc_data_clean$Taxon))

# k means on kmer features. Random initialization - seed required for reproducibility
kmLuc <- kmeans(kmer_scaled, centers = kL)
kmTax <- kmeans(kmer_scaled, centers = kT)

# View how many sequences in each cluster belong to each true class
km_luc_tab <- table(kmLuc$cluster, luc_data_clean$Luciferase_type)
km_tax_tab <- table(kmTax$cluster, luc_data_clean$Taxon)

# Plot heatmap to show the overall clustering structure and reveal broad patterns of similarity between clusters and biological labels

# Plot luciferase
pheatmap(
  km_luc_tab,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "K-means vs Luciferase Type",
)

# Plot taxon
pheatmap(
  km_tax_tab,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "K-means vs Taxon",
)

# Plot barplots to clearly quantify how many sequences of each label fall into each cluster, making cluster composition easy to interpret

# Preprocessing for Barplot to easily view clustering composition of both Luciferase and Taxon
df_luc <- as.data.frame(km_luc_tab) %>%
  rename(Cluster = Var1, Label = Var2, Count = Freq) %>%
  mutate(Type = "Luciferase")

df_tax <- as.data.frame(km_tax_tab) %>%
  rename(Cluster = Var1, Label = Var2, Count = Freq) %>%
  mutate(Type = "Taxon")

#Plot for Luciferase
p_luc <- ggplot(df_luc, aes(x = Cluster, y = Count, fill = Label)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  labs(title = "Luciferase Cluster Composition",
       x = "Cluster",
       y = "Number of Sequences",
       fill = "Luciferase Type") +
  scale_fill_viridis_d(option = "plasma") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right",
        legend.box = "vertical") +
  guides(fill = guide_legend(nrow = 6, byrow = FALSE))

# Plot for Taxon
p_tax <- ggplot(df_tax, aes(x = Cluster, y = Count, fill = Label)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  labs(title = "Taxon Cluster Composition",
       x = "Cluster",
       y = "Number of Sequences",
       fill = "Taxonomy") +
  scale_fill_viridis_d(option = "H") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right",
        legend.box = "vertical") +
  guides(fill = guide_legend(nrow = 6, byrow = FALSE))

# Combine the plots vertically
p_combined <- p_luc / p_tax

# Print the final combined figure
print(p_combined)




### ---------------- Silhouette Plot ----------------
# Silhouette width is computed to evaluate how well each sequence fits within its assigned cluster and how distinct the clusters are from one another

# Compute silhouette information
sil_luc <- silhouette(kmLuc$cluster, dist_matrix)
sil_tax <- silhouette(kmTax$cluster, dist_matrix)

# Average silhouette width
mean(sil_luc[, 3])   
mean(sil_tax[, 3])

# Preprocessing for Silhouette plot
# Function to tidy silhouette info
tidy_sil <- function(sil_obj, cluster_name) {
  sil_df <- data.frame(
    cluster = factor(sil_obj[, "cluster"]), # cluster for each seq
    sil_width = sil_obj[, "sil_width"], # sil width for each seq
    obs = 1:nrow(sil_obj) # index used for plotting in order
  )
  sil_df$group <- cluster_name 
  return(sil_df)
}

sil_luc_df <- tidy_sil(sil_luc, "Luciferase")
sil_tax_df <- tidy_sil(sil_tax, "Taxon")

# Combine for plotting
sil_combined <- bind_rows(sil_luc_df, sil_tax_df)

# Plot
ggplot(sil_combined, aes(x = obs, y = sil_width, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~group, ncol = 1, scales = "free_x") +
  scale_fill_viridis_d(option = "D") +
  theme_bw() +
  labs(
    x = "Sequence",
    y = "Silhouette Width",
    fill = "Cluster",
    title = "Silhouette Plot of K-means Clustering"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0)
  )



### ---------------- ARI Index ----------------

# Compute Adjusted Rand Index (ARI) to quantify how well the clustering matches known labels
# ARI ranges from -1 to 1, where 1 = perfect agreement, 0 = random agreement, negative = worse than random


# Perform hclust using a specific method (e.g., Ward.D2)
change_method = "ward.D2"
hclust_test <- hclust(dist_matrix, method = change_method)

# Number of unique  groups
KL <- length(unique(luc_data_clean$Luciferase_type))
KT <- length(unique(luc_data_clean$Taxon))

# Cut the hierarchical tree into K groups
test_clustersL <- cutree(hclust_test, k = KL)
test_clustersT <- cutree(hclust_test, k = KT)

# Compute ARI: compares clustering output to true labels
ari_test_luc <- adjustedRandIndex(test_clustersL, luc_data_clean$Luciferase_type)
ari_test_tax <- adjustedRandIndex(test_clustersT, luc_data_clean$Taxon)

ari_test_luc
ari_test_tax

