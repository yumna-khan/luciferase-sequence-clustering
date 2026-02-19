# Project Overview
This project explores whether alignment-free sequence representations (amino acid k-mers) can recover biologically meaningful structure in luciferase and photoprotein sequences.

Bioluminescent proteins are widely used in biomedical imaging, drug discovery, in vivo imaging, and protein–protein interaction studies. Because bioluminescence has evolved independently multiple times across the tree of life, luciferases are highly diverse in structure and evolutionary origin. This raises an important methodological question: can simple sequence-composition features capture functional or evolutionary signals?

Using unsupervised learning approaches, this study evaluates whether k-mer–based features cluster luciferase sequences more strongly by:
1. Functional type (luciferase classification)
or
2. Taxonomic group (evolutionary relationships)

This helps assess whether k-mer–based approaches capture dominant biological signals in luciferase evolution.

# Methodology
## 1. Dataset Description
The dataset was obtained from the National Library of Medicine and consists of 200 luciferase and photoprotein sequences, categorized into XII luciferase types and unknown classifications. For this study, the primary variables used were Luciferase type, Taxonomy, and NCBI accession numbers (used to retrieve protein sequences)

## 2. Data Exploration
- Inspected dataset structure
- Visualized luciferase type and taxon distributions

## 3. Preprocessing
- Filtered relevant columns
- Split multiple accession entries
- Corrected formatting inconsistencies
- Retrieved protein sequences from NCBI

## 4. Feature Engineering
- Constructed amino acid k-mer features (k = 2)
- Normalized/scaled feature matrix
- Computed Manhattan distance matrix

## 5. Unsupervised Learning & Evaluation
- UMAP dimensionality reduction
- Hierarchical clustering (dendrograms)
- K-means clustering
- Silhouette analysis

# Key Findings
Across all clustering approaches (UMAP, hierarchical clustering, and k-means), both luciferase type and taxonomy showed partial clustering patterns, but neither formed clearly separated groups. Some luciferase types (e.g., Type V) and some taxa (e.g., Insecta) displayed localized clustering, yet many clusters were heterogeneous.

Silhouette scores were low (0.28 for luciferase type and 0.37 for taxonomy), indicating weak separation and substantial overlap between groups. Although taxonomy performed slightly better, k-mer composition alone did not strongly distinguish functional or evolutionary categories.

These findings align with biological expectations: bioluminescence is a polyphyletic trait that evolved independently multiple times. Functional convergence and evolutionary divergence weaken clustering signals when using simple sequence composition features.

Overall, alignment-free k-mer methods were insufficient to produce strongly biologically meaningful clustering of luciferase sequences.

# Limitations
Several factors may have influenced the results of this analysis. Sampling was uneven across luciferase types and taxa, which may have biased cluster formation. Additionally, k-mer features capture only overall amino acid composition and do not account for structural domains or functionally important residues, limiting biological interpretability. There may also be inconsistencies in luciferase annotations, as labeling is not always standardized across datasets. Furthermore, unsupervised methods such as UMAP and k-means are sensitive to parameter choices, including dimensionality settings and cluster number selection, which can affect clustering outcomes. The overall weak separation observed in the clusters suggests that the detected patterns may reflect limitations of the feature representation rather than true underlying biological structure.

# Future Directions
Future work could improve upon this analysis in several ways. Applying supervised learning approaches would allow direct evaluation of how well sequence features predict luciferase type or taxonomy. Exploring alternative clustering methods, such as DBSCAN, spectral clustering, or semi-supervised techniques, may help detect more subtle patterns. Incorporating larger and more taxonomically diverse datasets could improve robustness and generalizability. Additionally, using higher-order k-mers or incorporating structural and functional domain features may better capture deeper evolutionary relationships and functional constraints within luciferase sequences.
