# R 4.4.2, Bioconductor 3.20 (dada2 1.34)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages("vegan")
BiocManager::install(c("dada2", "Biostrings", "ShortRead"), version = "3.20")
