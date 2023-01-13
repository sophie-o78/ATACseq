# Proteomics
The goal is to study protein mutations in various samples.

It is easier to use the IDE RStudio to display the graphs.

To run the code, several libraries/ packages must be installed: 
    
1. The dendextend package
2. The pheatmap package
3. The seriation package

They can be installed through the following command :
install.packages("name_package")

4. The limma package must be installed through BiocManager:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")

