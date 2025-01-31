# Tastemel-Jussila_2025
# This repository is supporting the manuscript titled "Context-Dependent and Gene-Specific Role of Chromatin Architecture Mediated by Histone Modifiers and Loop-extrusion Machinery," Melodi Tastemel†, Adam Jussila†, Bharath Saravanan, et. al. (2025)

Files contained here:

## Common Tools (directory)
A collection of tools and supporting files created by Bogdan Bintu during his time in the Xiaowei Zhuang lab and modified by Adam Jussila for the purposes of analyzing and understanding the chromatin tracing experiments performed in this study.

## sgRNA_scripts (directory)
Contains preparse.sh and parseSeq1_from_fastq.pl, which are used to take fastq files and extract the reads mapping to sgRNAs from the sgRNA table brie.lib.csv. All reads not containing sgRNA identifiers are discarded. Processed files from the output of this script are used to run MAGeCK for downstream analysis.

## chromatic_corr_BB_5_28_2020.pkl
A chromatic correction file used to correct for our computed chromatic effects for our specific microscope system used for these experiments.

## ChromTracer_Sox2.ipynb
The primary notebook used to analyze the chromatin tracing data from raw data through to 3D traces. Includes some additional analysis at the end which includes E-P distances, etc.

## clustering_embedding.ipynb
Notebook used for clustering using SnapATAC2. This uses some outputs from the ChromTracer_Sox2.ipynb, and generated UMAP embeddings seen in the manuscript.

## geometric_center_analysis.ipynb
Notebook containing centrality analyses as well as some correlation analyses, etc. This is for additional panels seen in the manuscript in many of the figures.

