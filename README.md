📌 Overview

This repository contains the complete analysis pipeline and scripts used in our study on Paroxysmal Nocturnal Haemoglobinuria (PNH), integrating bulk mRNA-seq, miRNA-seq, and single-cell RNA-seq (scRNA-seq) data to identify key regulatory axes and potential therapeutic agents.

🧬 Background

PNH is a clonal haematopoietic stem cell disorder characterized by:

1. Intravascular hemolysis
2. Thrombosis
3. Bone marrow failure
   
Although complement inhibitors are widely used in clinical practice, their efficacy is limited, underscoring the need for novel therapeutic targets and alternative treatment strategies. This study aims to provide a multi-omics–based theoretical foundation for PNH pathogenesis and drug repurposing.

🔬 Study Design and Workflow

1.PNH cell model construction

    CRISPR/RNP-mediated PIGA knockout (PIGA-KO) in THP-1 cells
    
2.Bulk transcriptomic analysis

    Differential expression analysis of mRNAs and miRNAs
    miRNA target prediction and validation filtering
    
3.miRNA–mRNA regulatory network construction

    Identification of key miRNAs and core target genes
    Pathway enrichment analysis
    
4.Single-cell RNA-seq integration

    Integration of PRJNA1061334 and GSE157344 datasets
    Cell type annotation, cell–cell communication analysis
    hdWGCNA module identification
    
5.Drug prediction and validation

    Drug repurposing using EpiMed, CMap, and DGIdb
    Molecular docking and MD simulations
6.Conceptual integration

    Omics → miRNA → Target genes → Drugs → Expected phenotypes
    
📁 mRNA/ — Bulk mRNA-seq Analysis

This folder contains scripts and results related to bulk mRNA differential expression and functional enrichment analysis based on the PIGA-KO PNH model.

📁 miRNA/ — miRNA Differential Expression Analysis

This folder includes analyses related to miRNA expression profiling and differential analysis.

📁 miRNA_Target_Gene/ — miRNA–Target Gene Integration

This directory focuses on miRNA target prediction and validation-based filtering, forming the basis of the miRNA–mRNA regulatory network.

📁 scRNAseq_analysis/ — Single-Cell RNA-seq Analysis

This folder contains scripts and outputs related to single-cell RNA-seq integration, cell-type identification, functional analysis, and network analysis.

📁 upstream/ — Upstream Input Gene Lists

This folder contains preprocessed gene/miRNA lists used as inputs for downstream analyses.

