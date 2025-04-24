# ASNPB: Aging-related Splicing-derived Neoantigen Prediction in Blood

This repository provides the full analysis pipeline for our manuscript:  
**"Aging-related alternative splicing drives neoantigen emergence revealed by transcriptome analysis of 1,255 human blood samples."**

It includes all essential scripts, data processing steps, figure generation code, and major intermediate results for reproducibility.

---

## 📦 Requirements

To reproduce the analysis, the following environments/tools are required:

- R (≥ 4.1.0)
- Python (≥ 3.8)
- Jupyter Notebook
- Linux (Ubuntu/CentOS recommended)
- Optional: [Git LFS](https://git-lfs.github.com/) for managing large files

---

## 📁 Repository Structure

### `1_Data_Filter/`
Contains scripts for dataset curation and metadata processing from the GEO database.  
Includes standardized filtering criteria and summary statistics of the selected samples.

### `2_DiffExpr_YoungVsOld/`
Differential gene expression analysis comparing the **young** and **old** cohorts.  
Includes:

- Batch effect correction  
- DEG identification (limma/DESeq2)
- GO and KEGG enrichment
- Exploration of links with immunosenescence signatures

### `3_AS_YoungVsOld/`
Analysis of **alternative splicing (AS)** differences between young and old groups.  
Covers:

- Identification of age-associated AS events (via rMATS)
- Event type breakdown (SE, RI, MXE, A3SS, A5SS)
- Gene-level and chromosomal distributions
- AS-associated splicing factor exploration

### `4_Neoantigen_Prediction/`
Complete pipeline for **neoantigen prediction** derived from aging-related AS events:

- Construction of altered splice isoforms
- Binding affinity prediction via NetMHCpan
- Generation of age-stratified peptide libraries (young + middle-aged)
- Validation using proteomic data from elderly blood samples

### `5_Figures/`

Jupyter notebooks used to generate the **main figures** in the manuscript, facilitating full reproducibility.

> ![Overview](/Overview.png)  
> _Overview of all main images_

### `6_Results/`
Key processed result files, such as:

- noCM_list.RData: Filtered list of candidate events not present in cancer-related CM group

- Differential expression matrices, rMATS results, NetMHC predictions, etc.

## 📧 Contact

📬 lishuhan1@stu.kust.edu.cn





