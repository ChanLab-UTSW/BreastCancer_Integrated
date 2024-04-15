# BreastCancer_Integrated

<p align="center">
  <img src="https://ars.els-cdn.com/content/image/1-s2.0-S2666379124001800-fx1_lrg.jpg" width="500" title="[Graphical Abstract](https://doi.org/10.1016/j.xcrm.2024.101511)">
</p>

This includes all of the code for the manuscript "A comprehensive single-cell breast tumor atlas defines cancer epithelial and immune cell heterogeneity and interactions predicting anti-PD-1 therapy response" (https://doi.org/10.1016/j.xcrm.2024.101511).



<hr style="border:2px solid gray">


Here is a breakdown of the file tree:

<details>
<summary> 1. **Preprocessing:** Contains the code for setting up and building the datasets, broken down based on technology type (scRNA-seq, bulk RNA-seq).</summary> 

  - **bulkRNA**: Includes subfolders for The Cancer Genome Atlas (TCGA), I-SPY2, and Cancer Cell Line Encyclopedia (CCLE) datasets.
  - **scRNA**: Includes subfolders for the integrated primary breast tumor atlas and Bassez et al. (anti-PD-1-treated) datasets.
  - **spatial**: Includes code for loading of Visium spatial datasets from Wu et al. 2021 and 10x Genomics.
       
</details>

<details>
<summary> 2. **Analysis**: Contains the code to generate all of the main figures in https://doi.org/10.1016/j.xcrm.2024.101511. </summary>

 - **Figure 1**: Creation of integrated primary breast tumor atlas and identification of NK cell subsets
 - **Figure 2**: Analysis of cancer epithelial cell intratumoral transcriptional heterogeneity, focusing on ERBB2/HER2 and TACSTD/TROP2 
 - **Figure 3**: Identification of 10 GEs to define cancer epithelial cell heterogeneity and immune cell interactions
 - **Figure 4**: Development of InteractPrint based on GE-immune interaction and predicting response to anti-PD-1 therapy 
    
</details>


</details>
