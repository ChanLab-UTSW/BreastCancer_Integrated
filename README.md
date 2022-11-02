# BreastCancer_Integrated


This includes all of the code for the manuscript 

*A comprehensive single-cell breast tumor atlas defines cancer epithelial and immune cell heterogeneity and interactions predicting anti-PD-1 therapy response*

published at the following link: https://www.biorxiv.org/content/10.1101/2022.08.01.501918v3.full (pre-print)



<hr style="border:2px solid gray">


Here is a breakdown of the file tree:

<details>
<summary>1. **Preprocessing**: contains the code for setting up and building the datasets.</summary> 
    <p>
    (1) Breakdown based on **technology type** (scRNA-seq, bulk RNA-seq, RPPA).
    (2) Then broken down by **name of dataset** ("scRNA" includes subfolders for the main primary atlas and Bassez's dataset, for instance, while "bulkRNA" includes subfolders for TCGA, iSPY, etc).
    </p>
    
</details>

<details>
<summary>2. **Analysis**: contains the code to generate all of the main and supplementary figures. Subfolders are separated by main figure.</summary>
    <p>
    a. **Figure 1**: Initial ERBB2/TACSTD2 exploratory analysis and main primary atlas QC
    b. **Figure 2**: NK subset and rNK analyses
    c. **Figure 3**: ITTH score generation and downstream analyses
    d. **Figure 4**: Validation on Bassez and iSPY
    </p>
</details>

<details>
<summary>3. **ReviewerAdditions**: contains the code for running completely new analyses/datasets suggested by the reviewers (Metabric addition, Gambardella et al. cell line dataset, etc.)</summary>

</details>
