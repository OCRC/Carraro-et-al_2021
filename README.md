This repository contains open-source code for reproducing the analysis of single-cell RNA sequencing data exposed in the manuscript 'Arcuate nucleus overexpression of Nhlh2 reduces body mass and attenuates obesity-associated anxiety/depression-like behavior', by Carraro et al., [published at the Journal of Neuroscience](https://doi.org/10.1523/JNEUROSCI.0222-21.2021).

In this study, we harnessed the data from [Campbell et al](https://doi.org/10.1038/nn.4495), consisting of ~20,000 single-cell transcriptomes from the arcuate nucleus and the median eminence (Arc-ME), and used it to identify Nhlh2 transcriptional targets.

# Retrieve data
  We retrieved gene expression and meta-data matrices from SingleCellPortal. You'll need to sign in with a Google account to continue:
  - [Gene expression data](https://singlecell.broadinstitute.org/single_cell/data/public/SCP97/a-molecular-census-of-arcuate-hypothalamus-and-median-eminence-cell-types?filename=expression.txt.gz)
  - [Metadata](https://singlecell.broadinstitute.org/single_cell/data/public/SCP97/a-molecular-census-of-arcuate-hypothalamus-and-median-eminence-cell-types?filename=meta.txt)

  Alternatively, download directly from the UNIX (Linux, MacOS) terminal:
  ```
  # Gene expression data
  wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93374/suppl/GSE93374_Merged_all_020816_DGE.txt.gz -O expression.txt.gz
  
  # Metadata
  wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93374/suppl/GSE93374_cell_metadata.txt.gz -O meta.txt.gz
  ```
  
  Please note that both these files must be in the same directory as the scripts (our working directory, or where you cloned the repo to), and should be named 'expression.txt.gz' and 'meta.txt.gz'.

# Single-cell analysis with dbMAP

  For this step, you'll need to have [dbMAP](https://github.com/davisidarta/dbMAP) installed, along with [Seurat](https://satijalab.org/seurat) for single-cell analysis in R.
  Main findings can be reproduced by executing the following line in the terminal:
  
  ```
  > R Zanesco_et_al_2022_Campbell_data_analysis.R
  ```
  Alternatively, you may open the `Zanesco_et_al_2022_Campbell_data_analysis.R` script in your favorite IDE (i.g. RStudio) and interactively navigate the code.
  
  
# Running Arboreto and PySCENIC

  For this step, you should download some databases: https://pyscenic.readthedocs.io/en/latest/installation.html#auxiliary-datasets
  
  To run arboreto, execute the following snippets in the terminal: 
  
  
```
arboreto_with_multiprocessing.py \
    CampNeurons.tsv \
    mm_mgi_tfs.txt \
    --method grnboost2 \
    --output adj.tsv \
    --num_workers 20 \
    --seed 777
```
    
    
 And
 
 ```
pyscenic ctx adj.tsv --annotations_fname SCENIC/resources/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname CampNeurons.loom --output arc_reg.csv --num_workers 64 --min_genes 10  databases/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather databases/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather databases/mm9-tss-centered-10kb-10species.mc9nr.feather databases/mm9-tss-centered-5kb-10species.mc9nr.feather databases/mm9-500bp-upstream-10species.mc9nr.feather
```

Additional code that could not be easily pipelined is also available at this repository.


Please send any questions and/or bug reports to davisidarta[at]fcm.unicamp.br


# Citation

  If these scripts are useful for you research, we ask that you cite: 
  
  Carraro, R; Nogueira, G; Sidarta-Oliveira, D _et al_. Arcuate nucleus overexpression of Nhlh2 reduces body mass and attenuates obesity-associated anxiety/depression-like behavior'. The Journal of Neuroscience.
  
  If dbMAP is useful for your research, please cite:
    
  Sidarta-Oliveira, D and Velloso, L. Comprehensive Visualization of High-Dimensional Single-Cell Data With Diffusion-Based Manifold Approximation and Projection (dbMAP).
  Available at SSRN: https://ssrn.com/abstract=3582067 or http://dx.doi.org/10.2139/ssrn.3582067
  
