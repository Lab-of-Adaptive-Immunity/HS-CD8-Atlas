#########################################################

# HS CD8+ Atlas Building 

#########################################################


  
  
A readme for constructing a Human Atlas of CD8+ cells.
The data-set in question is used, among others, in 
manuscript written by Tsyklauri et al., 2022.

**References to papers this data set is used in:**

Tsyklauri O, Chadimova T, Niederlova V, et al. Regulatory T cells suppress the formation of potent KLRK1 and IL-7R expressing effector CD8 T cells by limiting IL-2. Elife. 2023;12:e79342. doi:10.7554/eLife.79342
  
*Link*: https://elifesciences.org/articles/79342https://elifesciences.org/articles/79342

#########################################################
## Legal Information                                     
#########################################################

LICENSE: MIT License for provided scripts.

All scripts are distributed to ease the building of the Atlas,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
MIT License for more details.

Used data sets belong to 10X Genomics and are distributed
under Creative Commons Attribution license. See https://support.10xgenomics.com/ 
for detailed information.



#########################################################
## Requirements                                          
#########################################################
  
You need R 4.0.3 with following packages:
* Seurat 4.0.0 Creative Commons Attribution
* rmarkdown
* ggplot
* dplyr
* tibble

You also need tape archiver (tar) for extracting 10X files.

#########################################################
## Building a set                                        
#########################################################

You need to follow these steps:

1. ) Go to page [here](https://github.com/Lab-of-Adaptive-Immunity/HS-CD8-Atlas/edit/master/README.md) on our GitHub (if you're reading the readme you should be already here) and clone following:  
 git clone https://github.com/Lab-of-Adaptive-Immunity/HS-CD8-Atlas/
 
2. ) Go to downloaded directory:  
  cd HS-CD8-Atlas

3. ) Run bash script  Create_hierarchy_for_analysis.sh:  
    
    bash Create_hierarchy_for_analysis.sh   
    
    This should download all necessary directories and extract them.  
    All VDJ files will be also downloaded.  
  
4. ) On Rstudio or just in terminal, run following files:

    - CD8_HS_Atlas_small_datasets.Rmd
    - CD8_HS_Atlas_Donor_1.Rmd
    - CD8_HS_Atlas_Donor_2.Rmd
    - CD8_HS_Atlas_Donor_3.Rmd
    - CD8_HS_Atlas_Donor_4.Rmd
  
    In Rstudio, for each file use option knit or jsut run all chunks.
    For command line, open R and use:
  
    rmarkdown::render("file.Rmd", "output.html")
  
    where 'file.Rmd' is one of above files and 'output.html' is the name of knitted output file.
    You can name the output however you want, the file is purely informative.
  
5. ) Knit 'Donor_10X_Integration_2500_no_RBS.Rmd' in Rstudio (alternatively,
    run rmarkdown::render("Donor_10X_Integration_2500_no_RBS.Rmd", "Donor_10X_Integration_2500_no_RBS.html")) 
    
At the end you should have 'Donors_Integrated_with_MAIT_final.rds' and 
'Donors_Integrated_no_MAIT_final.rds' files in Donors_Dataset in cloned directory,
among other things. The former contains annotated data set with MAIT cells,
the latter does not have them. These are two data sets used in our paper.



#########################################################
## List of employed data sets and where to get them      
#########################################################

This is just a list of used data sets. You can and should use script 'Create_hierarchy_for_analysis.sh' to download them automatically.
The first link points to the page of the data set. The second link is the link to directly download raw count matrix.

LICENSE: these data-sets belong to 10X Genomics and are distributed under Creative Commons license. See https://support.10xgenomics.com/ 
for more informations.

Single Cell Gene Expression:

* 3' - 1: 	https://support.10xgenomics.com/single-cell-gene-expression/datasets/6.0.0/10k_PBMCs_TotalSeq_B_3p  
			https://cf.10xgenomics.com/samples/cell-exp/6.0.0/10k_PBMCs_TotalSeq_B_3p/10k_PBMCs_TotalSeq_B_3p_raw_feature_bc_matrix.tar.gz
* 3' - 2:		https://support.10xgenomics.com/single-cell-gene-expression/datasets/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex  
			https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_count_raw_feature_bc_matrix.tar.gz
* 3' - 3:		https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_v3  
			https://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_raw_feature_bc_matrix.tar.gz
* 3' - 4:		https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3  
			https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_raw_feature_bc_matrix.tar.gz
* 3' - 5:		https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k  
			https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_raw_gene_bc_matrices.tar.gz

Single Cell Immune Profiling:

* 5' - 1:		https://support.10xgenomics.com/single-cell-vdj/datasets/3.1.0/vdj_v1_hs_pbmc3  
			https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_raw_feature_bc_matrix.tar.gz
* 5' - 2:		https://support.10xgenomics.com/single-cell-vdj/datasets/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t  
			https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_count_raw_feature_bc_matrix.tar.gz
* 5' - 3:		https://support.10xgenomics.com/single-cell-vdj/datasets/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t  
			https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_count_raw_feature_bc_matrix.tar.gz			
* 5' - 4:		https://support.10xgenomics.com/single-cell-vdj/datasets/4.0.0/sc5p_v2_hs_PBMC_1k  
			https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_raw_feature_bc_matrix.tar.gz
* 5' - 5:		https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_hs_pbmc2_5gex_protein  
			https://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_5gex_protein/vdj_v1_hs_pbmc2_5gex_protein_raw_feature_bc_matrix.tar.gz

Donors:

* Donor 1:	https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.2/vdj_v1_hs_aggregated_donor1  
			https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor1/vdj_v1_hs_aggregated_donor1_filtered_feature_bc_matrix.tar.gz
* Donor 2:	https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.2/vdj_v1_hs_aggregated_donor2  
			https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor2/vdj_v1_hs_aggregated_donor2_filtered_feature_bc_matrix.tar.gz
* Donor 3:	https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.2/vdj_v1_hs_aggregated_donor3  
			https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor3/vdj_v1_hs_aggregated_donor3_filtered_feature_bc_matrix.tar.gz
* Donor 4:	https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.2/vdj_v1_hs_aggregated_donor4  
			https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor4/vdj_v1_hs_aggregated_donor4_filtered_feature_bc_matrix.tar.gz

VDJ - 5':

* 5' - 1:		https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_t_filtered_contig_annotations.csv
* 5' - 2: 	https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_all_contig_annotations.csv
* 5' - 3: 	https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_vdj_t_all_contig_annotations.csv
* 5' - 4:		https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_t_filtered_contig_annotations.csv
* 5' - 5:		https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_hs_pbmc2_t
			https://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_t/vdj_v1_hs_pbmc2_t_filtered_contig_annotations.csv

VDJ - Donors:

* Donor 1:	https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor1/vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv	
* Donor 2:	https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor2/vdj_v1_hs_aggregated_donor2_all_contig_annotations.csv
* Donor 3:	https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor3/vdj_v1_hs_aggregated_donor3_all_contig_annotations.csv
* Donor 4:	https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor4/vdj_v1_hs_aggregated_donor4_all_contig_annotations.csv



#########################################################
## MISCELLANEOUS                                         
#########################################################

Author: Juraj Michalik  
Date: 2021-08-24  
e-mail: juraj.michalik@img.cas.cz  

In case of any questions, feel free to use above e-mail adress.
