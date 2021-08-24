# Autor: Juraj Michalik
# Date: 08/14/2021
# What it is: Script to generate hierarchy and download + extract data from 10X site to create CD8+ HS Atlas
# License: GNU GPL v3

# This script is distributed to ease the building of the Atlas,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# See http://www.gnu.org/licenses/ for more information about license.

# The files to download are property of 10X genomics and subject
# to Creative Commons license. See https://support.10xgenomics.com/ for
# more informations

# Tested on Fedora and Ubuntu distributions

# bash script to create directories and create files
# make dirs

mkdir Initial_Data
mkdir General
mkdir Donors_Datasets

# download all the files into donor data-sets
# these should be stored in Initial_Data directory
cd Initial_Data

# See README for more information about the data
# GEX information

# 3' sequencing
wget https://cf.10xgenomics.com/samples/cell-exp/6.0.0/10k_PBMCs_TotalSeq_B_3p/10k_PBMCs_TotalSeq_B_3p_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_count_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_raw_gene_bc_matrices.tar.gz

# 5' sequencing
wget https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_count_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_count_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_5gex_protein/vdj_v1_hs_pbmc2_5gex_protein_raw_feature_bc_matrix.tar.gz

# large 5' data sets - be warned these are large and they may take long time to download
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor1/vdj_v1_hs_aggregated_donor1_filtered_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor2/vdj_v1_hs_aggregated_donor2_filtered_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor3/vdj_v1_hs_aggregated_donor3_filtered_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor4/vdj_v1_hs_aggregated_donor4_filtered_feature_bc_matrix.tar.gz

# extract them automatically (needs tape archiver)
tar -xzf 10k_PBMCs_TotalSeq_B_3p_raw_feature_bc_matrix.tar.gz --one-top-level
tar -xzf SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_count_raw_feature_bc_matrix.tar.gz --one-top-level
tar -xzf 5k_pbmc_v3_raw_feature_bc_matrix.tar.gz --one-top-level
tar -xzf pbmc_10k_protein_v3_raw_feature_bc_matrix.tar.gz --one-top-level
tar -xzf pbmc8k_raw_gene_bc_matrices.tar.gz --one-top-level

tar -xzf vdj_v1_hs_pbmc3_raw_feature_bc_matrix.tar.gz --one-top-level
tar -xzf sc5p_v2_hs_T_1k_multi_5gex_t_count_raw_feature_bc_matrix.tar.gz --one-top-level
tar -xzf sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_count_raw_feature_bc_matrix.tar.gz --one-top-level
tar -xzf sc5p_v2_hs_PBMC_1k_raw_feature_bc_matrix.tar.gz --one-top-level
tar -xzf vdj_v1_hs_pbmc2_5gex_protein_raw_feature_bc_matrix.tar.gz --one-top-level

tar -xzf vdj_v1_hs_aggregated_donor1_filtered_feature_bc_matrix.tar.gz --one-top-level
tar -xzf vdj_v1_hs_aggregated_donor2_filtered_feature_bc_matrix.tar.gz --one-top-level
tar -xzf vdj_v1_hs_aggregated_donor3_filtered_feature_bc_matrix.tar.gz --one-top-level
tar -xzf vdj_v1_hs_aggregated_donor4_filtered_feature_bc_matrix.tar.gz --one-top-level

# VDJ information

wget https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_t_filtered_contig_annotations.csv
wget https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_all_contig_annotations.csv
wget https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_vdj_t_all_contig_annotations.csv
wget https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_t_filtered_contig_annotations.csv
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_t/vdj_v1_hs_pbmc2_t_filtered_contig_annotations.csv

wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor1/vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv	
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor2/vdj_v1_hs_aggregated_donor2_all_contig_annotations.csv
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor3/vdj_v1_hs_aggregated_donor3_all_contig_annotations.csv
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor4/vdj_v1_hs_aggregated_donor4_all_contig_annotations.csv
