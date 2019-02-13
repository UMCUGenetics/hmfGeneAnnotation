#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/
mergeGeneStatusTables_R=$ROOT_DIR//scripts/pipeline/mergeGeneStatusTables.R

in_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/
out_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/analysis/hr_gene_def/
out_path=$out_dir/hmf_gene_statuses.txt.gz

guixr load-profile ~/.guix-profile --<< EOF
echo '\nMerging gene_statuses.txt.gz files'
Rscript $mergeGeneStatusTables_R $in_dir $out_dir/hmf_gene_statuses.txt.gz 'gene_statuses.txt.gz'

echo '\nMerging gene_statuses_expanded.txt.gz files'
Rscript $mergeGeneStatusTables_R $in_dir $out_dir/hmf_gene_statuses_expanded.txt.gz 'gene_statuses_expanded.txt.gz'
EOF
