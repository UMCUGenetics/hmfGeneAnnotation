#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/
in_dir=$base_dir/HMF_update/variants_subset/
out_dir=$base_dir/HMF_update/analysis/hr_gene_def/

source $base_dir/scripts_main/hmfVariantAnnotation/R/mergeGeneStatusTables.sh

## Condensed gene statuses
manifest_path=$out_dir/manifest_gene_statuses.txt
out_path=$out_dir/hmf_gene_statuses.txt
mergeGeneStatusTables $in_dir $out_path $manifest_path

## Expanded gene statuses
manifest_path=$out_dir/manifest_gene_statuses_expanded.txt
out_path=$out_dir/hmf_gene_statuses_expanded.txt
mergeGeneStatusTables $in_dir $out_path $manifest_path '_gene_statuses_expanded.txt'
