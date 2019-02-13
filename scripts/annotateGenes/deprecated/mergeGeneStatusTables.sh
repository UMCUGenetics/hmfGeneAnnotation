#!/bin/bash

mergeGeneStatusTables (){
	in_dir=$1
	out_path=$2
	manifest_path=$3
	suffix=${4:-_gene_statuses.txt}

	## Get paths to gene status files
	if [[ -f $manifest_path ]]; then 
		echo "File manifest exists. Skipping retrieving paths to files"
	else
		echo "Retrieving paths to files..."
		echo -e "path\tsample_name" > $manifest_path
		#suffix='_gene_statuses.txt'
		for i in $in_dir/*/gene_mut_profile/*$suffix; do
			sample_name=$(basename $i $suffix)
			echo -e "${i}\t${sample_name}"
		done >> $manifest_path
	fi

	## Merge tables
	rscript=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/R/mergeGeneStatusTables.R
guixr load-profile ~/.guix-profile/ <<EOF
Rscript $rscript $manifest_path $out_path
EOF

}

## Test
# base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/

# in_dir=$base_dir/HMF_update/variants_subset/
# manifest_path=$base_dir/HMF_update/analysis/hr_gene_def/manifest_gene_statuses.txt
# out_path=$base_dir/HMF_update/analysis/hr_gene_def/hmf_gene_statuses.txt

# mergeGeneStatusTables $in_dir $out_path $manifest_path

# mergeGeneStatusTables $base_dir/HMF_update/variants_subset/ \
# $base_dir/HMF_update/analysis/hr_gene_def/hmf_gene_statuses.txt \
# $base_dir/HMF_update/analysis/hr_gene_def/manifest_gene_statuses.txt
