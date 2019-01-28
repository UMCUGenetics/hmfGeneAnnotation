#!/bin/bash

rscript=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/R/detGeneStatuses.R
variants_subset_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/variants_subset/

job_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/jobs/detGeneStatuses/
mkdir -p $job_dir
cd $job_dir

for sample_dir in $variants_subset_dir/*; do
	sample_dir=$(realpath $sample_dir)
	#echo $sample_dir
	sample_name=$( echo $(basename $sample_dir) | cut -d '_' -f 2 )
	job_file=$job_dir/dgs_$sample_name

	if [[ -f $sample_dir/gene_mut_profile/done ]]; then
		echo Skipping $sample_dir
	else
		echo "guixr load-profile ~/.guix-profile/ --<<EOF" > $job_file
		echo "Rscript $rscript $sample_dir" >> $job_file
		echo "EOF" >> $job_file
		
		qsub -S /bin/bash -cwd -m ea -l h_rt=1:00:00 -l h_vmem=10G $job_file
		#echo ... $sample_dir
	fi
done
