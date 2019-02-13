#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

detGeneStatuses_R=$ROOT_DIR//scripts/annotateGenes/detGeneStatuses.R

annotateGenes(){
	out_dir=$1
	cnv=$2
	germ=$3
	som=$4
	purity=$5
	sample_name=${6:-''}
	interactive=${7:-0}

	job_dir=$out_dir/jobs/; mkdir -p $job_dir; cd $job_dir

	dgs_job=${job_dir}/dgs_${sample_name}.job
	if [[ -f $dgs_job ]]; then rm $dgs_job; fi
	touch $dgs_job

	log_dir=$out_dir/log/; mkdir -p $log_dir
	dgs_done=$log_dir/dgs.done
	if [[ -f $dgs_done ]]; then rm $dgs_done; fi

	echo "echo \> Determining gene statuses..." >> $dgs_job
	echo "guixr load-profile ~/.guix-profile --<< EOF" >> $dgs_job
	echo "Rscript $detGeneStatuses_R $out_dir $cnv $germ $som $purity $sample_name && touch $dgs_done" >> $dgs_job
	echo "EOF" >> $dgs_job

	if [[ $interactive -eq 1 ]]; then 
		sh $dgs_job
	# else
	# 	qsub -S /bin/bash -cwd -l h_rt=1:00:00 -l h_vmem=16G $dgs_job
	fi
}

# out_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020493R_CPCT02020493T/
# annotateGenes \
# $out_dir \
# $out_dir/CPCT02020493R_CPCT02020493T.purple.gene.cnv \
# $out_dir/varsigs_germ.txt.gz \
# $out_dir/varsigs_som.txt.gz \
# $out_dir/CPCT02020493R_CPCT02020493T.purple.purity \
# 1
