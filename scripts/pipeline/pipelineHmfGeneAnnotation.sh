#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

preProcHmfOutput_sh=$ROOT_DIR/scripts/preProcHmfOutput/preProcHmfOutput.sh; source $preProcHmfOutput_sh
annotateGenesExec_sh=$ROOT_DIR/scripts/annotateGenes/annotateGenesExec.sh; source $annotateGenesExec_sh

#========= Submit manifest =========#
hmf_data_dir=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/
manifest_path=$base_dir/HMF_update/manifest/hmf_file_manifest.txt
variants_dir=$base_dir/HMF_update/vcf_subset/

counter=0
cat $manifest_path | while read sample_name sample_dir germ_vcf_name som_vcf_name gene_cnv_name purity_name; do
	counter=$((counter+1))

	echo -e '\n'
	echo "###### [$counter] Submitting gene annotation pipeline for $sample_name ######"

	out_dir=$variants_dir/$sample_name; mkdir -p $out_dir

	#--------- subset gene cnv, copy purity ---------#
	purity_path=$hmf_data_dir/$sample_dir/$purity_name
	gene_cnv_path=$hmf_data_dir/$sample_dir/$gene_cnv_name

	purity_out=$out_dir/${sample_name}.purple.purity
	gene_cnv_out=$out_dir/${sample_name}.purple.gene.cnv
	
	procPurpleOutput $out_dir $sample_name $purity_path $gene_cnv_path
	ppo_job_name=ppo_${sample_name}.job
	ppo_job=$out_dir/jobs/$ppo_job_name

	if [[ ! -f $purity_path && ! -f $gene_cnv_out ]]; then
		qsub -S /bin/bash -cwd -l h_rt=1:00:00 -l h_vmem=16G -N $ppo_job_name $ppo_job
	else 
		echo "Skipping procPurpleOutput. $(basename $purity_path) and $(basename $gene_cnv_out) exist."
	fi

	#--------- germ/som vcfs to varsig table ---------#
	germ_vcf_path=$hmf_data_dir/$sample_dir/$germ_vcf_name
	som_vcf_path=$hmf_data_dir/$sample_dir/$som_vcf_name

	germ_out=$out_dir/varsigs_germ.txt.gz
	som_out=$out_dir/varsigs_som.txt.gz

	vcf2varsigTableExec $out_dir $sample_name $germ_vcf_path $som_vcf_path
	v2v_job_name=v2v_${sample_name}.job
	v2v_job=$out_dir/jobs/$v2v_job_name
	
	if [[ ! -f $germ_out && ! -f $som_out ]]; then
		qsub -S /bin/bash -cwd -l h_rt=6:00:00 -l h_vmem=16G -hold_jid $ppo_job_name -N $v2v_job_name $v2v_job
	else 
		echo "Skipping vcf2varsigTableExec. $(basename $germ_out) and $(basename $som_out) exist."
	fi

	#--------- Gene annotation ---------#
	gene_statuses_out=$out_dir/gene_statuses.txt.gz
	gene_statuses_expanded_out=$out_dir/gene_statuses_expanded.txt.gz

	annotateGenes $out_dir $gene_cnv_out $germ_out $som_out $purity_out $sample_name
	dgs_job_name=dgs_${sample_name}.job
	dgs_job=$out_dir/jobs/$dgs_job_name

	if [[ ! -f $gene_statuses_out && ! -f $gene_statuses_expanded_out ]]; then
		qsub -S /bin/bash -cwd -l h_rt=1:00:00 -l h_vmem=16G -hold_jid $v2v_job_name -N $dgs_job_name $dgs_job
	else 
		echo "Skipping annotateGenes. $(basename $gene_statuses_out) and $(basename $gene_statuses_expanded_out) exist."
	fi


	#if [[ $counter -eq 1 ]]; then break; fi
done






