#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

#========= Purple tables (purity, gene CNV) =========#
subsetGeneCnv_R=$ROOT_DIR/scripts/preProcHmfOutput/subsetGeneCnv.R

procPurpleOutput(){
	out_dir=$1
	sample_name=$2
	purple_purity=$3
	purple_gene_cnv=$4
	interactive=${5:-0}

	job_dir=$out_dir/jobs/; mkdir -p $job_dir; cd $job_dir

	ppo_job=${job_dir}/ppo_${sample_name}.job
	if [[ -f $ppo_job ]]; then rm $ppo_job; fi
	touch $ppo_job

	log_dir=$out_dir/log/; mkdir -p $log_dir

	ppo_cp_purity_done=$log_dir/ppo.cp_purity.done
	if [[ -f $ppo_cp_purity_done ]]; then rm $ppo_cp_purity_done; fi

	ppo_ss_gene_cnv_done=$log_dir/ppo.ss_gene_cnv.done
	if [[ -f $ppo_ss_gene_cnv_done ]]; then rm $ppo_ss_gene_cnv_done; fi

	echo "echo \> Copying purple purity to output dir..." > $ppo_job
	echo "cp $purple_purity $out_dir/${sample_name}.purple.purity && touch $ppo_cp_purity_done" >> $ppo_job

	echo "echo \> Retrieving ENSEMBL gene ids for gene CNV table and subsetting by selected genes..." >> $ppo_job
	echo "guixr load-profile ~/.guix-profile --<< EOF" >> $ppo_job
	echo "Rscript $subsetGeneCnv_R $purple_gene_cnv $out_dir/${sample_name}.purple.gene.cnv && touch $ppo_ss_gene_cnv_done" >> $ppo_job
	echo "EOF" >> $ppo_job

	if [[ $interactive -eq 1 ]]; then 
		sh $ppo_job
	fi
}

# procPurpleOutput \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020493R_CPCT02020493T/ \
# 'CPCT02020493R_CPCT02020493T' \
# /hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170416_HMFregCPCT_FR12245150_FR14064764_CPCT02020493/CPCT02020493T.purple.purity \
# /hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170416_HMFregCPCT_FR12245150_FR14064764_CPCT02020493/CPCT02020493T.purple.gene.cnv \
# 1

#========= Germline/somatic vcfs =========#
vcf2varsigTable_sh=$ROOT_DIR/scripts/preProcHmfOutput/vcf2varsigTable.sh

vcf2varsigTableExec(){
	out_dir=$1
	sample_name=$2
	germ_vcf_in=$3
	som_vcf_in=$4
	interactive=${5:-0}

	#out_dir=${variants_dir}/$sample_name; mkdir -p $out_dir
	#germ_vcf_in=$hmf_data_dir/$dir/$germ_vcf_name
	#som_vcf_in=$hmf_data_dir/$dir/$som_vcf_name

	if [[ -f $out_dir ]]; then 
		echo $out_dir not found
		return 1
	fi

	job_dir=$out_dir/jobs/; mkdir -p $job_dir; cd $job_dir

	v2v_job=${job_dir}/v2v_${sample_name}.job
	if [[ -f $v2v_job ]]; then rm $v2v_job; fi
	touch $v2v_job
	
	log_dir=$out_dir/log/; mkdir -p $log_dir
	v2v_germ_done=$log_dir/v2v.germ.done
	if [[ -f $v2v_germ_done ]]; then rm $v2v_germ_done; fi

	v2v_som_done=$log_dir/v2v.som.done
	if [[ -f $v2v_som_done ]]; then rm $v2v_som_done; fi

	echo source $vcf2varsigTable_sh > $v2v_job

	## Filter multi vcf
	germ_vcf_out=$out_dir/${sample_name}.germ.vcf.gz
	germ_txt_out=$out_dir/varsigs_germ.txt.gz

	echo "echo \> Filtering germline variants and subsetting by selected genes..." >> $v2v_job
	echo filterVcfs $germ_vcf_in $germ_vcf_out \"germline\" >> $v2v_job
	echo extractFields $germ_vcf_out $germ_txt_out \"germline\" >> $v2v_job
	
	## Filter somatic vcf
	som_vcf_out=$out_dir/${sample_name}.som.vcf.gz
	som_txt_out=$out_dir/varsigs_som.txt.gz
	
	echo "echo \> Filtering somatic variants and subsetting by selected genes..." >> $v2v_job
	echo filterVcfs $som_vcf_in $som_vcf_out \"somatic\" >> $v2v_job
	echo extractFields $som_vcf_out $som_txt_out \"somatic\" >> $v2v_job

	## Get clinical annotation
	echo "echo \> Getting clinical annotation for germline variants..." >> $v2v_job
	echo "addClinAnn $germ_txt_out && touch $v2v_germ_done" >> $v2v_job

	echo "echo \> Getting clinical annotation for somatic variants..." >> $v2v_job
	echo "addClinAnn $som_txt_out && touch $v2v_som_done" >> $v2v_job

	if [[ $interactive -eq 1 ]]; then 
		sh $v2v_job
	# else
	# 	qsub -S /bin/bash -cwd -l h_rt=6:00:00 -l h_vmem=16G $v2v_job
	fi
}

# vcf2varsigTableExec \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020493R_CPCT02020493T/ \
# 'CPCT02020493R_CPCT02020493T' \
# /hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170416_HMFregCPCT_FR12245150_FR14064764_CPCT02020493/170416_HMFregCPCT_FR12245150_FR14064764_CPCT02020493.annotated.vcf.gz \
# /hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170416_HMFregCPCT_FR12245150_FR14064764_CPCT02020493/CPCT02020493R_CPCT02020493T_post_processed_v2.2.vcf.gz \
# 1
