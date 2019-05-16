#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/
source $ROOT_DIR/loadPaths.sh

#========= Main =========#
pipeline (){
	
	#--------- Arguments ---------#
	debug=0
	skip_to_step=0
	xvf_ignore_type=0
	dgs_init_path=''

	while getopts 'o:b:n:p:c:g:s:t:i:k:d:' arg; do
	  case "${arg}" in
	  	o) out_dir=${OPTARG} ;;
	  	b) genes_bed=${OPTARG} ;;
		n) sample_name=${OPTARG} ;;
		p) purity_path=${OPTARG} ;;
		c) gene_cnv_path=${OPTARG} ;;
		g) germ_vcf_path=${OPTARG} ;;
		s) som_vcf_path=${OPTARG} ;;
		t) xvf_ignore_type=${OPTARG} ;;
		i) dgs_init_path=${OPTARG} ;;
		k) skip_to_step=${OPTARG} ;;
		d) debug=${OPTARG} ;;
	    *) printf "Usage: ...\n"; exit 1 ;;
	  esac
	done
	
	mkdir -p $out_dir

	OPTIND=1 ## Reset getopts

	## Make job and log dirs
	job_dir=$out_dir/jobs/; mkdir -p $job_dir; cd $job_dir
	log_dir=$out_dir/logs/; mkdir -p $log_dir

	#--------- Job submitter ---------#
	execJob (){		
		## Default values
		run_time='1:00:00'
		memory='1G'

		## Parse args
		while getopts 'i:o:p:c:t:m:w:' arg; do
		  case "${arg}" in
		  	i) input=${OPTARG} ;;
			o) output=${OPTARG} ;;
		  	p) job_prefix=${OPTARG} ;;
			c) command=${OPTARG} ;;
			t) run_time=${OPTARG} ;;
			m) memory=${OPTARG} ;;
			w) wait_job=${OPTARG} ;;
		    *) printf "Usage: ...\n"; exit 1 ;;
		  esac
		done

		## Sanity checks
		# if [[ -z $input ]]; then echo 'Warning: No input file was specified'; fi
		# if [[ -z $output ]]; then echo 'Warning: No output file was specified'; fi
		
		if [[ -z $command ]]; then echo 'Error: Please specify bash command'; return 1; fi
		if [[ -z $job_prefix ]]; then echo 'Error: Please specify job prefix'; return 1; fi

		job_file=$job_dir/${job_prefix}_${sample_name}.job
		done_file=$log_dir/${job_prefix}.done

		## Print strings to job
		echo "#!/bin/bash" > $job_file
		echo -e "$command && touch $done_file" >> $job_file

		## Replace INPUT/OUTPUT with real paths. 
		sed -i "s|@INPUT|$input|g" $job_file
		sed -i "s|@OUTPUT|$output|g" $job_file

		## Submit
		if [[ ( -f $output || -d $output ) && -f $done_file ]]; then
			echo "SKIPPING: $(basename $job_file). $(basename $output) and $(basename $done_file) exist"
		else
			if [[ $debug -eq 1 ]]; then 
				echo "RUNNING: $(basename $job_file)"
				sh $job_file
			else 
				qsub -S /bin/bash -cwd -l h_rt=$run_time -l h_vmem=$memory -hold_jid $wait_job \
				-N $(basename $job_file) $job_file
			fi
		fi

		## Reset getopts
		OPTIND=1
	}

	#--------- Main ---------#
	purity_out=$out_dir/${sample_name}.purple.purity
	gene_cnv_ss=$out_dir/${sample_name}.purple.gene.cnv

	if [[ $skip_to_step -le 1 ]]; then
		echo -e "\n#========= Subset gene cnv; Copy purity =========#"
		if [[ ! -f $purity_out ]]; then	
			echo 'Copying purple purity file'; cp $purity_path $purity_out
		else
			echo 'SKIPPING: Copying purple purity file'
		fi
		
		execJob -i $gene_cnv_path -o $gene_cnv_ss -p ssgc -c \
"{ 
guixr load-profile ~/.guix-profile --<<EOF
Rscript $subsetGeneCnv_R @INPUT @OUTPUT $genes_bed
EOF
}" ## Braces are required so that && works after EOF
	fi

	
	som_vcf_ss=$out_dir/${sample_name}.som.vcf.gz
	germ_vcf_ss=$out_dir/${sample_name}.germ.vcf.gz
	if [[ $skip_to_step -le 2 ]]; then 
		echo -e "\n#========= Filter vcfs for gene coords =========#"
		execJob -i $som_vcf_path -o $som_vcf_ss -p fvS -m 8G -t 1:00:00 \
		-c "source $filterVcf_sh; filterVcf @INPUT @OUTPUT $genes_bed 'somatic'"

		execJob -i $germ_vcf_path -o $germ_vcf_ss -p fvG -m 8G -t 1:00:00 \
		-c "source $filterVcf_sh; filterVcf @INPUT @OUTPUT $genes_bed 'germline'"
	fi


	som_txt_ss=$out_dir/${sample_name}.som.txt.gz
	germ_txt_ss=$out_dir/${sample_name}.germ.txt.gz
	
	if [[ $xvf_ignore_type -eq 1 ]]; then
		xvfS_mode='ignore'
		xvfG_mode='ignore'
	else
		xvfS_mode='somatic'
		xvfG_mode='germline'
	fi

	if [[ $skip_to_step -le 3 ]]; then
		echo -e "\n#========= Extract relevant vcf fields into txt =========#"
		execJob -i $som_vcf_ss -o $som_txt_ss -p xvfS -w fvS_${sample_name}.job -m 8G -t 1:00:00 \
		-c "source $extractVcfFields_sh; extractVcfFields @INPUT @OUTPUT $xvfS_mode"

		execJob -i $germ_vcf_ss -o $germ_txt_ss -p xvfG -w fvG_${sample_name}.job -m 8G -t 1:00:00 \
		-c "source $extractVcfFields_sh; extractVcfFields @INPUT @OUTPUT $xvfG_mode" 
	fi


	varsig_dir=$out_dir/varsig/; mkdir -p $varsig_dir
	
	clinsig_som_txt=$out_dir/varsig/clinsig_som.txt.gz
	clinsig_germ_txt=$out_dir/varsig/clinsig_germ.txt.gz
	# if [[ $skip_to_step -le 4 ]]; then
	# 	echo -e "\n#========= ClinVar/ENIGMA annotation =========#"
	# 	execJob -i $som_txt_ss -o $clinsig_som_txt -p gvsS -w xvfS_${sample_name}.job \
	# 	-c "source $ROOT_DIR/loadPaths.sh; source $getClinSig_sh; getClinSig @INPUT @OUTPUT"
		
	# 	execJob -i $germ_txt_ss -o $clinsig_germ_txt -p gvsG -w xvfG_${sample_name}.job \
	# 	-c "source $ROOT_DIR/loadPaths.sh; source $getClinSig_sh; getClinSig @INPUT @OUTPUT"
	# fi


	cadd_som_txt=$out_dir/varsig/cadd_som.txt.gz
	cadd_germ_txt=$out_dir/varsig/cadd_germ.txt.gz
	# if [[ $skip_to_step -le 5 ]]; then
	# 	echo -e "\n#========= CADD annotation =========#"
	# 	execJob -i $som_txt_ss -o $cadd_som_txt -p gcaS -w xvfS_${sample_name}.job -t 3:00:00 \
	# 	-c "$getCaddAnn_py -i @INPUT -o @OUTPUT"

	# 	execJob -i $germ_txt_ss -o $cadd_germ_txt -p gcaG -w xvfG_${sample_name}.job -t 3:00:00 \
	# 	-c "$getCaddAnn_py -i @INPUT -o @OUTPUT"
	# fi

	cap_som_txt=$out_dir/varsig/cap_som.txt.gz
	cap_germ_txt=$out_dir/varsig/cap_germ.txt.gz
	# if [[ $skip_to_step -le 6 ]]; then
	# 	echo -e "\n#========= MCAP/SCAP annotation =========#"
	# 	execJob -i $som_txt_ss -o $cap_som_txt -p gCAPaS -w xvfS_${sample_name}.job -t 1:00:00 \
	# 	-c "$getCapAnn_py -i @INPUT -o @OUTPUT"

	# 	execJob -i $germ_txt_ss -o $cap_germ_txt -p gCAPaG -w xvfG_${sample_name}.job -t 1:00:00 \
	# 	-c "$getCapAnn_py -i @INPUT -o @OUTPUT"
	# fi

	gnomad_som_txt=$out_dir/varsig/gnomad_som.txt.gz
	gnomad_germ_txt=$out_dir/varsig/gnomad_germ.txt.gz
	# if [[ $skip_to_step -le 7 ]]; then
	# 	echo -e "\n#========= GNOMAD annotation =========#"
	# 	execJob -i $som_txt_ss -o $gnomad_som_txt -p gGNOMADaS -w xvfS_${sample_name}.job -t 2:00:00 \
	# 	-c "$getGnomadAnn_py -i @INPUT -o @OUTPUT"

	# 	execJob -i $germ_txt_ss -o $gnomad_germ_txt -p gGNOMADaG -w xvfG_${sample_name}.job -t 5:00:00 \
	# 	-c "$getGnomadAnn_py -i @INPUT -o @OUTPUT"
	# fi

	hotspots_som_txt=$out_dir/varsig/hotspots_som.txt.gz
	hotspots_germ_txt=$out_dir/varsig/hotspots_germ.txt.gz
	# if [[ $skip_to_step -le 8 ]]; then
	# 	echo -e "\n#========= Hotspot annotation =========#"
	# 	execJob -i $som_txt_ss -o $hotspots_som_txt -p ghS -w xvfS_${sample_name}.job -t 2:00:00 \
	# 	-c "source $ROOT_DIR/loadPaths.sh; source $detIsHotspotMut_sh; detIsHotspotMut @INPUT @OUTPUT"

	# 	execJob -i $germ_txt_ss -o $hotspots_germ_txt -p ghG -w xvfG_${sample_name}.job -t 5:00:00 \
	# 	-c "source $ROOT_DIR/loadPaths.sh; source $detIsHotspotMut_sh; detIsHotspotMut @INPUT @OUTPUT"


	som_varsig_txt=$varsig_dir/${sample_name}_varsigs_som.txt.gz
	germ_varsig_txt=$varsig_dir/${sample_name}_varsigs_germ.txt.gz
	if [[ $skip_to_step -le 9 ]]; then
		echo -e "\n#=========Merge variant significance with variant txt =========#"
		execJob -o $som_varsig_txt -p mvsS -w "gvsS_${sample_name}.job" -c \
		"paste <(zcat $som_txt_ss) <(zcat $clinsig_som_txt) <(zcat $hotspots_som_txt) | gzip -c > $som_varsig_txt"

		execJob -o $germ_varsig_txt -p mvsG -w "gvsG_${sample_name}.job" -c \
		"paste <(zcat $germ_txt_ss) <(zcat $clinsig_germ_txt) <(zcat $hotspots_germ_txt) | gzip -c > $germ_varsig_txt"
	fi

	
	if [[ $skip_to_step -le 10 ]]; then
		echo -e "\n#========= Determine gene statuses =========#"
		gene_statuses_dir=$out_dir/gene_statuses/; mkdir -p $gene_statuses_dir
		execJob -p dgs -m 8G \
		-w "ssgc_${sample_name}.job,mvsS_${sample_name}.job,mvsG_${sample_name}.job" \
		-o $gene_statuses_dir \
		-i "$gene_cnv_ss $germ_varsig_txt $som_varsig_txt $purity_out $genes_bed $dgs_init_path" \
		-c \
"{ 
guixr load-profile ~/.guix-profile --<<EOF
Rscript $detGeneStatuses_R @OUTPUT @INPUT
EOF
}"
	fi

}

# #========= Exec =========#
# ## base paths
# hmf_data_dir=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/

# base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/
# variants_dir=$base_dir/HMF_update/vcf_subset/

# ## Exec
# out_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020386R_CPCT02020386T/
# #out_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts/preProcHmfOutput/pipeline_test/
# sample_name='CPCT02020386R_CPCT02020386T'

# purity_path=$hmf_data_dir/161108_HMFregCPCT_FR12244761_FR13272009_CPCT02020386/CPCT02020386T.purple.purity
# gene_cnv_path=$hmf_data_dir/161108_HMFregCPCT_FR12244761_FR13272009_CPCT02020386/CPCT02020386T.purple.gene.cnv

# germ_vcf_path=$hmf_data_dir/161108_HMFregCPCT_FR12244761_FR13272009_CPCT02020386/161108_HMFregCPCT_FR12244761_FR13272009_CPCT02020386.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf.gz
# som_vcf_path=$hmf_data_dir/161108_HMFregCPCT_FR12244761_FR13272009_CPCT02020386/CPCT02020386R_CPCT02020386T_post_processed_v2.2.vcf.gz

# pipeline $out_dir $sample_name $purity_path $gene_cnv_path $germ_vcf_path $som_vcf_path

# #========= Submit manifest =========#
# hmf_data_dir=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/

# base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/
# manifest_path=$base_dir/HMF_update/manifest/hmf_file_manifest.txt
# variants_dir=$base_dir/HMF_update/vcf_subset/

# counter=0
# cat $manifest_path | while read sample_name sample_dir germ_vcf_name som_vcf_name gene_cnv_name purity_name; do
# 	counter=$((counter+1))

# 	echo -e "\n########## [$counter] Submitting gene annotation pipeline for $sample_name ##########"
	
# 	out_dir=$variants_dir/$sample_name; mkdir -p $out_dir

# 	#--------- inputs ---------#
# 	bed_path=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/data/gene_selection/genes.bed
# 	purity_path=$hmf_data_dir/$sample_dir/$purity_name
# 	gene_cnv_path=$hmf_data_dir/$sample_dir/$gene_cnv_name

# 	germ_vcf_path=$hmf_data_dir/$sample_dir/$germ_vcf_name
# 	som_vcf_path=$hmf_data_dir/$sample_dir/$som_vcf_name

# 	#--------- submit ---------#
# 	pipeline -o $out_dir -b $bed_path -n $sample_name -p $purity_path -c $gene_cnv_path -g $germ_vcf_path -s $som_vcf_path \
# 	-d 0 -k 7

# 	#if [[ $counter -eq 1 ]]; then break; fi
# done


