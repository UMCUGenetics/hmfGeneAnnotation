var_sig_db_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/data/variant_significance/

clinvar_txt=${var_sig_db_dir}/clinvar_20181217.txt
enigma_txt=${var_sig_db_dir}/enigma_variants_20181221.txt

makeVarSigTable (){

	## args
	vcf_in=$1
	file_out=$2
	mode=$3

	#file_out_temp=$(dirname $file_out)/temp.txt

	## write header
	echo -e "chrom\tpos\tref\talt\tsnpeff_eff\tsnpeff_gene\tgene_id\ttranscript_id\thgvs_c\tclinvar_sig\tenigma_sig\ttumor_gt_ref\ttumor_gt_alt\ttumor_ad_ref\ttumor_ad_alt" > $file_out

	## main
	#zcat $vcf_in | grep -vE '^#' | grep 'BRCA2' | while read chrom pos id ref alt qual filter info format s1 s2; do
	zcat $vcf_in | grep -vE '^#' | while read chrom pos id ref alt qual filter info format s1 s2; do


		## Find in variant significance databases
		getVarSigInDb (){
			db_txt=$1
			## id chrom pos ref alt sig; requires this column order in db files
			## 172799 17 41260985 G C Benign
			sig=$(grep -P "${chrom}\t${pos}\t${ref}\t${alt}\t" $db_txt | cut -f6) ## Included terminal \t so that '... G C' does not match with '... G CC'
			
			if [[ ! -z $sig ]]; then echo "$sig";
			else echo "NA";	fi
		}

		## Extract format values from tumor samples and split into columns
		##
		## Multi-vcf example
		## FORMAT  CPCT02450014R   CPCT02450014T
		## GT:AD:DP:GQ:PL  0/1:11,13:24:99:506,0,419       0/1:57,59:116:99:2236,0,2160
		##
		## Somatic vcf example
		## FORMAT  CPCT02450014T
		## GT:AD:DP        0/1:74,31:122
		##

		extractFormat (){
			if [[ $mode == "germline" ]]; then 
				tumor=$s2
			elif [[ $mode == "somatic" ]]; then 
				tumor=$s1
			else
				echo Available modes are: \"somatic\", \"germline\"
			fi

			gt=$(echo $tumor | cut -d ':' -f1)
			ad=$(echo $tumor | cut -d ':' -f2)

			gt1=$(echo $gt | cut -d '/' -f1); gt2=$(echo $gt | cut -d '/' -f2)
			ad1=$(echo $ad | cut -d ',' -f1); ad2=$(echo $ad | cut -d ',' -f2)

			echo -e "${gt1}\t${gt2}\t${ad1}\t${ad2}"
		}

		## vcf lines with multiple variants (i.e. many ALT sequences sep by comma) will be ignored with this grep approach. 
		## These multi-variants are likely subclonal. Difficult to use these for determining gene deficiencies.
		if [[ $(echo $alt | grep -c ,) -ne 1 ]]; then
			
			## Get SNV/indel clinical significance
			clinvar_sig=$(getVarSigInDb $clinvar_txt)
			enigma_sig=$(getVarSigInDb $enigma_txt)

			## Get variant type and gene from SnpEff annotation. 
			## Note: variants outside of gene regions will not have ANN.
			ann=$(echo $info | grep -oE 'ANN=.+')

			#ann_ss=$(echo $ann | cut -d '|' -f 2,4,5,7,10 | sed 's/|/\t/g') ## cut will return fields in original order!!
			## 2: snpeff_eff
			## 4: snpeff_gene
			## 5: gene_id
			## 7: transcript_id
			## 10: hgvs_c

			snpeff_eff=$(echo $ann | cut -d '|' -f 2)
			snpeff_gene=$(echo $ann | cut -d '|' -f 4)
			gene_id=$(echo $ann | cut -d '|' -f 5)
			transcript_id=$(echo $ann | cut -d '|' -f 7)
			hgvs_c=$(echo $ann | cut -d '|' -f 10)

			## Below statements are in case HR gene coords from the bed file include some regions outside genes.
			## Use R later to re-find gene name at location with the bed file if necessary
			if [[ -z $var_gene ]]; then var_gene="NA"; fi
			if [[ -z $var_type ]]; then var_type="NA"; fi

			## Extract format values
			format_values=$(extractFormat)

			## Output final line
			#echo -e "${chrom}\t${pos}\t${ref}\t${alt}\t${snpeff_gene}\t${gene_id}\t${transcript_id}\t${snpeff_eff}\t${hgvs_c}\t${clinvar_sig}\t${enigma_sig}\t${format_values}"
			echo -e "${chrom}\t${pos}\t${ref}\t${alt}\t${snpeff_eff}\t${snpeff_gene}\t${gene_id}\t${transcript_id}\t${hgvs_c}\t${clinvar_sig}\t${enigma_sig}\t${format_values}"
		fi
	done >> $file_out
}

## Test
# vcf_in=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/CPCT02010419R_CPCT02010419T.multi.vcf.gz
# #vcf_in=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/CPCT02010419R_CPCT02010419T_post_processed_v2.2.vcf.gz

# file_out=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/germline_varsigs.txt

# makeVarSigTable $vcf_in $file_out "germline"

# vcf_in=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/CPCT02010419R_CPCT02010419T.multi.vcf.gz
# #vcf_in=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/CPCT02010419R_CPCT02010419T_post_processed_v2.2.vcf.gz

# file_out=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/variants_subset/CPCT02010419R_CPCT02010419T/germline_varsigs.txt

# makeVarSigTable $vcf_in $file_out "germline"
