## Load dependencies. Get dir of script, then path to loadPaths.sh
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )/loadPaths.sh 

## main
subsetVcfByCoords (){
	#vcf_in=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/180504_HMFregCPCT_FR16981715_FR16983900_CPCT02080224/180504_HMFregCPCT_FR16981715_FR16983900_CPCT02080224.annotated.vcf.gz
	vcf_in=$1
	vcf_out=$2
	mode=$3
	bed=${4:-/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/data/gene_coords/hr_gene_coords_ensembl.bed}

	#temp_vcf=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/bash/test.vcf.gz
	#vcf_out=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/bash/test2.vcf.gz
	
	temp_vcf=$(dirname $vcf_out)/temp.vcf.gz

	echo Keeping only variants in $(basename $bed)
	zcat $vcf_in | $java -jar $snpsift intervals $bed | gzip -c > $temp_vcf

	
	if [[ $mode == "somatic" ]]; then
		echo Keeping only PASS variants
		zcat $temp_vcf | $java -jar $snpsift filter "(FILTER = 'PASS')" | gzip -c > $vcf_out
	elif [[ $mode == "germline" ]]; then
		echo Keeping only PASS variants and removing variants with GT 0/0 in blood sample
		zcat $temp_vcf | $java -jar $snpsift filter \
		"(FILTER = 'PASS') & !(GEN[0].GT = '0/0') & !(GEN[0].GT = './.')" | \
		gzip -c > $vcf_out
	else 
		echo Available modes are: \"somatic\", \"germline\"
	fi
	
	echo Removing temp vcf
	rm $temp_vcf

	echo Done
}

## test
subsetVcfByCoords \
/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/180426_HMFregCPCT_FR13997274_FR16982076_CPCT02020719/180426_HMFregCPCT_FR13997274_FR16982076_CPCT02020719.annotated.vcf.gz \
/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/bash/test.vcf.gz \
"germline"