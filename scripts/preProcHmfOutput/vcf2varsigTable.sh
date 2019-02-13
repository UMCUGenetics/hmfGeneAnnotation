#!/bin/bash

## Load dependencies. Get dir of script, then path to loadPaths.sh
#wd="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

java=$ROOT_DIR/dep/jre1.8.0_191/bin/java
snpSift=$ROOT_DIR/dep/snpEff/SnpSift.jar
addClinAnn_R=$ROOT_DIR/scripts/preProcHmfOutput/addClinAnn.R

#========= Filter vcfs  =========#
## by PASS and selected gene coords. For germline vcf, also remove somatic variants
filterVcfs(){
	vcf_in=$1
	vcf_out=$2
	mode=$3
	GENES_BED=${4:-$ROOT_DIR/data/gene_selection/genes.bed}

	# vcf_in=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/180426_HMFregCPCT_FR13997274_FR16982076_CPCT02020719/180426_HMFregCPCT_FR13997274_FR16982076_CPCT02020719.annotated.vcf.gz
	# vcf_out=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020719R_CPCT02020719T/CPCT02020719R_CPCT02020719T.germ.vcf.gz
	# mode='germline'
	# bed=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/data/gene_selection/genes.bed

	if [[ $mode == "somatic" ]]; then
		filt_string="(FILTER='PASS')"
	elif [[ $mode == "germline" ]]; then
		filt_string="(FILTER='PASS') & !(GEN[0].GT='0/0') & !(GEN[0].GT='./.')"
	else 
		echo Available modes are: \"somatic\", \"germline\"
	fi

	zcat $vcf_in | \
	$java -jar $snpSift intervals $GENES_BED | \
	$java -jar $snpSift filter "$filt_string" | \
	gzip -c > $vcf_out
}

# filterVcfs \
# /hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170320_HMFregCPCT_FR13274499_FR10244814_CPCT02020459/170320_HMFregCPCT_FR13274499_FR10244814_CPCT02020459.annotated.vcf.gz \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.vcf.gz \
# 'germline'

#========= Convert filtered vcf to table =========#
extractFields(){
	vcf_in=$1
	txt_out=$2
	mode=$3

	# vcf_in=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020719R_CPCT02020719T/CPCT02020719R_CPCT02020719T.germ.vcf.gz
	# txt_out=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020719R_CPCT02020719T/CPCT02020719R_CPCT02020719T.germ.txt.gz
	# mode='germline'

	if [[ $mode == "somatic" ]]; then
		tumor_col=0
	elif [[ $mode == "germline" ]]; then
		tumor_col=1
	else 
		echo Available modes are: \"somatic\", \"germline\"
	fi

	$java -jar $snpSift extractFields $vcf_in \
		CHROM POS REF ALT \
		ANN[0].EFFECT \
		ANN[0].GENE \
		ANN[0].GENEID \
		ANN[0].FEATUREID \
		ANN[0].HGVS_C \
		GEN[$tumor_col].GT \
		GEN[$tumor_col].AD[0] GEN[$tumor_col].AD[1] | \
	sed -e "1s/.*/chrom\tpos\tref\talt\tsnpeff_eff\tsnpeff_gene\tensembl_gene_id\tensembl_transcript_id\thgvs_c\ttumor_gt\ttumor_ad_ref\ttumor_ad_alt/" | \
	gzip -c > $txt_out
}

# extractFields \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.vcf.gz \
# /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.txt.gz \
# 'germline'

#========= Add annotation from clinvar and enigma =========#
addClinAnn(){
	variants_txt=$1
guixr load-profile ~/.guix-profile --<< EOF
Rscript $addClinAnn_R $variants_txt
EOF
}

#addClinAnn /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02020459R_CPCT02020459T/CPCT02020459R_CPCT02020459T.germ.txt.gz.test


