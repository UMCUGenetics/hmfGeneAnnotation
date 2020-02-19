#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/
source $ROOT_DIR/loadPaths.sh

# #========= Args =========#
# DEBUG=false
# ## GENES_BED already has a default path from loadPaths.sh

# while [[ $# -gt 0 ]]
# do
# key="$1"
# case $key in
# 	-o|--out_dir)
# 	OUT_DIR="$2"; shift; shift;;

# 	-n|--sample_name)
# 	SAMPLE_NAME="$2"; shift; shift;;

# 	-c|--gene_cnv)
# 	GENE_CNV="$2"; shift; shift;;

# 	-s|--som_vcf)
# 	SOM_VCF="$2"; shift; shift;;

# 	-g|--germ_vcf)
# 	GERM_VCF="$2"; shift; shift;;

# 	-b|--genes_bed)
# 	GENES_BED="$2"; shift; shift;;

# 	-i|--dgs_ini)
# 	DGS_INI="$2"; shift; shift;;

# 	--default)
# 	DEBUG=true; shift;;
# 	*)    # unknown option

# esac
# done

# if [[ ! -d $OUT_DIR ]]; then
# 	echo "Error: Parent dir does not exist: $OUT_DIR"
# 	exit 1
# fi

#echo $GENES_BED

#========= Init =========#
## Test args
DGS_INI=$detGeneStatuses_ini

parent_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CUPs_classifier_data/gene_ann/
SAMPLE_NAME=CPCT02010422T
OUT_DIR=$parent_dir/$SAMPLE_NAME/

som_data_dir=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR-104/data/somatics/161123_HMFregCPCT_FR12244764_FR13275498_CPCT02010422/
GENE_CNV=$som_data_dir/CPCT02010422T.purple.cnv.gene.tsv
SOM_VCF=$som_data_dir/CPCT02010422T.purple.somatic.vcf.gz
GERM_VCF=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR-104/data/germline/161123_HMFregCPCT_FR12244764_FR13275498_CPCT02010422/161123_HMFregCPCT_FR12244764_FR13275498_CPCT02010422.annotated.vcf.gz

if [[ ! -d $(dirname $OUT_DIR) ]]; then
	"Error: Cannot make output dir. Parent dir does not exist: $(dirname $OUT_DIR)"
	exit 1
fi

mkdir -p $OUT_DIR
job_dir=$OUT_DIR/jobs/; mkdir -p $job_dir

jobNameToId (){ 
	out=$(squeue --noheader --format %i --name $1)
	if [[ -z $out ]]; then
		echo 1
	else
		echo $out
	fi
}

#========= Subset input files =========#
processing_dir=$OUT_DIR/proc/; mkdir -p $processing_dir

## Output paths
gene_cnv=$processing_dir/${SAMPLE_NAME}.purple.gene.cnv.txt

som_vcf_ss=$processing_dir/${SAMPLE_NAME}.som.vcf.gz
som_txt=$processing_dir/${SAMPLE_NAME}.som.txt.gz

germ_vcf_ss=$processing_dir/${SAMPLE_NAME}.germ.vcf.gz
germ_txt=$processing_dir/${SAMPLE_NAME}.germ.txt.gz

## Make job
subset_input_job=$job_dir/subsetInput_${SAMPLE_NAME}.job
subset_input_done=${subset_input_job}.done

if [[ ! -f $subset_input_done ]]; then

echo "#!/bin/bash

guixr load-profile ~/.guix-profile --<<EOF

Rscript $subsetGeneCnv_R $GENE_CNV $gene_cnv $GENES_BED &&

source $filterVcf_sh &&
filterVcf $SOM_VCF $som_vcf_ss $GENES_BED 'somatic' &&
filterVcf $GERM_VCF $germ_vcf_ss $GENES_BED 'germline' &&

source $extractVcfFields_sh &&
extractVcfFields $som_vcf_ss $som_txt &&
extractVcfFields $germ_vcf_ss $germ_txt &&

touch $subset_input_done

EOF
" > $subset_input_job

sbatch --time=00:30 --mem=4G \
--job-name=$(basename $subset_input_job) --output=${subset_input_job}.o \
$subset_input_job

fi

#========= Annotate variants =========#
clinsig_som_txt=$processing_dir/clinsig_som.txt.gz
clinsig_germ_txt=$processing_dir/clinsig_germ.txt.gz

som_txt_ann=$processing_dir/${SAMPLE_NAME}.som.ann.txt.gz
germ_txt_ann=$processing_dir/${SAMPLE_NAME}.germ.ann.txt.gz

ann_variants_job=$job_dir/annVariants_${SAMPLE_NAME}.job
ann_variants_done=${ann_variants_job}.done

if [[ ! -f $ann_variants_done ]]; then

echo "#!/bin/bash

$getClinSig_py -i $som_txt -o $clinsig_som_txt &&
$getClinSig_py -i $germ_txt -o $clinsig_germ_txt &&

paste <(zcat $som_txt) <(zcat $clinsig_som_txt) | gzip -c > $som_txt_ann &&
paste <(zcat $germ_txt) <(zcat $clinsig_germ_txt) | gzip -c > $germ_txt_ann &&

touch $ann_variants_done
" > $ann_variants_job

sbatch --time=00:30 --mem=1G \
--job-name=$(basename $ann_variants_job) --output=${ann_variants_job}.o \
--dependency=afterany:$(jobNameToId $(basename $subset_input_job)) \
$ann_variants_job

fi

#========= Determine gene biallelic status =========#
gene_statuses_dir=$OUT_DIR/gene_statuses/; mkdir -p $gene_statuses_dir

det_gene_statuses_job=$job_dir/detGeneStatuses_${SAMPLE_NAME}.job
det_gene_statuses_done=${det_gene_statuses_job}.done

if [[ ! -f $det_gene_statuses_done ]]; then

echo "#!/bin/bash
guixr load-profile ~/.guix-profile --<<EOF
Rscript $detGeneStatuses_R $gene_statuses_dir $gene_cnv $germ_txt_ann $som_txt_ann $GENES_BED $DGS_INI &&
touch $det_gene_statuses_done
EOF
" > $det_gene_statuses_job

sbatch --time=00:30 --mem=8G \
--job-name=$(basename $det_gene_statuses_job) --output=${det_gene_statuses_job}.o \
--dependency=afterany:$(jobNameToId $(basename $ann_variants_job)) \
$det_gene_statuses_job

fi

