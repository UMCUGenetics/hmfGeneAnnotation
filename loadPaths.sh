#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

#========= Dependencies =========#
JAVA=$ROOT_DIR/dep/jre1.8.0_191/bin/java
SNPSIFT=$ROOT_DIR/dep/snpEff/SnpSift.jar

#========= Data =========#
GENES_BED=$ROOT_DIR/data/gene_selection/genes.bed

CLINVAR_DB=$ROOT_DIR/data/variant_significance/clinvar/clinvar_20181217_ss.txt.gz
ENIGMA_DB=$ROOT_DIR/data/variant_significance/enigma/enigma_variants_20181221.txt.gz

#========= Scripts =========#
#--------- preProcHmfOutput ---------#
subsetInputFiles_dir=$ROOT_DIR/scripts/subsetInputFiles/

extractVcfFields_sh=$subsetInputFiles_dir/extractVcfFields.sh
filterVcf_sh=$subsetInputFiles_dir/filterVcf.sh
subsetGeneCnv_R=$subsetInputFiles_dir/subsetGeneCnv.R

#--------- getAnnotations ---------#
annotateVariants_dir=$ROOT_DIR/scripts/annotateVariants/

addSigAnn_sh=$annotateVariants_dir/addSigAnn.sh
getClinSig_sh=$annotateVariants_dir/getClinSig.sh
getClinSig_py=$annotateVariants_dir/getClinSig.py

#--------- Combine all data ---------#
detGeneStatuses_R=$ROOT_DIR/scripts/pipeline/detGeneStatuses_exec.R
detGeneStatuses_ini=$ROOT_DIR/scripts/pipeline/detGeneStatuses_ini.R

