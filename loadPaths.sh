#!/bin/bash

ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

#========= Dependencies =========#
JAVA=$ROOT_DIR/dep/jre1.8.0_191/bin/java
SNPSIFT=$ROOT_DIR/dep/snpEff/SnpSift.jar

#========= Data =========#
GENES_BED=$ROOT_DIR/data/gene_selection/genes.bed
CLINVAR_DB=$ROOT_DIR/data/variant_significance/clinvar/clinvar_20181217_ss.txt
ENIGMA_DB=$ROOT_DIR/data/variant_significance/enigma/enigma_variants_20181221.txt
CADD_DB_SNV=/hpc/cog_bioinf/common_dbs/CADD/whole_genome_SNVs_inclAnno.tsv.gz
SCAP_DB=$ROOT_DIR/data/variant_significance/SCAP/scap_v1_0.sorted.txt.gz
MCAP_DB=$ROOT_DIR/data/variant_significance/MCAP/mcap_v1_3.txt.gz

#========= Scripts =========#
## subset databases
subsetCadd_sh=$ROOT_DIR/scripts/prepDb/subsetCadd.sh

## preProcHmfOutput
preProcHmfOutput_dir=$ROOT_DIR/scripts/preProcHmfOutput/

extractVcfFields_sh=$preProcHmfOutput_dir/extractVcfFields.sh
filterVcf_sh=$preProcHmfOutput_dir/filterVcf.sh
procPurpleOutput_sh=$preProcHmfOutput_dir/procPurpleOutput.sh
subsetGeneCnv_R=$preProcHmfOutput_dir/subsetGeneCnv.R

## getAnnotations
getAnnotations_dir=$ROOT_DIR/scripts/getAnnotations

addSigAnn_sh=$getAnnotations_dir/addSigAnn.sh
getClinSig_sh=$getAnnotations_dir/getClinSig.sh
getCaddAnn_py=$getAnnotations_dir/getCaddAnn.py
getCapAnn_py=$getAnnotations_dir/getCapAnn.py

## annotateGenes
annotateGenesExec_sh=$ROOT_DIR/scripts/annotateGenes/annotateGenesExec.sh
detGeneStatuses_R=$ROOT_DIR//scripts/annotateGenes/detGeneStatuses.R

