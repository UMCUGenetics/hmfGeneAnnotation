#!/bin/bash

module load tabix/1.7 ## tabix contains bgzip module
ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

#--------- Clinvar ---------#
cd $ROOT_DIR/data/variant_significance/clinvar/

bgzip -c clinvar_20181217_ss.txt > clinvar_20181217_ss.txt.gz
tabix -b 2 -e 2 clinvar_20181217_ss.txt.gz

#--------- ENIGMA ---------#
cd $ROOT_DIR/data/variant_significance/enigma/

tail -n +2 enigma_variants_20181221.txt | sort -V | bgzip -c > enigma_variants_20181221.txt.gz
tabix -b 2 -e 2 -S 1 enigma_variants_20181221.txt.gz