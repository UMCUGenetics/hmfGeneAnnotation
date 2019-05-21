#!/bin/bash

module load tabix/1.7 ## tabix contains bgzip module
ROOT_DIR=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/

#--------- Clinvar ---------#
cd $ROOT_DIR/data/variant_significance/HMF_hotspots/

bgzip -c KnownHotspots.tsv > KnownHotspots.tsv.gz
tabix -b 2 -e 2 KnownHotspots.tsv.gz