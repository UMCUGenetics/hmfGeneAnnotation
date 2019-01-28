#!/bin/bash

rscript=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/R/mergeGeneStatusTables.R

guixr load-profile ~/.guix-profile/ <<EOF
Rscript $rscript
EOF
