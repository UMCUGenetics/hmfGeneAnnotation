#!/bin/bash

module load tabix/1.7 ## tabix contains bgzip module

cd /hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/data/variant_significance/MCAP/

## Convert to bzip
gunzip mcap_v1_3.txt.gz
bgzip mcap_v1_3.txt

tabix -b 2 -e 2 mcap_v1_3.txt.gz