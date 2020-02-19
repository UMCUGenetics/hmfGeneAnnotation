#!/usr/bin/env Rscript
## Run on hpc
library(devtools)
load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/')

args <- commandArgs(trailingOnly=T)

detGeneStatuses(
   out.dir = args[1],
   cnv.path = args[2],
   germ.path = args[3],
   som.path = args[4],
   genes.bed.path = args[5],
   ini.path = args[6]
)

