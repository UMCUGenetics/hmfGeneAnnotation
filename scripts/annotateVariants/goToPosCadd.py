#!/usr/bin/env python

import tabix
import os
import subprocess

ROOT_DIR='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'
cadd_db_path=subprocess.check_output('cat '+ ROOT_DIR + 'loadPaths.sh | grep CADD_DB_SNV', shell=True).strip('\n').split('=')[1]
cadd_db=tabix.open(cadd_db_path)

#chrom='13'; pos='32954050'
chrom='8'; pos='145738767'


tabix_pos="{0}:{1}-{1}".format(chrom, pos)
tabix_output=list(cadd_db.querys(tabix_pos))

for i in tabix_output:
	print "{} {} {} {} {} {}".format(i[0],i[1],i[2],i[3],i[-1],i[-21])