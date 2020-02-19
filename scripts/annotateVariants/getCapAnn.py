#!/usr/bin/env python

import gzip
import tabix
import os
import subprocess
import argparse

#========= Get path to CADD database =========#
ROOT_DIR='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'

def getPath (var_name):
	string=subprocess.check_output('cat '+ ROOT_DIR + "loadPaths.sh | grep {}".format(var_name), shell=True).strip('\n').split('=')[1]
	string=string.replace('$ROOT_DIR',ROOT_DIR)
	return string

mcap_db_path=getPath('MCAP_DB')
scap_db_path=getPath('SCAP_DB')

# mcap_db=tabix.open(mcap_db_path)
# tabix_pos="{0}:{1}-{1}".format('1', '69091')

# print list(mcap_db.querys(tabix_pos))

#========= Main =========#
def getCapAnn (in_txt_path, out_txt_path):
	## 
	mcap_db=tabix.open(mcap_db_path)
	scap_db=tabix.open(scap_db_path)

	in_txt=gzip.open(in_txt_path, "rb")

	try:
		os.remove(out_txt_path)
	except OSError:
		pass

	out_txt=gzip.open(out_txt_path, "wb")

	## Write header
	out_txt.writelines('cap_score' + '\t' + 'cap_type' + '\n')

	with in_txt:
		next(in_txt) ## skip header row

		for line in in_txt:
			line=line.strip("\n").split("\t")
			chrom=line[0]
			pos=line[1]
			ref=line[2]
			alt=line[3]
			snpeff_eff=line[4]
			
			tabix_pos="{0}:{1}-{1}".format(chrom, pos)
			
			if len(ref)!=1 or len(alt)!=1: ## skip indels
				cap_output=['NA','NA']
			else:
				## make query
				if 'splice' in snpeff_eff:
					tabix_output=list(scap_db.querys(tabix_pos))
				else:
					tabix_output=list(mcap_db.querys(tabix_pos))

				## return row with correct ref/alt combination
				if len(tabix_output)==0: ## deal with no chrom/pos matches
					cap_output=['NA','NA']
				else:
					tabix_row=0
					for i in tabix_output:
						if i[2]==ref and i[3]==alt: 
							tabix_row=i
							break

					if tabix_row==0: ## deal with no ref/alt matches
						cap_output=['NA','NA']
					else:
						## Return correct columns depending on scap/mcap
						cap_score=tabix_row[4].replace(" ", "") ## Remove trailing spaces
						if cap_score=='COMMON': cap_score='0'

						if 'splice' in snpeff_eff: cap_type=tabix_row[5].replace(" ", "")
						else: cap_type='mcap'

						cap_output=[cap_score,cap_type]

			#print cap_output
			out_txt.writelines('\t'.join(cap_output) + '\n')

	in_txt.close()
	out_txt.close()

# in_txt_path='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02050135R_CPCT02050135T/CPCT02050135R_CPCT02050135T.som.txt.gz'
# out_txt_path='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts/getAnnotations/test.txt.gz'
# getCapAnn (in_txt_path, out_txt_path)


#========= Exec =========#
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Return MCAP/SCAP scores')

	parser.add_argument("-i", "--in_txt_path", help="Input txt.gz", required=True)
	parser.add_argument("-o", "--out_txt_path", help="Output txt.gz", required=True)

	args = parser.parse_args()

	getCapAnn(args.in_txt_path, args.out_txt_path)

