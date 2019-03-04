#!/usr/bin/env python

import gzip
import tabix
import os
import subprocess
import argparse

#========= Get path to CADD database =========#
ROOT_DIR='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'
cadd_db_path=subprocess.check_output('cat '+ ROOT_DIR + 'loadPaths.sh | grep CADD_DB_SNV', shell=True).strip('\n').split('=')[1]

#========= Main =========#
def getCaddAnn (in_txt_path, out_txt_path):
	## 
	cadd_db=tabix.open(cadd_db_path)
	in_txt=gzip.open(in_txt_path, "rb")

	try:
		os.remove(out_txt_path)
	except OSError:
		pass

	out_txt=gzip.open(out_txt_path, "wb")

	## Write header
	out_txt.writelines('cadd_phred' + '\n')

	with in_txt:
		next(in_txt) ## skip header row

		for line in in_txt:
			line=line.strip("\n").split("\t")
			chrom=line[0]
			pos=line[1]
			ref=line[2]
			alt=line[3]
			
			tabix_pos="{0}:{1}-{1}".format(chrom, pos)
			tabix_output=list(cadd_db.querys(tabix_pos))

			cadd_phred='NA' ## set output to NA in case no matches are found

			if len(ref)==1 and len(alt)==1: ## skip indels
				for i in tabix_output:
					if i[2]==ref and i[4]==alt: 
						cadd_phred=i[-1]
						break

			## output
			#print "{} {} {} {} {} {} {} {}".format(ref, alt, tabix_entry[0],tabix_entry[1],tabix_entry[2],tabix_entry[4],tabix_entry[-1],tabix_entry[-21])
			#print cadd_phred
			out_txt.writelines(cadd_phred + '\n')

	in_txt.close()
	out_txt.close()

# in_txt_path='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02050135R_CPCT02050135T/CPCT02050135R_CPCT02050135T.som.txt.gz'
# out_txt_path='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts/preProcHmfOutput/test.txt.gz'
# getCaddAnn (in_txt_path, out_txt_path)


#========= Exec =========#
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Return single column table of CADD phred scores')

	parser.add_argument("-i", "--in_txt_path", help="Input txt.gz", required=True)
	parser.add_argument("-o", "--out_txt_path", help="Output txt.gz", required=True)

	args = parser.parse_args()

	getCaddAnn(args.in_txt_path, args.out_txt_path)

