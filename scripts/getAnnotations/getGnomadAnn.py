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

gnomad_db_dir=getPath('GNOMAD_DB_DIR')

#========= Main =========#
def getGnomadAnn (in_txt_path, out_txt_path):
	## 
	# mcap_db=tabix.open(mcap_db_path)
	# scap_db=tabix.open(scap_db_path)

	in_txt=gzip.open(in_txt_path, "rb")

	try:
		os.remove(out_txt_path)
	except OSError:
		pass

	## Open output file and write header
	out_txt=gzip.open(out_txt_path, "wb")
	out_txt.writelines('gnomad_filter'+'\t'+'gnomad_af'+'\n')	

	with in_txt:
		next(in_txt) ## skip header row

		## gnomad db is split into separate files by chrom. Open the correct chrom file when chrom changes
		chrom_counter=0
		# gnomad_db_path='{0}/gnomad.genomes.r2.1.sites.chr1.vcf.bgz'.format(gnomad_db_dir)
		# gnomad_db=tabix.open(gnomad_db_path)

		for line in in_txt:
			line=line.strip("\n").split("\t")
			chrom=line[0]
			pos=line[1]
			ref=line[2]
			alt=line[3]
			
			if chrom_counter != chrom:
				chrom_counter=chrom
				gnomad_db_path='{0}/gnomad.genomes.r2.1.sites.chr{1}.vcf.bgz'.format(gnomad_db_dir,chrom)
				gnomad_db=tabix.open(gnomad_db_path)

			tabix_pos="{0}:{1}-{1}".format(chrom, pos)
			tabix_output=list(gnomad_db.querys(tabix_pos))

			if len(tabix_output)==0:
				gnomad_out=['NA','NA']
			else:
				tabix_row=0
				for i in tabix_output:
					if i[3]==ref and i[4]==alt: 
						tabix_row=i
						break

				if tabix_row==0: ## deal with no ref/alt matches
					gnomad_out=['NA','NA']
				else:
					gnomad_filter=tabix_row[6]
					gnomad_af=tabix_row[7].split(";")[2].replace('AF=','')

					gnomad_out=[gnomad_filter,gnomad_af]

			#print gnomad_out
			out_txt.writelines('\t'.join(gnomad_out) + '\n')

	in_txt.close()
	out_txt.close()

# in_txt_path='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02050135R_CPCT02050135T/CPCT02050135R_CPCT02050135T.som.txt.gz'
# out_txt_path='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts/getAnnotations/test.txt.gz'
# getGnomadAnn (in_txt_path, out_txt_path)


#========= Exec =========#
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Return GNOMAD annotations')

	parser.add_argument("-i", "--in_txt_path", help="Input txt.gz", required=True)
	parser.add_argument("-o", "--out_txt_path", help="Output txt.gz", required=True)

	args = parser.parse_args()

	getGnomadAnn(args.in_txt_path, args.out_txt_path)

