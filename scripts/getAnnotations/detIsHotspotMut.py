#!/usr/bin/env python

import gzip
import tabix
import os
import subprocess
import argparse

#========= Get path to CADD database =========#
ROOT_DIR='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'

def getPath (grep_string):
	cmd='cat '+ ROOT_DIR + 'loadPaths.sh | grep ' + grep_string
	path=subprocess.check_output('cat '+ ROOT_DIR + 'loadPaths.sh | grep ' + grep_string, shell=True)
	path=path.strip('\n').split('=')[1].split(' ')[0]
	path=path.replace('$ROOT_DIR',ROOT_DIR)
	return path

hotspots_db_path=getPath('^HOTSPOTS_DB')

#print hotspots_db_path

#========= Main =========#
def detIsHotspotMut(in_txt_path, out_txt_path):
	
	#--------- Prep in/out txt files---------#
	in_txt=gzip.open(in_txt_path, "rb")

	try:
		os.remove(out_txt_path)
	except OSError:
		pass

	out_txt=gzip.open(out_txt_path, "wb")
	out_txt.writelines('is_hotspot_mut' + '\n')

	#--------- Load database ---------#
	hotspots_db=tabix.open(hotspots_db_path)

	#--------- Main ---------#
	with in_txt:
		next(in_txt) ## skip header row

		counter=0
		for line in in_txt:
			counter+=1
			#if counter==50: break

			line=line.strip("\n").split("\t")
			chrom=line[0]
			pos=line[1]
			ref=line[2]
			alt=line[3]
			
			tabix_pos="{0}:{1}-{1}".format(chrom, pos)
			tabix_out=list(hotspots_db.querys(tabix_pos))
			is_hotspot_mut=0

			if len(tabix_out)==1:
				i=tabix_out[0]
				if i[2]==ref and i[3]==alt:
					is_hotspot_mut=1
			else:
				for i in tabix_out:
					if i[2]==ref and i[3]==alt:
						is_hotspot_mut=1
						break

			#print clinvar_sig
			out_txt.writelines(str(is_hotspot_mut) +'\n')

# in_txt_path='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/CPCT02050135R_CPCT02050135T/CPCT02050135R_CPCT02050135T.som.txt.gz'
# out_txt_path='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts/getAnnotations/test.txt.gz'

# detIsHotspotMut(in_txt_path, out_txt_path)

#========= Exec =========#
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Return a column containing 0/1 (whether a variant is in the HMF hotspot list)')

	parser.add_argument("-i", "--in_txt_path", help="Input txt.gz", required=True)
	parser.add_argument("-o", "--out_txt_path", help="Output txt.gz", required=True)

	args = parser.parse_args()

	getClinSig(args.in_txt_path, args.out_txt_path)
