library(stringr)
options(stringsAsFactors=F)

base_dir <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
if(dir.exists('/Users/lnguyen/')){
   base_dir <- paste0('/Users/lnguyen/', base_dir)
}

###
diplotypes_ss <- read.delim(paste0(base_dir,'/HMF_update/scripts/annotate_genes/hmf_gene_diplotypes_max_brca.txt.gz'))
head(diplotypes_ss)

###
hmf_coords <- diplotypes_ss[,c('sample','a2.chrom','a2.pos','a2.hgvs_c','diplotype_origin')]
colnames(hmf_coords) <- gsub('a2[.]','',colnames(hmf_coords))

hmf_coords <- hmf_coords[grepl('>',hmf_coords$hgvs_c),]

ref_alt <- do.call(rbind,strsplit(
   str_extract_all(hmf_coords$hgvs_c,'\\w>\\w$',simplify=T)[,1],
   '>'
))

hmf_coord_strings <- paste(hmf_coords$chrom, hmf_coords$pos, ref_alt[,1], ref_alt[,2], sep='_')

###
intersect_clinvar_hmfHotpots <- read.delim(paste0(base_dir,'/scripts_main/hmfGeneAnnotation/data/variant_significance/HMF_hotspots/int_BRCA12_CDK12_TP53_clin_hmf.tsv'))

hmfHotpots_coord_strings <- with(intersect_clinvar_hmfHotpots,{
   paste(chrom,pos,ref,alt,sep='_')
})

###
hmf_coords[hmf_coord_strings %in% hmfHotpots_coord_strings,]