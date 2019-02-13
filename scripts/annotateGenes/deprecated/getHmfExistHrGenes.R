library(SPARQL)
library(stringr)
library(magrittr)
options(stringsAsFactors = F)
base_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/'
meta_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/metadata/'

#========= Retrieve all genes in hmf dataset =========#
# query <-
# "
# PREFIX col: <http://sparqling-genomics/table2rdf/Column/>
# 
# SELECT DISTINCT STRAFTER(STR(?gene), \"http://sparqling-genomics/Gene/\")
# FROM <http://hmfpatients/somaticvariant>
# WHERE {
#    ?row col:gene ?gene
# }
# "
# 
# hmf_exist_genes <- SPARQL(url = endpoint, curl_args = auth_options, query = query)
# 
# hmf_exist_genes <- as.character(hmf_exist_genes$results)
# hmf_exist_genes <- sort(hmf_exist_genes)
# 
# write.table(hmf_exist_genes,'/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/scripts/mltoolkit/random_forest_training_v2/hmf_update_annotation/hmf_exist_genes.txt',
#             sep = '\t', col.names = F, row.names = F)

#========= Get gene id's for gene names =========#
hmf_exist_genes <- read.table(paste0(meta_dir, '/exist_genes/hmf_exist_genes.txt'),sep='\t',stringsAsFactors = F)
hmf_exist_genes <- hmf_exist_genes[,1]

#========= Gene lists =========#
#hmf_exist_genes[grep('^RAD.*',hmf_exist_genes)] %>% sort

hr_genes <- read.table(
   paste0(base_dir,'/data/hr_genes.txt'),
   sep='\t', header = T, check.names = F, stringsAsFactors = F
)

# hr_genes_ensembl <- read.table(
#    paste0(base_dir,'/data/gene_coords/hr_gene_coords_ensembl.bed'),
#    sep='\t', header = T, check.names = F, stringsAsFactors = F
# )

hr_gene_aliases <- str_split(hr_genes$string, '[[:space:]]*,[[:space:]]*')
hr_genes_in_hmf <- lapply(hr_gene_aliases, function(i){ i[i %in% hmf_exist_genes] })
# names(hr_genes_in_hmf) <- hr_genes$gene_id

write.table(unlist(hr_genes_in_hmf, use.names = F), 
            paste0(meta_dir, '/exist_genes/hmf_exist_hr_genes.txt'),
            row.names = F, col.names = F, quote = F)