library(SPARQL)
library(stringr)

base_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/'

#========= Get germline SNV/indels from relevant genes and donors in hmf dataset =========#
#--------- sparql params ---------#
#! Remember to run 'ssh -L 8892:localhost:8892 fedor13' before making a query
endpoint <- "http://localhost:8892/sparql-auth"
auth_options <- curlOptions(userpwd = "lnguyen:yidogyemUkhebTod", post = 1L)

hr_exon_coords <- read.table(paste0(base_dir, '/data/hr_gene_exon_coords.txt'), header = T)

#--------- Main ---------#
getGermlineSnvIndels <- function(chrom, start, end, genes, limit=10){
   # i = 1
   # chrom=hr_exon_coords_split[[i]]$chrom
   # start=hr_exon_coords_split[[i]]$start
   # end=hr_exon_coords_split[[i]]$end
   # genes=hr_exon_coords_split[[i]]$gene
   
   mk_filter_string <- function(chrom, start, end){ 
      paste0(
         '?chrom = rdf_biosemantics:',chrom,' && ',
         '?pos >= ',start,' && ',
         '?pos <= ',end
      )
   }
   
   ## If filter_string is too long, query will break
   filter_string <- mk_filter_string(chrom, start, end)
   
   if(length(chrom) > 1){
      filter_string <- paste0('(',
                              paste(filter_string, collapse = ') || ('),
                              ')'
      )
   }
   
   limit_string <- ifelse(!is.null(limit), paste0('LIMIT ', limit), '')
   
   query <- paste0("
   ## Main prefixes
   PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
   PREFIX sg: <http://sparqling-genomics/>
   PREFIX bh: <http://biohackathon.org/resource/faldo#>
   PREFIX vcf2rdf: <http://sparqling-genomics/vcf2rdf/>
   PREFIX vcf2rdf_format: <http://sparqling-genomics/vcf2rdf/FormatItem/>
   PREFIX vcf2rdf_variant: <http://sparqling-genomics/vcf2rdf/VariantCall/>
   
   ## Extra prefixes for formating objects
   PREFIX sg_Sample: <http://sparqling-genomics/Sample/>
   PREFIX rdf_biosemantics: <http://rdf.biosemantics.org/data/genomeassemblies/hg19#>
   PREFIX sg_Sequence: <http://sparqling-genomics/vcf2rdf/Sequence/>
   
   ## Selection variables which are referenced in the FROM statement
   ## Use STRAFTER(STR(), STR()) to format the output
   #SELECT ?vc ?p ?o ## Exploring vcf colnames
   SELECT 
   STRAFTER(STR(?vc), STR(vcf2rdf:)) AS ?variant_id
   STRAFTER(STR(?sample), STR(sg_Sample:)) AS ?sample
   STRAFTER(STR(?chrom), STR(rdf_biosemantics:)) AS ?chrom
   ?pos
   STRAFTER(STR(?ref),STR(sg_Sequence:)) AS ?ref
   STRAFTER(STR(?alt),STR(sg_Sequence:)) AS ?alt
   STRAFTER(STR(?gt),STR(vcf2rdf:)) AS ?gt
   
   
   FROM <http://hmf/germline> {
   ?vc sg:sample ?sample.
   ?vc bh:reference ?chrom.
   ?vc bh:position ?pos.
   ?vc vcf2rdf_variant:REF ?ref.
   ?vc vcf2rdf_variant:ALT ?alt.
   ?vc vcf2rdf_format:GT ?gt.
   ?vc <http://sparqling-genomics/vcf2rdf/FormatItem/GT> ?gt.
   
   ## Gene coord filter
   FILTER (",filter_string,")
   }
   ",limit_string,"
   
   "
   )
   
   out <- SPARQL(url = endpoint, curl_args = auth_options, query = query)$results
   if(nrow(out != 0)){
      out$gene <- unique(genes)
   } else {
      out <- NA
   }
   return(out)
}

counter <- 0
n_exons <- nrow(hr_exon_coords)

#--------- Exec ---------#
# for(i in 1:nrow(hr_exon_coords)){
# #for(i in 1){
#    
#    chrom <- hr_exon_coords[i,]$chrom
#    start <- hr_exon_coords[i,]$start
#    end <- hr_exon_coords[i,]$end
#    gene <- hr_exon_coords[i,]$gene
#    exon_id <- hr_exon_coords[i,]$exon_id
#    
#    counter <<- counter + 1
#    message('Retrieving germline variants for: ', counter,'/',n_exons,' | ',gene,', ',exon_id)
#    
#    out_path <- paste0(base_dir,'/data/snv_indel_germline/split/',gene,'_',exon_id,'.txt')
#    
#    if(!file.exists(out_path)){
#       out <- getGermlineSnvIndels(chrom, start, end, gene, limit=NULL)
#       if(is.data.frame(out)){
#          write.table(out, out_path, sep = '\t')
#       } else {
#          file.create(out_path)
#       }
#    }
# }

#========= Merge germline SNV/indel tables =========#
# germline_snv_indel_files <- list.files(paste0(base_dir,'/data/snv_indel_germline/split/'),full.names=T)
# 
# counter <- 0
# n_files <- length(germline_snv_indel_files)
# germline_snv_indel <- lapply(germline_snv_indel_files, function(i){
#    counter <<- counter + 1
#    message('Reading ', counter,'/',n_files)
#    df <- try(read.table(i, sep = '\t', header = T, stringsAsFactors = F))
#    if(inherits(df, "try-error")){ return(NULL) }
#    else { return(df) }
# })
# 
# germline_snv_indel <- do.call(rbind, germline_snv_indel)
# germline_snv_indel <- germline_snv_indel[order(germline_snv_indel$sample),]

merged_germline_snv_indel_path <- paste0(base_dir,'/data/snv_indel_germline/merged.txt')
#write.table(germline_snv_indel, merged_germline_snv_indel_path)

#========= Get SnpEff variant type from the above variants =========#
germline_snv_indel <- read.table(merged_germline_snv_indel_path, stringsAsFactors = F)

## bind donor hash
donor_hash <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/sparqling_genomics/scripts/donor_hashes.txt')
colnames(donor_hash) <- c('donor_hash','n3_file')
germline_snv_indel$donor_hash <- str_remove(germline_snv_indel$variant_id, '-V\\d+')
germline_snv_indel <- merge(germline_snv_indel, donor_hash, all.x=T, by='donor_hash')
germline_snv_indel$donor_hash <- NULL

## get snpEff variant type
counter <- 0
n_variants <- nrow(germline_snv_indel)
Map(function(n3_file, variant_id){
   counter <<- counter + 1
   cat('\r','Getting variant ',counter,'/',n_variants)
   # message()
   # system(paste0("zcat | ",n3_file," "))
},germline_snv_indel$n3_file[1:10], germline_snv_indel$variant_id[1:10])