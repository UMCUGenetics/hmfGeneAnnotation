#library(mltoolkit)
library(SPARQL)
library(stringr)
options(stringsAsFactors = F)
base_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfVariantAnnotation/'
meta_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/metadata/'

#========= Misc functions =========#
merge_by_rownames <- function(...){
   df <- merge(..., by = 'row.names')
   
   rownames(df) <- df[,1]
   df[,1] <- NULL
   
   return(df)
}

write.clipboard <- function(data){
   clip <- pipe("pbcopy", "w")                       
   write.table(data, file=clip)                               
   close(clip)
}

read.clipboard <- function(sep = '\t', header = T, ...){
   read.table(pipe("pbpaste"), sep = sep, header = header)
}

#========= SPARQL args =========#
#! Remember to run 'ssh -L 8890:localhost:8890 fedor13' before making a query
endpoint <- "http://localhost:8890/sparql-auth"
auth_options <- curlOptions(userpwd = "lnguyen:yidogyemUkhebTod")


# ## check for unfound genes
# idx_hr_genes_in_hmf <- lapply(hr_genes_in_hmf, function(i){ grep(i, df_hr_genes$string) }) %>% 
#    unlist() %>% unique() %>% 
#    df_hr_genes$gene[.]
# 
# df_hr_genes$gene %>% .[!(. %in% idx_hr_genes_in_hmf)]

#========= Select patients =========#
# hmf_annotation <- read.table('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/metadata/brca_annotation/prelim_output_HMF_DR10-update-fix_20180828.txt', 
#                              header = T, sep = '\t', row.names = 'Sample_ID_1')

hmf_response <- read.table(paste0(meta_dir, '/brca_annotation/hmf_response.txt'))
selected_patients <- rownames(hmf_response)

#========= Retrieve somatic SNV/INDELs from relevant genes and donors in hmf dataset =========#
getSomaticSnvIndels <- function(){
   query_somatic <-
      "
   PREFIX col: <http://sparqling-genomics/table2rdf/Column/>
   
   PREFIX sampleid_pf: <http://sparqling-genomics/Sample/>
   PREFIX gene_pf: <http://sparqling-genomics/Gene/>
   PREFIX filter_pf: <http://sparqling-genomics/FilterItem/>
   PREFIX type_pf: <http://sparqling-genomics/VariantCallType/>
   
   SELECT
   STRAFTER(STR(?sampleid), STR(sampleid_pf:)) AS ?sampleid
   STRAFTER(STR(?gene), STR(gene_pf:)) AS ?gene
   STRAFTER(STR(?filter), STR(filter_pf:)) AS ?filter
   STRAFTER(STR(?type), STR(type_pf:)) AS ?type
   ?worstcodingeffect
   
   FROM <http://hmfpatients/somaticvariant>
   
   WHERE {
   ?row col:sampleid ?sampleid;
   col:gene ?gene;
   col:filter ?filter;
   col:type ?type;
   col:worstcodingeffect ?worstcodingeffect.
   
   FILTER (?sampleid IN (__sampleid_uri__))
   FILTER(?gene IN (__gene_uri__))
   FILTER (?filter = filter_pf:PASS)
   FILTER (?type IN (type_pf:SNP, type_pf:INDEL))
   FILTER (STR(?worstcodingeffect) IN (\"NONSENSE_OR_FRAMESHIFT\",\"MISSENSE\"))
   }
   "
   
   snv_indel_selected_patients <- list()
   n_patients <- length(selected_patients)
   counter <- 0
   for(i in selected_patients){
      counter <- counter + 1
      message('Retrieving data for: ', i, ' | ', counter, '/', n_patients)

      gene_uri <- paste( paste0('gene_pf:', hr_genes_in_hmf), collapse = ',')
      subquery_somatic <- str_replace(query_somatic, '__gene_uri__', gene_uri)

      sampleid_uri <- paste( paste0('sampleid_pf:', i), collapse = ',')
      subquery_somatic <- str_replace(subquery_somatic, '__sampleid_uri__', sampleid_uri)

      out <- SPARQL(url = endpoint, curl_args = auth_options, query = subquery_somatic)$results
      somatic_selected_patients[[i]] <- out
   }

   somatic_selected_patients <- do.call(rbind,somatic_selected_patients)
   rownames(somatic_selected_patients) <- NULL
   
   write.table(somatic_selected_patients,
               paste0(base_dir,'/hrd_gene_annotation/sparql_mc_hr_genes_somatic_snv_indel.txt'),
               sep = '\t')
}

#========= Retrieve (somatic) gene CN and donors in hmf dataset =========#
query_cnv <-
"
PREFIX col: <http://sparqling-genomics/table2rdf/Column/>
PREFIX sampleid_pf: <http://sparqling-genomics/Sample/>
PREFIX gene_pf: <http://sparqling-genomics/Gene/>

SELECT
STRAFTER(STR(?sampleid),\"http://sparqling-genomics/Sample/\") AS ?sampleid
STRAFTER(STR(?gene), \"http://sparqling-genomics/Gene/\") AS ?gene
?mincopynumber
?maxcopynumber
?minminoralleleploidy

FROM <http://hmfpatients/genecopynumber>
WHERE{
?row col:sampleid ?sampleid;
col:gene ?gene;
col:mincopynumber ?mincopynumber;
col:maxcopynumber ?maxcopynumber;
col:minminoralleleploidy ?minminoralleleploidy.

FILTER(?sampleid IN (__sampleid_uri__))
FILTER(?gene IN (__gene_uri__))
}
"

# cnv_selected_patients <- list()
# n_patients <- length(selected_patients)
# counter <- 0
# for(i in selected_patients){
#    counter <- counter + 1
#    message('Retrieving data for: ', i, ' | ', counter, '/', n_patients)
#
#    gene_uri <- paste( paste0('gene_pf:', hr_genes_in_hmf), collapse = ',')
#    sub_query_cnv <- str_replace(query_cnv, '__gene_uri__', gene_uri)
#
#    sampleid_uri <- paste( paste0('sampleid_pf:', i), collapse = ',')
#    sub_query_cnv <- str_replace(sub_query_cnv, '__sampleid_uri__', sampleid_uri)
#
#    out <- SPARQL(url = endpoint, curl_args = auth_options, query = sub_query_cnv)$results
#    cnv_selected_patients[[i]] <- out
# }
#
# cnv_selected_patients <- do.call(rbind,cnv_selected_patients)
# rownames(cnv_selected_patients) <- NULL

# write.table(cnv_selected_patients,
#             paste0(base_dir,'/hrd_gene_annotation/sparql_mc_hr_genes_cnv.txt'),
#             sep = '\t')

