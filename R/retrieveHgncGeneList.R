#' Retrieve gene identifiers from HGNC
#'
#' @param hgnc.url URL to HGNC database, where the output contains the columns: HGNC ID, Approved
#' symbol, Previous symbols, Synonyms, ENSEMBL gene id
#' @param out.path Export path
#' @param export.as.rdata Save as RData file?
#' @param verbose Show messages?
#' @param ... Arguments that can be passed to biomaRt::useMart()
#'
#' @return A dataframe of gene identifers from HGNC
#' @export
#'

require(biomaRt)

retrieveHgncGeneList <- function(
   ## excludes withdrawn symbols
   hgnc.url='https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit',
   out.path=NULL,
   export.as.rdata=F,
   verbose=T,
   ...
){
   if(verbose){ message('Downloading HGNC gene list...') }
   genes_hgnc <- read.delim(hgnc.url,check.names=F, comment.char='#')
   colnames(genes_hgnc) <- gsub(' ','_',tolower(colnames(genes_hgnc)))
   
   mart <- useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org',...)
   
   empty_ensg <- nchar(genes_hgnc$ensembl_gene_id) == 0
   if(any(empty_ensg) & verbose){
      message(sprintf(
         '%s/%s entries were found without ENSGs. Attempting to retrieve ENSGs using biomaRt...', 
         sum(empty_ensg), nrow(genes_hgnc)
      ))
   }
   
   ## Split dataframe into missing/non-missing ENSGs
   genes_hgnc_with_ensg <- genes_hgnc[!empty_ensg,]
   genes_hgnc_no_ensg <- genes_hgnc[empty_ensg,]
   
   ## Retrieve ENSGs
   biomart_out <- getBM(
      mart=useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org'),
      attributes=c('hgnc_symbol','ensembl_gene_id'),
      filters='hgnc_symbol',
      values=genes_hgnc_no_ensg$approved_symbol,
      verbose=F
   )
   
   new_ensg <- biomart_out$ensembl_gene_id[match(genes_hgnc_no_ensg$approved_symbol, biomart_out$hgnc_symbol)]
   
   if(any(is.na(new_ensg)) & verbose){
      message(sprintf(
         'ENSGs for %s/%s entries could be retrieved using biomaRt...', 
         sum(!is.na(new_ensg)), length(new_ensg)
      ))
   }
   
   ## Return
   new_ensg[is.na(new_ensg)] <- ''
   genes_hgnc_no_ensg$ensembl_gene_id <- new_ensg
   
   GENES_HGNC <- rbind(genes_hgnc_with_ensg, genes_hgnc_no_ensg)
   GENES_HGNC <- GENES_HGNC[match(GENES_HGNC$hgnc_id, genes_hgnc$hgnc_id),]
   
   if(is.null(out.path)){
      return(GENES_HGNC)
   } else {
      message(sprintf('Exporting table to: %s', out.path))
      if(export.as.rdata & verbose){
         save(GENES_HGNC, file=out.path)
      } else {
         write.table(GENES_HGNC, out.path, sep='\t',row.names=F,quote=F)
      }
   }
}






