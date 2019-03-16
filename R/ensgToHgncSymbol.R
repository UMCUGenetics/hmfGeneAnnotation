#' Retrieve HGNC symbols from a vector of ENSEMBL gene ids
#'
#' @description This function will attempt to retrieve HGNC gene symbols from a list of gene names
#' from the HGNC database.
#'
#' @param ensembl.gene.ids A character vector of ENSEMBL gene ids
#' @param genes.hgnc A HGNC genes list as a dataframe
#' @param na.rough.fix Convert NA values to '' ?
#'
#' @return A character vector of HGNC symbols
#' @export
#'
ensgToHgncSymbol <- function(ensembl.gene.ids, genes.hgnc=GENES_HGNC, na.rough.fix=F){
   #ensembl.gene.ids <- input$germ$ensembl_gene_id

   if(is.null(genes.hgnc)){ stop('Please provide a HGNC genes list') }

   ## Use unique ids to speed up search
   unique_ensg <- toupper(unique(ensembl.gene.ids))
   unique_ensg <- unique_ensg[nchar(unique_ensg)!=0] ## Rm unknown ensg ids

   gene_ids <- data.frame(
      ensembl_gene_id=unique_ensg,
      hgnc_symbol=toupper(genes.hgnc$approved_symbol)[match(unique_ensg, genes.hgnc$ensembl_gene_id)]
   )

   na_ensg <- gene_ids[is.na(gene_ids$hgnc_symbol),'ensembl_gene_id']
   if(length(na_ensg)!=0){
      warning(sprintf(
         'HGNC symbols could not be retrieved for %s ENSGs: %s',
         length(na_ensg),
         if(length(na_ensg) <= 10){
            paste(na_ensg, collapse=', ')
         } else {
            paste0(paste(na_ensg[1:10], collapse=', '), '...')
         }

      ))
   }

   ## Convert NA's to empty strings
   if(na.rough.fix){
      gene_ids$hgnc_symbol[is.na(gene_ids$hgnc_symbol)] <- ''
   }

   ## Return the non-unique vector of HGNC ids
   return(
      gene_ids$hgnc_symbol[match(ensembl.gene.ids, gene_ids$ensembl_gene_id)]
   )
}
