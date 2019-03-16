#' Retrieve ENSEMBL gene ids from a vector of gene names
#'
#' @description This function will attempt to retrieve ENSEMBL gene ids from a list of gene names
#' from the HGNC database.
#' @param gene.names A character vector of gene names
#' @param genes.hgnc A HGNC genes list as a dataframe
#' @param na.rough.fix Convert NA values to '' ?
#'
#' @return A character vector of ENSEMBL gene ids
#' @export
#'
geneNamesToEnsg <- function(gene.names, genes.hgnc=GENES_HGNC, na.rough.fix=F){
   #gene.names = mut_profile$germ$snpeff_gene

   if(is.null(genes.hgnc)){ stop('Please provide a HGNC genes list') }

   unique_gene_names <- toupper(unique(gene.names))
   gene_ids <- data.frame(
      gene_names = unique_gene_names,
      ensembl_gene_id = genes.hgnc$ensembl_gene_id[match(unique_gene_names, toupper(genes.hgnc$approved_symbol))]
   )

   ## df of matches found by directly querying HGNC (approved) symbol
   gene_ids_match <- gene_ids[
      nchar(gene_ids$gene_names) != 0 ## rows with no gene_names
      & !is.na(gene_ids$ensembl_gene_id) ## rows with match but no ensembl_gene_id
   ,]

   gene_ids_match$ensembl_gene_id[nchar(gene_ids_match$ensembl_gene_id)==0] <- NA ## Set empty strings to NA (these have an ENSG bu to HGNC symbol)

   ## df of non matches
   gene_ids_non_match <- gene_ids[is.na(gene_ids$ensembl_gene_id),]
   #nrow(gene_ids_match) + nrow(gene_ids_non_match)

   ## Create a list of alternative gene names. Names of each list item is the ENSG
   alt_symbols <- apply(
      cbind(
         strsplit(genes.hgnc$previous_symbols,', '),
         strsplit(genes.hgnc$synonyms,', ')
      ),
      1, unlist
   )
   names(alt_symbols) <- genes.hgnc$ensembl_gene_id

   alt_symbols <- alt_symbols[
      nchar(names(alt_symbols))!=0 & unlist(lapply(alt_symbols, length)) != 0
   ] ## Rm items with no corresponding ENSG, or alternative gene names

   alt_symbols <- lapply(alt_symbols, toupper)

   ## Scan alternative gene names
   gene_ids_non_match$ensembl_gene_id <- unlist(lapply(gene_ids_non_match$gene_names, function(i){
      ensembl_gene_id <- NA
      for(j in 1:length(alt_symbols)){
         if( i[1] %in% alt_symbols[[j]] ){
            ensembl_gene_id <- names(alt_symbols[j])
            break
         }
      }
      return(ensembl_gene_id)
   }))

   ## Return
   gene_ids_final <- rbind(gene_ids_match,gene_ids_non_match)

   na_genes <- gene_ids_final[is.na(gene_ids_final$ensembl_gene_id),'gene_names']
   if(length(na_genes)!=0){
      message(sprintf(
         'ENSEMBL gene ids could not be retrieved for %s genes. %s',
         length(na_genes),
         paste(na_genes, collapse=', ')
      ))
   }

   out <- gene_ids_final$ensembl_gene_id[match(gene.names, gene_ids_final$gene_names)]
   names(out) <- gene.names

   if(na.rough.fix){
      out[is.na(out)] <- ''
   }

   return(out)
}

