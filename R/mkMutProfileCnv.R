#' Create annotations to the HMF CNV table
#'
#' @description Primarily used to determine full gene loss and LOH
#'
#' @param df.cnv HMF *.gene.cnv table
#' @param full.gene.loss.max.max.copy.number The max max_copy_number for a gene to be considered to 
#' be completely lost (deep deletion)
#' @param trunc.max.min.copy.number The max min_copy_number to for a gene to be considered to have
#' a truncation
#' @param loh.max.min.minor.allele.ploidy The maximum min_minor_allele_ploidy for a gene to be
#' considered to have loss of heterozygosity
#' @param rm.non.selected.genes Remove genes not specified by the user in genes.bed
#' @param genes.bed A bed file which contains ENSEMBL gene ids
#'
#' @return The original input dataframe with annotation columns
#' @export
#'

mkMutProfileCnv <- function(
   df.cnv,
   
   ## Cutoffs
   full.gene.loss.max.max.copy.number = CUTOFFS$full.gene.loss.max.max.copy.number, ## full_gene_loss
   trunc.max.min.copy.number = CUTOFFS$trunc.max.min.copy.number,
   loh.max.min.minor.allele.ploidy = CUTOFFS$loh.max.min.minor.allele.ploidy, ## loh

   ## Scoring
   scoring=SCORING_CNV,

   ## Misc
   rm.non.selected.genes=T, genes.bed=NULL,
   verbose=T
){
   if(nrow(df.cnv) == 0){
      if(verbose){ warning('No rows are present in the gene cnv table. Returning NA...') }
      stop(return(NA))
   }

   #df.cnv=input$cnv
   
   #--------- Main ---------#
   out <- cbind(
      df.cnv,
      full_gene_loss=0,
      trunc=0,
      loh=0
   )
   
   if(verbose){ message('Determining the presence of CNV events...') }
   out <- within(out,{
      full_gene_loss[ min_copy_number <= full.gene.loss.max.max.copy.number ] <- 1
      trunc[ max_copy_number <= trunc.max.min.copy.number ] <- 1
      loh[ min_minor_allele_ploidy <= loh.max.min.minor.allele.ploidy ] <- 1
   })
   
   if(verbose){ message('Getting max CNV effect...') }
   ## Get max CNV eff based on ordering in cnv_eff_names_order
   cnv_eff_names_order <- c('full_gene_loss','trunc','loh')
   out_ss <- as.matrix(out[cnv_eff_names_order])
   out$cnv_eff <- cnv_eff_names_order[ max.col(out_ss, ties.method='first') ]
   
   ## Set if all cnv_eff are 0, set cnv_eff to 'none' (instead of 'full_gene_loss')
   out$cnv_eff[ rowSums(out_ss)==0 ] <- 'none' 
   
   if(verbose){ message('Assigning biallele score to CNV events...') }
   cnv_scores <- unname( scoring[match(out$cnv_eff, names(scoring))] )
   cnv_scores <- do.call(rbind, lapply(cnv_scores, function(i){
      if(is.null(i)){ c(0,0) }
      else { i }
   }))
   colnames(cnv_scores) <- c('a1.score','a2.score')
   
   out <- cbind(out, cnv_scores)
   
   #--------- Overwrite HGNC symbol (some hgnc symbols from the purple.gene.cnv file are outdated) ---------#
   if(verbose){ message('Getting HGNC symbols...') }
   out$hgnc_symbol <- ensgToHgncSymbol(out$ensembl_gene_id)
   
   #--------- Remove ENSG not provide by user (in genes bed file) ---------#
   if(rm.non.selected.genes){
      if(is.null(genes.bed)){
         warning('No genes.bed was not specified. Skipping removing non user selected genes')
      } else {
         if(verbose){ message('Removing non user selected genes...') }
         out <- out[out$ensembl_gene_id %in% genes.bed$ensembl_gene_id,]
      }
   }

   return(out)
}
