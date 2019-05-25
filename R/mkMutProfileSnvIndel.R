#' Create annotations to the germline/somatic tables originating from HMF germline/somatic vcfs
#'
#' @description Primarily used to determine
#'
#' @param df.snv.indel A germline/somatic table originating from HMF germline/somatic vcfs
#' @param rm.non.selected.genes Remove genes not specified by the user in genes.bed
#' @param genes.bed A bed file which contains ENSEMBL gene ids
#' @param gene.identifier Can be 'snpeff_gene' or 'ensembl_gene_id'. If 'snpeff_gene', then ENSEMBL
#' gene ids will be retrieved
#' @param keep.only.first.eff Only keep first snpeff_eff (items are separated by '&')?
#' @param verbose Show messages?
#'
#' @return The original input dataframe with annotation columns
#' @export
#'

mkMutProfileSnvIndel <- function(
   df.snv.indel, 
   mode,
   scoring=SCORING_MUT,
   rm.non.selected.genes=T, 
   genes.bed=NULL,
   gene.identifier='ensembl_gene_id',
   keep.only.first.eff=T,
   verbose=T
){
   
   #--------- Sanity checks ---------#
   if(nrow(df.snv.indel)==0){ stop(return(NA)) }
   #if(!(mode %in% c('germline','somatic'))){ stop("Please specify mode: 'germline', 'somatic'") }
   #if(length(tumor.purity) != 1){ stop("tumor.purity must be a single numeric value") }
   
   #df.snv.indel=input$germ

   #--------- Score annotations ---------#
   if(verbose){ message('Getting snpEff scores...') }
   snpeff_eff <- unlist(lapply(df.snv.indel$snpeff_eff, function(i){
      if(grepl('&',i)){ i <- strsplit(i,'&')[[1]][1] }
      return(i)
   }))
   snpeff_score <- unname(scoring$snpeff)[ match(snpeff_eff,names(scoring$snpeff)) ]
   snpeff_score[is.na(snpeff_score)] <- 0 ## If annotation not found in snpeff table, return 0

   if(verbose){ message('Getting ClinVar scores...') }
   clinvar_score <- unname(scoring$clinvar)[ match(df.snv.indel$clinvar_sig,names(scoring$clinvar)) ]
   clinvar_score[is.na(clinvar_score)] <- 0

   if(verbose){ message('Getting ENIGMA scores...') }
   enigma_score <- unname(scoring$enigma)[ match(df.snv.indel$enigma_sig,names(scoring$enigma)) ]
   enigma_score[is.na(enigma_score)] <- 0

   sig_scores <- cbind(snpeff_score, clinvar_score, enigma_score)

   #--------- Calculate max score and which database it came from ---------#
   if(verbose){ message('Calculating max scores...') }
   db_names <- sapply(strsplit(colnames(sig_scores),'_'), `[`, 1, USE.NAMES=F)
   
   max_sig_scores <- t(apply(sig_scores,1,function(i){
      max_score <- max(i)
      if(max_score==0){ max_score_origin <- 'none' }
      else{
         max_score_origin <- db_names[which(i==max_score)]
         max_score_origin <- paste(max_score_origin, collapse=',')
      }
      return(c(max_score, max_score_origin))
   }))
   max_sig_scores <- as.data.frame(max_sig_scores)
   colnames(max_sig_scores) <- c('max_score','max_score_origin')
   max_sig_scores$max_score <- as.integer(max_sig_scores$max_score)
   
   ### Modifiy out ###
   out <- cbind(
      df.snv.indel,
      sig_scores,
      max_sig_scores
   )
   
   #subset(out,ensembl_gene_id=='ENSG00000175279')

   #--------- Get ENSG if if not exist ---------#
   ### Modifiy out ###
   if(gene.identifier=='snpeff_gene'){
      out$ensembl_gene_id <- geneNamesToEnsg(out$snpeff_gene)

      if(sum(is.na(out$ensembl_gene_id))!=0){
         message(sprintf(
            'Could not determine ENSEMBL gene id for %s/%s variants. Removing these...',
            sum(is.na(out$ensembl_gene_id)),
            nrow(out)
         ))
         out <- out[!is.na(out$ensembl_gene_id),]
      }
   }

   #--------- Get HGNC symbol ---------#
   ### Modifiy out ###
   out <- insColAfter(
      out,
      ensgToHgncSymbol(out$ensembl_gene_id),
      after = 'ensembl_gene_id',
      colname = 'hgnc_symbol'
   )

   #--------- Remove ENSG not provide by user (in genes bed file) ---------#
   ### Modifiy out ###
   if(rm.non.selected.genes){
      if(is.null(genes.bed)){
         warning('No genes.bed was not specified. Skipping removing non user selected genes')
      } else {
         if(verbose){ message('Removing non user selected genes...') }
         out <- out[out$ensembl_gene_id %in% genes.bed$ensembl_gene_id,]
      }
   }

   #--------- Deal with multiple sneff effs ---------#
   ### Modifiy out ###
   if(keep.only.first.eff){
      if(verbose){ message('Keeping only first snpeff_eff...') }
      out$snpeff_eff <- gsub('&.+$','',out$snpeff_eff)
   }

   #--------- Export ---------#
   return(out)
}
