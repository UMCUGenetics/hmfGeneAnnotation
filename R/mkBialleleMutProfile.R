#' Performs pairwise joining of the cnv, germ, and som mut_profile tables
#'
#' @description Rows are joined by ENSEMBL gene id. All rows are kept (natural join).
#'
#' @param mut_profile mut_profile object. A list of of cnv, germ and som tables.
#' @param verbose Show messages?
#'
#' @return A list containing cnv_germ, cnv_som, and germ_som tables
#' @export
#'
mkBialleleMutProfile <- function(mut_profile, verbose=T){

   SEL_COLS <- list(
      all = c('ensembl_gene_id', 'hgnc_symbol'),

      cnv = c('full_gene_loss','loh','cn_break_in_gene','amp'),

      germ = c(
         'chrom','pos','hgvs_c',
         #'ExAC_AC','ExAC_AF','gnomad_filter','gnomad_af',
         'snpeff_eff','clinvar_sig','enigma_sig','max_score','max_score_origin',
         #'cadd_phred', 'cap_score','cap_type',
         'adj_tumor_ad_ref','adj_tumor_ad_alt','alt_exists','ref_loss'
      ),

      som = c(
         'chrom','pos','hgvs_c',
         #'ExAC_AC','ExAC_AF','gnomad_filter','gnomad_af',
         'snpeff_eff','clinvar_sig','enigma_sig','max_score','max_score_origin',
         #'cadd_phred', #'cap_score','cap_type',
         'adj_tumor_ad_ref','adj_tumor_ad_alt','alt_exists','ref_loss'
      )
   )

   mkValidDf <- function(index){
      col_names <- c(SEL_COLS[['all']], SEL_COLS[[index]])

      ## For df in gene_mut_profile with no variants, create df of NAs for errorless cbind
      if(!is.data.frame(mut_profile[[index]])){
         df <- data.frame(matrix( nrow=length(union_gene_ids), ncol=length(col_names) ))
         colnames(df) <- col_names

      } else {
         df <- mut_profile[[index]][c(SEL_COLS[['all']], SEL_COLS[[index]])]
      }

      return(df)
   }

   out <- list()

   ## NOTE: slightly different for each combination
   if(verbose){ message('Merging mutation profiles and converting logicals to scores:') }
   if(verbose){ message('  cnv + germline...') }
   out$cnv_germ <- (function(){
      ## Merge
      df_cnv <- mkValidDf('cnv')
      df_germ <- mkValidDf('germ')
      df_m <- merge(df_cnv, df_germ, by=SEL_COLS$all, all=T)


      ## Score
      df_m$full_gene_loss[df_m$full_gene_loss==1] <- SCORING$full_gene_loss
      df_m$loh[df_m$loh==1] <- SCORING$loh
      df_m$cn_break_in_gene[df_m$cn_break_in_gene==1] <- SCORING$cn_break_in_gene

      df_m$ref_loss[df_m$ref_loss==1] <- SCORING$germ.ref_loss
      df_m$alt_exists[df_m$alt_exists==1] <- SCORING$germ.alt_exists

      return(df_m)
   })()

   if(verbose){ message('  cnv + som...') }
   out$cnv_som <- (function(){
      ## Merge
      df_cnv <- mkValidDf('cnv')
      df_som <- mkValidDf('som')
      df_m <- merge(df_cnv, df_som, by=SEL_COLS$all, all=T)

      ## Score
      df_m$full_gene_loss[df_m$full_gene_loss==1] <- SCORING$full_gene_loss
      df_m$loh[df_m$loh==1] <- SCORING$loh
      df_m$cn_break_in_gene[df_m$cn_break_in_gene==1] <- SCORING$cn_break_in_gene

      df_m$ref_loss[df_m$ref_loss==1] <- SCORING$som.ref_loss
      df_m$alt_exists[df_m$alt_exists==1] <- SCORING$som.alt_exists

      return(df_m)
   })()

   if(verbose){ message('  germ + som...') }
   out$germ_som <- (function(){
      ## Prepate for merging
      df_cnv <- mkValidDf('cnv')[c(SEL_COLS$all,'cn_break_in_gene')]
      df_germ <- getGeneMaxEff(mkValidDf('germ'))
      df_som <- getGeneMaxEff(mkValidDf('som'))
      
      ## Score
      df_cnv$cn_break_in_gene[df_cnv$cn_break_in_gene==1] <- SCORING$cn_break_in_gene

      df_germ$ref_loss[df_germ$ref_loss==1] <- SCORING$germ.ref_loss
      df_germ$alt_exists[df_germ$alt_exists==1] <- SCORING$germ.alt_exists
   
      df_som$ref_loss[df_som$ref_loss==1] <- SCORING$som.ref_loss
      df_som$alt_exists[df_som$alt_exists==1] <- SCORING$som.alt_exists

      ## Merge
      df_m <- merge(df_cnv, df_germ, by=SEL_COLS$all, all=T)
      df_m <- merge(df_m, df_som, by=SEL_COLS$all, all=T)

      ## Convert suffix to prefix
      colnames(df_m) <- gsub("(.*).([xy])$", "\\2.\\1", colnames(df_m))

      ## Prefix with germ or som
      colnames(df_m) <- gsub('^x[.]','germ.',colnames(df_m))
      colnames(df_m) <- gsub('^y[.]','som.',colnames(df_m))

      return(df_m)
   })()

   ## In all dataframes, set NA to 0 to avoid downstream arithmetic errors
   out <- lapply(out, function(i){
      i[is.na(i)] <- 0
      return(i)
   })

   return(out)
}
