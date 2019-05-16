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
   
   ## Select columns and reassign NA values
   SEL_COLS <- list(
      all = list(
         ensembl_gene_id='none', 
         hgnc_symbol='none'
      ),
      
      cnv = list(
         full_gene_loss=0,
         loh=0,
         cn_break_in_gene=0,
         amp=0
      ),
      
      germ = list(
         chrom=0,
         pos=0,
         hgvs_c='none',
         
         snpeff_eff='none',
         clinvar_sig='none',
         enigma_sig='none',
         max_score=0,
         max_score_origin='none',
         
         is_hotspot_mut=0,
         adj_tumor_ad_ref=0,
         adj_tumor_ad_alt=0,
         alt_exists=0,
         ad_diff_score=0,
         ref_loss=0
      )
   )
   SEL_COLS$som <- SEL_COLS$germ
   
   ## NA lookup table
   SEL_COLS_uniq <- unlist(unname(SEL_COLS), recursive=F)
   SEL_COLS_uniq <- SEL_COLS_uniq[!duplicated(names(SEL_COLS_uniq))]
   
   ## Apply column selection
   mut_profile_ss <- lapply(names(mut_profile), function(i){
      col_names <- c(names(SEL_COLS[['all']]), names(SEL_COLS[[i]]))
      df <- mut_profile[[i]][col_names]
      return(df)
   })
   names(mut_profile_ss) <- names(mut_profile)
   
   ## Init output
   out <- list()

   ## NOTE: slightly different for each combination
   if(verbose){ message('Merging mutation profiles and converting logicals to scores:') }
   
   applyScoringCnvMut <- function(mode){
      #mode='som'
      ## Merge
      df <- merge(
         mut_profile_ss$cnv, 
         mut_profile_ss[[mode]], 
         by=names(SEL_COLS$all), all=T
      )
      
      ## Fill NAs
      for(i in colnames(df)){
         df[,i][ is.na(df[,i]) ] <- SEL_COLS_uniq[[i]]
      }
      
      df_split <- list(
         full_gene_loss=df[df$full_gene_loss==1,],
         not_full_gene_loss=df[df$full_gene_loss==0,]
      )
      
      ## Full gene losses need to be dealt with separately
      if(nrow(df_split$full_gene_loss)!=0){
         df_split$full_gene_loss <- unique(df_split$full_gene_loss)
         
         df_split$full_gene_loss <- within(df_split$full_gene_loss,{
            full_gene_loss <- SCORING$full_gene_loss
            loh <- 0
            max_score <- 5
            max_score_origin <- 'full_gene_loss'
         })
      }
      
      ## CNV + mut
      df_split$not_full_gene_loss <- within(df_split$not_full_gene_loss,{
         loh[loh==1] <- SCORING$loh
         cn_break_in_gene[cn_break_in_gene==1] <- SCORING$cn_break_in_gene
         is_hotspot_mut[is_hotspot_mut==1] <- SCORING[[paste0(mode,'.is_hotspot_mut')]]
         ref_loss[ref_loss==1] <- SCORING[[paste0(mode,'.ref_loss')]]
         alt_exists[alt_exists==1] <- SCORING[[paste0(mode,'.alt_exists')]]
      })
      
      df <- do.call(rbind, df_split)
      rownames(df) <- NULL
      return(df)
   }
   
   if(verbose){ message('  cnv + germ...') }
   out$cnv_germ <- applyScoringCnvMut('germ')
   
   if(verbose){ message('  cnv + som...') }
   out$cnv_som <- applyScoringCnvMut('som')
   
   if(verbose){ message('  germ + som...') }
   out$germ_som <- (function(){
      ## Prepate for merging
      df_cnv <- mut_profile_ss$cnv[c(names(SEL_COLS$all),'cn_break_in_gene')]
      df_germ <- getGeneMaxEff(mut_profile_ss$germ)
      df_som <- getGeneMaxEff(mut_profile_ss$som)
      
      ## Score
      df_cnv$cn_break_in_gene[df_cnv$cn_break_in_gene==1] <- SCORING$cn_break_in_gene
      
      applyScoringGermSomMut <- function(df, mode){
         #df[is.na(df)] <- 0
         
         within(df,{
            is_hotspot_mut[is_hotspot_mut==1] <- SCORING[[paste0(mode,'.is_hotspot_mut')]]
            ref_loss[ref_loss==1] <- SCORING[[paste0(mode,'.ref_loss')]]
            alt_exists[alt_exists==1] <- SCORING[[paste0(mode,'.alt_exists')]]
         })
      }
      
      df_germ <- applyScoringGermSomMut(df_germ,'germ')
      df_som <- applyScoringGermSomMut(df_som,'som')

      ## Merge
      df <- merge(df_cnv, df_germ, by=names(SEL_COLS$all), all=T)
      df <- merge(df, df_som, by=names(SEL_COLS$all), all=T)

      ## Convert suffix to prefix
      colnames(df) <- gsub("(.*).([xy])$", "\\2.\\1", colnames(df))

      ## Prefix with germ or som
      colnames(df) <- gsub('^x[.]','germ.',colnames(df))
      colnames(df) <- gsub('^y[.]','som.',colnames(df))
      
      ## Fill NAs
      root_colnames <- gsub('germ[.]|som[.]','',colnames(df))
      for(i in 1:ncol(df)){
         df[,i][ is.na(df[,i]) ] <- SEL_COLS_uniq[[root_colnames[i]]]
      }
      
      return(df)
   })()

   return(out)
}
