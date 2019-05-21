#' Per gene, get row(s) where score is max
#'
#' @param df A dataframe containing ensembl_gene_id
#' @param mode Can be: 'cnv_germ', 'cnv_som', 'germ_som'.
#' @param simplify.snpeff.eff Simplify snpeff_eff (e.g. stop_gained, stop_loss; becomes nonsense)
#'
#' @return A subsetting dataframe
#' @export
#'
getGeneDiplotypes <- function(df, mode, simplify.snpeff.eff=T){
   
   #df=biall_mut_profile$cnv_som
   #df=biall_mut_profile$germ_som
   
   SEL_COLS <- list(
      common = list(
         ensembl_gene_id='none', 
         hgnc_symbol='none',
         hit_score=0.0,
         hit_score_boosted=0,
         hit_type='none',
         def_type='none',
         is_def=0,
         cn_break_in_gene=0
      ),
      
      allele = list(
         chrom=0,
         pos=0,
         hgvs_c='none',
         
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
   SEL_COLS_uniq <- unlist(unname(SEL_COLS), recursive=F)

   modes <- c('cnv_germ','cnv_som','germ_som')
   if(!(mode %in% modes)){
      stop(paste0("Please specify a valid mode: 'cnv_germ', 'cnv_som', 'germ_som'"))
   }

   simplifySnpeffEff <- function(v){
      v_new <- SNPEFF_SIMPLE_ANN_LOOKUP[match(v, SNPEFF_SIMPLE_ANN_LOOKUP$ann), 'ann_s2']
      v_new[is.na(v_new)] <- 'none'
      return(v_new)
   }

   if(mode=='cnv_germ' | mode=='cnv_som'){
      origin_string <- if(mode=='cnv_germ'){ 'germ' } else if(mode=='cnv_som'){ 'som' }
      
      diplotypes <- as.data.frame(do.call(rbind, Map(function(full_gene_loss, loh, snpeff_eff){
         
         if(full_gene_loss==SCORING$full_gene_loss){
            a1 <- 'full_gene_loss'
            a2 <- 'full_gene_loss'
            a1.origin <- a2.origin <- 'cnv'
         } else if(loh==SCORING$loh) {
            a1 <- 'loh'
            a1.origin <- 'cnv'
            a2 <- snpeff_eff
            a2.origin <- origin_string
         } else {
            a1 <- 'none'
            a1.origin <- 'cnv'
            a2 <- snpeff_eff
            a2.origin <- origin_string
         }
         
         return(c(
            a1,a2,a1.origin,a2.origin
         ))
      }, df$full_gene_loss, df$loh, df$snpeff_eff)))
      colnames(diplotypes) <- c('a1','a2','a1.origin','a2.origin')
      
      #--------- Allele 1 ---------#
      ## Initiate empty dataframe with same nrows as input df
      out_a1 <- (function(){
         col_names <- c(paste0('a1.', names(SEL_COLS$allele)))

         out_a1 <- data.frame(matrix(nrow=nrow(df),ncol=length(col_names)))
         colnames(out_a1) <- col_names

         out_a1 <- cbind(out_a1, diplotypes[,c('a1','a1.origin')])

         return(out_a1)
      })()
      
      ## Add chrom info
      out_a1$a1.chrom <- df$chrom
      
      ## Assign score for loh/full_gene_loss
      out_a1$a1.max_score <- ifelse(out_a1$a1 %in% c('full_gene_loss','loh'),5,0)

      #--------- Allele 2 ---------#
      out_a2 <- df[,names(SEL_COLS$allele)]
      colnames(out_a2) <- paste0('a2.',colnames(out_a2))
      
      out_a2 <- cbind(out_a2, diplotypes[,c('a2','a2.origin')])
      
      out_a2$a2 <- 
         if(simplify.snpeff.eff){
            simplifySnpeffEff(diplotypes$a2)
         } else {
            diplotypes$a2
         }
      
      ## Export
      out <- cbind(df[,names(SEL_COLS$common)], diplotype_origin = mode, out_a1, out_a2)
   }

   if(mode=='germ_som'){
      out <- cbind(
         df[,names(SEL_COLS$common)],
         diplotype_origin = mode,
         
         df[,paste0('germ.',names(SEL_COLS$allele))],
         a1=df$germ.snpeff_eff,
         a1.origin='germ',
         
         df[,paste0('som.',names(SEL_COLS$allele))],
         a2=df$som.snpeff_eff,
         a2.origin='som'
      )

      colnames(out) <- gsub('germ','a1',colnames(out))
      colnames(out) <- gsub('som','a2',colnames(out))

      if(simplify.snpeff.eff){
         out$a1 <- simplifySnpeffEff(out$a1)
         out$a2 <- simplifySnpeffEff(out$a2)
      }
   }
   
   ## Fill in NAs
   root_colnames <- gsub('a1[.]|a2[.]','',colnames(out))
   fill_cols <- which(root_colnames %in% names(SEL_COLS_uniq))
   for(i in fill_cols){
      out[,i][ is.na(out[,i]) ] <- SEL_COLS_uniq[[root_colnames[i]]]
   }

   return(out)
}
