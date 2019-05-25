####################################################################################################
BIALLELE_COLS <- list(
   
   common=list(
      #sample='none',
      ensembl_gene_id='none',
      hgnc_symbol='none'
   ),
   
   mut_profile_origin=list(
      diplotype_origin='none'#,
      # a1.origin='none',
      # a2.origin='none'
   ),
   
   allele1=list(
      a1.chrom='none',
      a1.pos=0,
      a1.hgvs_c='none',
      a1='none',
      a1.max_score=0,
      a1.max_score_origin='none'
   ),
   
   allele2=list(
      a2.chrom='none',
      a2.pos=0,
      a2.hgvs_c='none',
      a2='none',
      a2.max_score=0,
      a2.max_score_origin='none'
   )
)

getBialleleCols <- function(groups=NULL, as=NULL){
   l <- 
      if(is.null(groups)){ unname(BIALLELE_COLS) }
      else { unname(BIALLELE_COLS[groups]) }
   
   l <- unlist(l, recursive=F)
   
   ## Return
   if(as=='names'){ names(l) } 
   else if(as=='data.frame'){ as.data.frame(l) } 
   else { l }
}

####################################################################################################
#' Merge cnv and germ/som mut profiles into a table the describes the events of allele 1 and 2
#'
#' @param mut.profile.cnv CNV mut profile table
#' @param mut.profile.mut germ or som mut profile table
#' @param mut.origin Indicate whether allele 2 originates from 'germ' or 'som'. Used to create
#' the diplotype_origin column
#' @param verbose Print progress messages?
#'
#' @return A table the describes the events of allele 1 and 2
#' @export
mkGeneDiplotypesCnvMut <- function(mut.profile.cnv, mut.profile.mut, mut.origin, verbose=T){
   # mut.profile.cnv=mut_profile$cnv
   # mut.profile.mut=mut_profile$som
   # mut.origin='som'
   # mut.profile.mut=mut_profile$germ
   # mut.origin='germ'
   
   if(!(mut.origin %in% c('germ','som'))){ stop("mut.origin must be 'germ' or 'som'") }
   
   union_genes <- unique(c( mut.profile.cnv$ensembl_gene_id, mut.profile.mut$ensembl_gene_id ))
   
   if(verbose){
      pb <- txtProgressBar(min=0, max=length(union_genes), initial=0, style=3, width=100)
      counter <- 0
   }
   
   l <- lapply(union_genes, function(i){
      if(verbose){
         counter <<- counter+1
         setTxtProgressBar(pb,counter)
      }

      #i="ENSG00000225178"
      #i='ENSG00000012048' ## BRCA1
      #i='ENSG00000139618' ## BRCA2
      #i='ENSG00000141510' ## TP53
      #i='ENSG00000067606' ## PRKCZ
      #i='ENSG00000175279' ## CENPS
      #i='test_null'
      
      df_cnv <- mut.profile.cnv[mut.profile.cnv$ensembl_gene_id==i,]
      df_mut <- mut.profile.mut[mut.profile.mut$ensembl_gene_id==i,]
      
      if(nrow(df_cnv)==0 & nrow(df_mut)==0){
         return( getBialleleCols(as='data.frame') )
      }
      
      #========= Full gene loss / truncation cases =========#
      if( nrow(df_cnv)!=0 && df_cnv$cnv_eff %in% c('full_gene_loss','trunc') ){
         ## df_cnv should only ever have one row
         out <- getBialleleCols(as='data.frame')
         
         ## Scores
         out$a1.max_score <- df_cnv$a1.score
         out$a2.max_score <- df_cnv$a2.score
         
         ## Variant effect
         out$a1 <- out$a2 <- df_cnv$cnv_eff
         
         out <- within(out,{
            a1.max_score_origin <- a2.max_score_origin <- 'cnv'
            #a1.origin <- a2.origin <- 'cnv'
            diplotype_origin <- 'cnv_cnv'
         })
         
         ## Common
         common_cols <- getBialleleCols('common',as='names')
         out[common_cols] <- df_cnv[common_cols]
         
         return(out)
      }
      
      #========= Main =========#
      common_cols <- getBialleleCols('common',as='names')
      
      #--------- CNV (left hand side) ---------#
      out_cnv <- getBialleleCols(c('common','mut_profile_origin','coords','allele1'),as='data.frame')
      
      if(nrow(df_cnv)!=0){
         out_cnv$a1.chrom <- df_cnv$chrom
         out_cnv$a1.max_score <- df_cnv$a1.score
         out_cnv$a1 <- df_cnv$cnv_eff
         out_cnv$a1.max_score_origin <- 'cnv'
      }
      
      ## fill origin cols
      out_cnv$diplotype_origin <- paste0('cnv_',mut.origin)
      # out_cnv <- within(out_cnv, {
      #    diplotype_origin <- paste0('cnv_',mut.origin)
      #    a1.origin <- 'cnv'
      #    a2.origin <- mut.origin
      # })
      
      #--------- germ or som (right hand side) ---------#
      out_mut <- getBialleleCols(c('common','coords','allele2'),as='data.frame')
      
      if(nrow(df_mut)!=0){
         ## Initialize empty df
         out_mut <- getBialleleCols(c('common','coords','allele2'),as='data.frame')
         out_mut <- out_mut[rep(1,nrow(df_mut)),]
         out_mut[common_cols] <- df_mut[common_cols]
         
         ## Fill metadata
         meta_cols <- getBialleleCols(c('common','coords'),as='names')
         out_mut[meta_cols] <- df_mut[meta_cols]
         
         ## Fill allele data
         allele_cols_out <- getBialleleCols('allele2',as='names')
         
         allele_cols_source <- allele_cols_out
         allele_cols_source[allele_cols_source=='a2'] <- 'snpeff_eff'
         allele_cols_source <- gsub('a2[.]','',allele_cols_source)
         
         out_mut[allele_cols_out] <- df_mut[allele_cols_source]
      }
      
      ## Add gene names
      out_cnv[common_cols] <- if(nrow(df_cnv)!=0){
         df_cnv[common_cols]
      } else {
         df_mut[1,common_cols]
      }
      
      #--------- Join ---------#
      #out <- merge(out_cnv, out_mut, all=T)
      out <- cbind(out_cnv, out_mut[!(colnames(out_mut) %in% common_cols)])
      
      ## column check
      #ncol(out)
      #length(getBialleleCols(as='names'))
      
      return(out)
   })
   
   return( do.call(rbind,l) )
}

# df_cnv_mut <- mergeCnvMut(mut_profile$cnv,mut_profile$germ,'germ')
# View(df_cnv_mut)
# table(df_cnv_mut$a1)

####################################################################################################
#' Merge germ and som mut profiles into a table the describes the events of allele 1 and 2
#'
#' @description Per gene, only the most pathogenic germ/som variant pairs are merged. This
#' prevents the table from growing exponentially.
#'
#' @param mut.profile.cnv CNV mut profile table
#' @param mut.profile.mut germ or som mut profile table
#' @param mut.origin Indicate whether allele 2 originates from 'germ' or 'som'. Used to create
#' the diplotype_origin column
#' @param verbose Print progress messages?
#'
#' @return A table the describes the events of allele 1 and 2
#' @export
mkGeneDiplotypesGermSom <- function(mut.profile.germ, mut.profile.som, verbose=T){
   #mut.profile.germ=mut_profile$germ
   #mut.profile.som=mut_profile$som
   
   union_genes <- unique(c( mut.profile.germ$ensembl_gene_id, mut.profile.som$ensembl_gene_id ))
   
   if(verbose){
      pb <- txtProgressBar(min=0, max=length(union_genes), initial=0, style=3, width=100)
      counter <- 0
   }
   
   l <- lapply(union_genes, function(i){
      if(verbose){
         counter <<- counter+1
         setTxtProgressBar(pb,counter)
      }
      
      #print(i)
      
      #i='ENSG00000114315'
      #i='ENSG00000012048' ## BRCA1
      #i='ENSG00000139618' ## BRCA2
      #i='ENSG00000141510' ## TP53
      #i='ENSG00000067606' ## PRKCZ
      #i='test_null'
      
      df_germ <- mut.profile.germ[mut.profile.germ$ensembl_gene_id==i,]
      df_som <- mut.profile.som[mut.profile.som$ensembl_gene_id==i,]
      
      #========= Col names =========#
      common_cols <- getBialleleCols('common',as='names')
      
      sel_cols <- list(
         common=getBialleleCols('common',as='names'),
         a1_out=getBialleleCols('allele1',as='names'),
         a2_out=getBialleleCols('allele2',as='names')
      )
      
      sel_cols$a1_source <- (function(){
         v <- sel_cols$a1_out
         v[v=='a1'] <- 'snpeff_eff'
         v <- gsub('a1[.]','',v)
         return(v)
      })()
      
      sel_cols$a2_source <- (function(){
         v <- sel_cols$a2_out
         v[v=='a2'] <- 'snpeff_eff'
         v <- gsub('a2[.]','',v)
         return(v)
      })()
      
      #========= Germ =========#
      out_germ <- getBialleleCols(c('common','mut_profile_origin','coords','allele1'),as='data.frame')
      
      if(nrow(df_germ)!=0){
         df_germ <- df_germ[which.max(df_germ$max_score),]
         # out_germ[common_cols] <- df_germ[common_cols]
         out_germ[sel_cols$a1_out] <- df_germ[sel_cols$a1_source]
      }
      
      ## fill origin cols
      out_germ$diplotype_origin <- 'germ_som'
      # out_germ <- within(out_germ, {
      #    diplotype_origin <- 'germ_som'
      #    a1.origin <- 'germ'
      #    a2.origin <- 'som'
      # })
      
      #========= Som =========#
      out_som <- getBialleleCols(c('common','coords','allele2'),as='data.frame')
      
      if(nrow(df_som)!=0){
         df_som <- df_som[which.max(df_som$max_score),]
         # out_som[common_cols] <- df_som[common_cols]
         out_som[sel_cols$a2_out] <- df_som[sel_cols$a2_source]
      }
      
      #========= Add gene names =========#
      out_germ[common_cols] <- 
         if(nrow(df_germ)!=0){
            df_germ[common_cols]
         } else {
            df_som[common_cols]
         }
      
      ##
      out <- cbind(out_germ, out_som[!(colnames(out_som) %in% common_cols)])
      
      return(out)
      #ncol(out)
      #length(getBialleleCols(as='names'))
   })
   
   return( do.call(rbind,l) )
}
#df_germ_som <- mergeGermSom(mut_profile$germ, mut_profile$som)
#View(df_germ_som)

#df <- rbind(df_cnv_mut, df_germ_som)
