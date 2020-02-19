#' Determines the reading frame in case of one or more frameshifts
#'
#' @param gene.diplotypes Table containing the mutation info of alleles 1 and 2
#' @param mut.profile.germ Germ mut profile table
#' @param mut.profile.som Som mut profile table
#'
#' @return A dataframe containing info of: the number of detected indels, the lengths of each
#' frameshift, the reading frame (sum of the frameshift lengths), and a logical indicating whether
#' the gene is out of frame
#' @export
#'
#' @examples
detReadingFrame <- function(gene.diplotypes, mut.profile.germ, mut.profile.som){
   # gene.diplotypes=gene_diplotypes_max
   # mut.profile.germ=mut_profile$germ
   # mut.profile.som=mut_profile$som
   
   has_fs <- with(gene.diplotypes,{ a1=='frameshift_variant' | a2=='frameshift_variant'})
   fs_genes <- gene.diplotypes[has_fs,'hgnc_symbol']
   
   ## Only look back into mut_profile if gene has frameshift
   fs_info <- do.call(rbind,lapply(fs_genes, function(i){
      #i='BRCA2'
      mut_profile_ss <- lapply(list(mut.profile.germ, mut.profile.som), function(j){
         j[with(j,{ hgnc_symbol==i & snpeff_eff=='frameshift_variant'}),]
      })
      # mut_profile_ss[[1]]$mut_origin <- 'germ'
      # mut_profile_ss[[2]]$mut_origin <- 'som'
      
      df <- rbind(mut_profile_ss[[1]], mut_profile_ss[[2]])
      
      ##
      fs_lengths <- nchar(df$alt) - nchar(df$ref)
      reading_frame <- sum(fs_lengths)
      if(reading_frame %% 3 == 0){
         out_of_frame <- 0
      } else {
         out_of_frame <- 1
      }
      
      ##
      fs_origins <- c(
         rep('germ',nrow(mut_profile_ss[[1]])),
         rep('som',nrow(mut_profile_ss[[2]]))
      )
      
      ## 
      fs_positions <- c(
         mut_profile_ss[[1]]$pos,
         mut_profile_ss[[2]]$pos
      )

      return(data.frame(
         n_fs=nrow(df),
         fs_lengths=paste(fs_lengths,collapse=','),
         fs_origins=paste(fs_origins, collapse=','),
         fs_positions=paste(fs_positions,collapse=','),
         reading_frame,
         out_of_frame
      ))
   }))
   
   ## Initialize output
   out <- data.frame(
      n_frameshifts=0,
      fs_lengths='',
      fs_origins='',
      fs_positions='',
      reading_frame=0,
      out_of_frame=0
   )
   out <- out[rep(1,nrow(gene.diplotypes)),]
   
   ## Fill in rows (genes) which have frameshift
   out[has_fs,] <- fs_info
   
   return(out)
}
