#' Create annotations to the HMF CNV table
#'
#' @description Primarily used to determine full gene loss and LOH
#'
#' @param df.cnv HMF *.gene.cnv table
#' @param min.cn.diff.in.gene Placeholder
#' @param max.max.copy.number Placeholder
#' @param min.min.copy.number Placeholder
#' @param max.min.minor.allele.ploidy Placeholder
#' @param rm.non.selected.genes Remove genes not specified by the user in genes.bed
#' @param genes.bed A bed file which contains ENSEMBL gene ids
#'
#' @return The original input dataframe with annotation columns
#' @export
#'

# sample_name='CPCT02070023R_CPCT02070023TII'
# out_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/',sample_name)
#
# input_paths <- list(
#    cnv = paste0(out_dir,'/',sample_name,'.purple.gene.cnv'),
#    germ = paste0(out_dir,'/varsig/',sample_name,'_varsigs_germ.txt.gz'),
#    som = paste0(out_dir,'/varsig/',sample_name,'_varsigs_som.txt.gz'),
#    purity = paste0(out_dir,'/',sample_name,'.purple.purity')
# )
#
# input <- list(
#    cnv = read.delim(input_paths$cnv),
#    germ = read.delim(input_paths$germ),
#    som = read.delim(input_paths$som),
#    purity = read.table(input_paths$purity, skip=1)[,1]
# )
# df.cnv <- input$cnv

mkMutProfileCnv <- function(
   df.cnv,
   min.cn.diff.in.gene = CUTOFFS$min.cn.diff.in.gene,
   max.max.copy.number = CUTOFFS$max.max.copy.number,
   min.min.copy.number = CUTOFFS$min.min.copy.number,
   max.min.minor.allele.ploidy = CUTOFFS$max.min.minor.allele.ploidy,
   rm.non.selected.genes=T, genes.bed=NULL,
   verbose=T
){
   if(nrow(df.cnv) == 0){
      # col_names <- c('chrom','start','end','hgnc_symbol','ensembl_gene_id',
      #   'min_copy_number','max_copy_number','min_minor_allele_ploidy',
      #   'full_gene_loss','loh','cn_diff_in_gene','cn_break_in_gene','amp'
      #   )
      #
      # stop(return(
      #    setNames(data.frame(matrix(ncol=length(col_names), nrow=0)), col_names)
      # ))
      if(verbose){ warning('No rows are present in the gene cnv table. Returning NA...') }
      stop(return(NA))
   }

   if(verbose){ message('Annotating gene cnv table...') }
   out <- do.call(rbind, lapply(1:nrow(df.cnv), function(i){
      #i=1
      row <- df.cnv[i,]

      ## Initiate scoring vector
      full_gene_loss=0
      loh=0
      amp=0
      cn_diff_in_gene=0
      cn_break_in_gene=0

      ## 'full gene loss' and 'loh' scores
      if(row$max_copy_number < max.max.copy.number){ full_gene_loss <- 1 }
      if(row$min_minor_allele_ploidy < max.min.minor.allele.ploidy){ loh <- 1 }

      cn_diff_in_gene <- abs(row$max_copy_number - row$min_copy_number) ## max should be > min, but do abs() just in case
      if(cn_diff_in_gene >= min.cn.diff.in.gene){ cn_break_in_gene <- 1 }

      if(row$min_copy_number >= min.min.copy.number){ amp <- 1 }

      cbind(
         row,
         full_gene_loss, loh, cn_diff_in_gene, cn_break_in_gene, amp
      )
   }))

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
