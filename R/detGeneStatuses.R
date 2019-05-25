#' Determine gene biallelic status
#'
#' @param out.dir Directory to write output
#' @param cnv.path Purple gene CNV table
#' @param germ.path Germline variant table
#' @param som.path Somatic variant table
#' @param purity.path Purple purity table
#' @param genes.bed.path Bed file of gene selection
#' @param ini.path An R script that provides the settings for this detGeneStatuses.
#'
#' @return Write the tables biall_mut_profile_*, gene_diplotypes, and gene_diplotypes_max to out.dir
#' @export
#'
detGeneStatuses <- function(
   out.dir, 
   cnv.path, 
   germ.path, 
   som.path, 
   purity.path, 
   genes.bed.path, 
   ini.path
){
   
   #========= Inputs =========#
   options(stringsAsFactors=F)
   
   ## Testing
   # sample_name='CPCT02070055T' ## BRCA2 full gene loss
   # sample_name='CPCT02010419T' ## BRCA2 LOH + stop gained
   # sample_name='CPCT02010708T' ## BRCA1 LOH + frameshift
   # sample_name='CPCT02050231T' ## BRCA1 LOH + missense
   # in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_DR010_DR047/vcf_subset/',sample_name)
   # out.dir=paste0(in_dir,'/gene_statuses/')
   # ini.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts/pipeline/detGeneStatuses_ini.R'
   # setwd("~/Documents")

   # input_paths <- list(
   #    cnv = paste0(in_dir,'/',sample_name,'.purple.gene.cnv'),
   #    germ = paste0(in_dir,'/varsig/',sample_name,'_varsigs_germ.txt.gz'),
   #    som = paste0(in_dir,'/varsig/',sample_name,'_varsigs_som.txt.gz'),
   #    purity = paste0(in_dir,'/',sample_name,'.purple.purity')
   # )
   # 
   # input <- list(
   #    cnv = read.delim(input_paths$cnv),
   #    germ = read.delim(input_paths$germ),
   #    som = read.delim(input_paths$som),
   #    purity = read.table(input_paths$purity, skip=1)[,1]
   # )
   
   ## Real
   source(ini.path)
   
   input <- list(
      cnv = read.delim(cnv.path),
      germ = read.delim(germ.path),
      som = read.delim(som.path),
      purity = read.table(purity.path, skip=1)[,1]
   )
   
   #genes.bed.path <- paste0(ROOT_DIR, '/data/gene_selection/genes.bed')
   genes_bed <- read.delim(genes.bed.path, check.names=F)
   
   #========= Make output dir =========#
   if(!dir.exists(out.dir)){ dir.create(out.dir) }
   
   #========= Add preliminary annotations =========#
   mut_profile_paths <- list(
      cnv=paste0(out.dir,'/mut_profile_cnv.txt.gz'),
      germ=paste0(out.dir,'/mut_profile_germ.txt.gz'),
      som=paste0(out.dir,'/mut_profile_som.txt.gz')
   )
   
   if(OPTIONS$overwrite.mut.profile==F & all(sapply(mut_profile_paths,file.exists))){
      mut_profile <- lapply(mut_profile_paths,read.delim)
   
   } else {
      mut_profile <- list()
      
      if(OPTIONS$verbose){ message('\n## Making mutation profile...') }
      if(OPTIONS$verbose){ message('> cnv...') }
      mut_profile$cnv <- mkMutProfileCnv(
         input$cnv,
         rm.non.selected.genes=T,
         genes.bed=genes_bed,
         verbose=OPTIONS$verbose
      )
      
      if(OPTIONS$verbose){ message('> som...') }
      mut_profile$som <- mkMutProfileSnvIndel(
         input$som,
         gene.identifier=OPTIONS$gene.identifier,
         rm.non.selected.genes=T,
         genes.bed=genes_bed,
         keep.only.first.eff=OPTIONS$keep.only.first.eff
      )
      
      if(OPTIONS$verbose){ message('> germ...') }
      mut_profile$germ <- mkMutProfileSnvIndel(
         input$germ,
         gene.identifier=OPTIONS$gene.identifier,
         rm.non.selected.genes=T,
         genes.bed=genes_bed,
         keep.only.first.eff=OPTIONS$keep.only.first.eff
      )
      
      if(OPTIONS$verbose){ message('\n## Exporting mutation profiles...') }
      for(i in names(mut_profile)){
         write.table(
            mut_profile[[i]],
            #gzfile(paste0(out.dir,'/mut_profile_',i,'.txt.gz')),
            gzfile(mut_profile_paths[[i]]),
            sep='\t', quote=F, row.names=F
         )
      }
   }
   
   #========= Make gene diplotypes =========#
   if(OPTIONS$verbose){ message('\n## Making gene diplotypes tables...') }
   l_gene_diplotypes <- list()
   
   if(OPTIONS$verbose){ message('> cnv_som...') }
   l_gene_diplotypes$cnv_som <- mkGeneDiplotypesCnvMut(
      mut_profile$cnv, mut_profile$som, 'som',
      verbose=OPTIONS$verbose
   ) ## Prioritize somatic mutations
   
   if(OPTIONS$verbose){ message('> cnv_germ...') }
   l_gene_diplotypes$cnv_germ <- mkGeneDiplotypesCnvMut(
      mut_profile$cnv, mut_profile$germ, 'germ',
      verbose=OPTIONS$verbose
   )
   
   if(OPTIONS$verbose){ message('> germ_som...') }
   l_gene_diplotypes$germ_som <- mkGeneDiplotypesGermSom(
      mut_profile$germ, mut_profile$som,
      verbose=OPTIONS$verbose
   )

   # subset(l_gene_diplotypes$cnv_som, hgnc_symbol %in% c('BRCA1','BRCA2'))
   # View(subset(l_gene_diplotypes$cnv_germ, hgnc_symbol %in% c('BRCA1','BRCA2')))
   # lapply(l_gene_diplotypes, function(i){ which(apply(i,1,anyNA)) })
   
   #========= Calculate hit scores and determine gene max effect =========#
   if(OPTIONS$verbose){ message('\n## Merging diplotype origins into one table...') }
   gene_diplotypes <- do.call(rbind, l_gene_diplotypes)
   
   if(OPTIONS$verbose){ message('\n## Calculating hit_scores...') }
   #gene_diplotypes$hit_score <- calcHitScores(gene_diplotypes)
   gene_diplotypes <- cbind(
      gene_diplotypes, 
      calcHitScores(gene_diplotypes, DIPLOTYPE_ORIGIN_RANK)
   )
   
   #head(gene_diplotypes)
   
   if(OPTIONS$verbose){ message('\n## Determining most pathogenic diplotype per gene...') }
   gene_diplotypes_max <- (function(){
      df <- getGeneDiplotypeMaxEff(gene_diplotypes, colname='hit_score')
      df <- df[order(df$hgnc_symbol),]
      #df$hit_score_boosted <- NULL
      return(df)
   })()
   
   # subset(gene_diplotypes_max, hgnc_symbol %in% c('BRCA1','BRCA2'))
   # subset(gene_diplotypes, hgnc_symbol %in% c('BRCA1','BRCA2'))
   # table(gene_diplotypes_max$a1)
   
   if(OPTIONS$verbose){ message('\n## Exporting gene diplotype tables...') }
   write.table(
      gene_diplotypes,
      gzfile(paste0(out.dir,'/gene_diplotypes.txt.gz')),
      sep='\t', quote=F, row.names=F
   )
   
   write.table(
      gene_diplotypes_max,
      gzfile(paste0(out.dir,'/gene_diplotypes_max.txt.gz')),
      sep='\t', quote=F, row.names=F
   )
   
}

# sample_name='CPCT02070023R_CPCT02070023TII'
# in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/',sample_name)
# 
# out.dir=paste0(in_dir,'/gene_statuses/')
# 
# input_paths <- list(
#    cnv = paste0(in_dir,'/',sample_name,'.purple.gene.cnv'),
#    germ = paste0(in_dir,'/varsig/',sample_name,'_varsigs_germ.txt.gz'),
#    som = paste0(in_dir,'/varsig/',sample_name,'_varsigs_som.txt.gz'),
#    purity = paste0(in_dir,'/',sample_name,'.purple.purity')
# )
# 
# genes.bed.path <- paste0(ROOT_DIR, '/data/gene_selection/genes.bed')
# 
# detGeneStatuses(
#    out.dir,
#    input_paths$cnv,
#    input_paths$germ,
#    input_paths$som,
#    input_paths$purity,
#    genes.bed.path
# )











