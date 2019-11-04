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
   #purity.path, 
   genes.bed.path, 
   ini.path
){
   
   #========= Inputs =========#
   options(stringsAsFactors=F)
   source(ini.path)
   
   # ## Testing
   # sample_name='CPCT02070055T' ## BRCA2 full gene loss
   # sample_name='CPCT02010419T' ## BRCA2 LOH + stop gained
   # sample_name='CPCT02010708T' ## BRCA1 LOH + frameshift
   # sample_name='CPCT02050231T' ## BRCA1 LOH + missense
   # sample_name='CPCT02010753T' ## BRCA2 LOH + 2 somatic frameshifts
   # in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/HMF_DR010_DR047/vcf_subset/',sample_name)
   # out.dir=paste0(in_dir,'/gene_statuses/')
   # ini.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts/pipeline/detGeneStatuses_ini.R'

   # sample_name='P11_A2960_A236'
   # in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/Rotterdam_Patient_Samples/gene_annotation/Rotterdam_Patient_Samples/',sample_name)
   # out.dir=paste0(in_dir,'/gene_statuses/')
   # ini.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/Rotterdam_Patient_Samples/scripts/annotate_genes/run_pipeline/detGeneStatuses_ini.R'
   
   # sample_name='SBT-3.1_organoid'
   # in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/Ovarian_Organoids_Chris/gene_ann//',sample_name)
   # out.dir=paste0(in_dir,'/gene_statuses/')
   
   # input_paths <- list(
   #    cnv = paste0(in_dir,'/',sample_name,'.purple.gene.cnv'),
   #    germ = paste0(in_dir,'/varsig/',sample_name,'_varsigs_germ.txt.gz'),
   #    som = paste0(in_dir,'/varsig/',sample_name,'_varsigs_som.txt.gz')
   # )
   # 
   # input <- list(
   #    cnv = read.delim(input_paths$cnv),
   #    germ = read.delim(input_paths$germ),
   #    som = read.delim(input_paths$som)
   # )
   
   ## Real
   input <- list(
      cnv = read.delim(cnv.path),
      germ = read.delim(germ.path),
      som = read.delim(som.path)
   )
   
   #genes.bed.path <- paste0(ROOT_DIR, '/data/gene_selection/genes.bed')
   #genes.bed.path <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/Rotterdam_Patient_Samples/scripts/annotate_genes/genes.bed'
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
   
   if(OPTIONS$verbose){ message('> som_som...') }
   l_gene_diplotypes$som_som <- mkGeneDiplotypesMutMut(
      mut_profile$som, mut_profile$som, 'som_som',
      verbose=OPTIONS$verbose
   )
   
   if(OPTIONS$verbose){ message('> germ_som...') }
   l_gene_diplotypes$germ_som <- mkGeneDiplotypesMutMut(
      mut_profile$germ, mut_profile$som, 'germ_som',
      verbose=OPTIONS$verbose
   )

   # subset(l_gene_diplotypes$cnv_som, hgnc_symbol %in% c('BRCA1','BRCA2'))
   # View(subset(l_gene_diplotypes$cnv_germ, hgnc_symbol %in% c('BRCA1','BRCA2')))
   # lapply(l_gene_diplotypes, function(i){ which(apply(i,1,anyNA)) })
   
   #========= Calculate hit scores and determine gene max effect =========#
   if(OPTIONS$verbose){ message('\n## Merging diplotype origins into one table...') }
   gene_diplotypes <- do.call(rbind, l_gene_diplotypes)
   
   if(OPTIONS$verbose){ message('\n## Calculating hit_scores...') }
   gene_diplotypes <- insColAfter(
      gene_diplotypes,
      calcHitScores(gene_diplotypes, DIPLOTYPE_ORIGIN_RANK),
      after='hgnc_symbol'
   )
   
   if(OPTIONS$verbose){ message('\n## Determining biallelic hit type...') }
   gene_diplotypes <- insColAfter(
      gene_diplotypes,
      detHitType(gene_diplotypes),
      after='diplotype_origin',
      colname='hit_type'
   )
   
   if(OPTIONS$verbose){ message('\n## Determining most pathogenic diplotype per gene...') }
   gene_diplotypes_max <- (function(){
      df <- getGeneDiplotypeMaxEff(gene_diplotypes, colname='hit_score')
      df <- df[order(df$hgnc_symbol),]
      #df$hit_score_boosted <- NULL
      return(df)
   })()
   
   # View(subset(gene_diplotypes_max, diplotype_origin %in% c('germ_som','som_som')))
   # subset(gene_diplotypes, is.na(hit_score_boosted))
   # subset(gene_diplotypes_max, hgnc_symbol %in% c('BRCA1','BRCA2'))
   # subset(gene_diplotypes, hgnc_symbol %in% c('BRCA1','BRCA2'))
   # table(gene_diplotypes_max$a1)
   
   #========= Determine reading frame/add indel info =========#
   gene_diplotypes_max <- cbind(
      gene_diplotypes_max,
      detReadingFrame(gene_diplotypes_max, mut_profile$germ, mut_profile$som)
   )
   
   #========= Export diplotype tables =========#
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











