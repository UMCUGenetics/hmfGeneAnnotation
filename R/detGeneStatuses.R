#' Determine gene biallelic status
#'
#' @param out.dir Directory to write output
#' @param cnv.path Purple gene CNV table
#' @param germ.path Germline variant table
#' @param som.path Somatic variant table
#' @param purity.path Purple purity table
#' @param genes.bed.path Bed file of gene selection
#' @param init.path An R script that provides the settings for this detGeneStatuses.
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
   init.path
){
   
   #========= Inputs =========#
   options(stringsAsFactors=F)
   
   # # Testing
   # sample_name='CPCT02070023R_CPCT02070023TII'
   # sample_name='CPCT02040280R_CPCT02040280T'
   # sample_name='CPCT02020493R_CPCT02020493T'
   # sample_name='CPCT02010399R_CPCT02010399T'
   # sample_name='CPCT02010352R_CPCT02010352T'
   # sample_name='CPCT02070055R_CPCT02070055T'
   # in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/',sample_name)
   # out.dir=paste0(in_dir,'/gene_statuses/')
   # init.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts/pipeline/detGeneStatuses_init.R'

   # sample_name='P10_A2988_A228p'
   # in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/Rotterdam_Patient_Samples/vcf_subset/',sample_name)
   # out.dir=paste0(in_dir,'/gene_statuses/')
   # init.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/Rotterdam_Patient_Samples/scripts/annotate_genes/detGeneStatuses_init.R'
   
   # sample_name='PD11394'
   # in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/BRCA_EU/vcf_subset/',sample_name)
   # out.dir=paste0(in_dir,'/gene_statuses/')
   # init.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/BRCA_EU/scripts/annotate_genes/detGeneStatuses_init.R'
   
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
   # print(out.dir)
   # print(cnv.path)
   # print(germ.path)
   # print(som.path)
   # print(genes.bed.path)
   # print(init.path)
   
   source(init.path)
   
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
      if(OPTIONS$verbose){ message('  cnv...') }
      mut_profile$cnv <- mkMutProfileCnv(
         input$cnv,
         rm.non.selected.genes=T,
         genes.bed=genes_bed,
         verbose=OPTIONS$verbose
      )
      
      if(OPTIONS$verbose){ message('  germ...') }
      mut_profile$germ <- mkMutProfileSnvIndel(
         input$germ,
         gene.identifier=OPTIONS$gene.identifier,
         ignore.additional.evidence=OPTIONS$ignore.additional.evidence,
         tumor.purity=input$purity,
         mode='germline',
         rm.non.selected.genes=T,
         genes.bed=genes_bed
      )
      
      if(OPTIONS$verbose){ message('  som...') }
      mut_profile$som <- mkMutProfileSnvIndel(
         input$som,
         gene.identifier=OPTIONS$gene.identifier,
         ignore.additional.evidence=OPTIONS$ignore.additional.evidence,
         tumor.purity=input$purity,
         mode='somatic',
         rm.non.selected.genes=T,
         genes.bed=genes_bed
      )
      
      if(OPTIONS$verbose){ message('Exporting mutation profiles...') }
      for(i in names(mut_profile)){
         write.table(
            mut_profile[[i]],
            #gzfile(paste0(out.dir,'/mut_profile_',i,'.txt.gz')),
            gzfile(mut_profile_paths[[i]]),
            sep='\t', quote=F, row.names=F
         )
      }
   }
   
   #========= Calculate hit scores =========#
   biall_mut_profile_paths <- list(
      cnv_germ=paste0(out.dir,'/biall_mut_profile_cnv_germ.txt.gz'),
      cnv_som=paste0(out.dir,'/biall_mut_profile_cnv_som.txt.gz'),
      germ_som=paste0(out.dir,'/biall_mut_profile_germ_som.txt.gz')
   )
   
   biall_types <- c('cnv_som','cnv_germ','germ_som')
   
   if(OPTIONS$overwrite.biall.mut.profile==F & all(sapply(biall_mut_profile_paths,file.exists))){
      biall_mut_profile <- lapply(biall_mut_profile_paths,read.delim)
      names(biall_mut_profile) <- biall_types
   
   } else {
      if(OPTIONS$verbose){ message('\n## Making pairwise joins of cnv, germ, and som tables...') }
      biall_mut_profile <- mkBialleleMutProfile(mut_profile)
      
      if(OPTIONS$verbose){ message('\n## Determining hit_scores...') }
      biall_mut_profile <- lapply(biall_types, function(i){
         if(OPTIONS$verbose){ message(sprintf('  %s...', sub('_',' + ',i))) }
         calcHitScore(
            biall_mut_profile[[i]], mode=i,
            min.hit.score.filter=CUTOFFS$min.hit.score,
            ignore.additional.evidence=OPTIONS$ignore.additional.evidence
         )
      })
      names(biall_mut_profile) <- biall_types
      
      if(OPTIONS$verbose){ message('  Exporting tables...') }
      for(i in names(biall_mut_profile)){
         write.table(
            biall_mut_profile[[i]],
            gzfile(biall_mut_profile_paths[[i]]),
            sep='\t', quote=F, row.names=F
         )
      }
   }
   
   #========= Gene diplotypes =========#
   if(OPTIONS$verbose){ message('\n## Determining gene diplotypes...') }
   
   ## Placing cnv_som before cnv_germ prioritizes loh+som events being displayed after overall
   ## getGeneMaxEff step
   gene_diplotypes <- do.call(rbind, lapply(biall_types, function(i){
      if(OPTIONS$verbose){
         message(sprintf('  %s...',sub('_',' + ',i)))
      }
      getGeneDiplotypes(biall_mut_profile[[i]], i, simplify.snpeff.eff=T)
   }))

   if(OPTIONS$verbose){ message('\n## Determining most pathogenic diplotype per gene...') }
   gene_diplotypes_max <- (function(){
      df <- getGeneMaxEff(gene_diplotypes, colname='hit_score_boosted', show.n.max=T)
      df <- df[order(df$hgnc_symbol),]
      #df$hit_score_boosted <- NULL
      return(df)
   })()

   # subset(gene_diplotypes_max, hgnc_symbol %in% c('BRCA1','BRCA2'))
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











