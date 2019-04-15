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
   
   ## Testing
   # sample_name='CPCT02070023R_CPCT02070023TII'
   # in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_update/vcf_subset/',sample_name)
   # init.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/scripts/pipeline/detGeneStatuses_init.R'

   # sample_name='P10_A2988N_A228pT'
   # in_dir=paste0('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/Rotterdam_Patient_Samples/vcf_subset/',sample_name)
   # out.dir=paste0(in_dir,'/gene_statuses/')
   # init.path='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/Rotterdam_Patient_Samples/scripts/annotate_genes/detGeneStatuses_init.R'
   # 
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
   mut_profile <- list()
   
   if(OPTIONS$verbose){ message('\n## Making CNV mutation profile...') }
   mut_profile$cnv <- mkMutProfileCnv(
      input$cnv,
      rm.non.selected.genes=T,
      genes.bed=genes_bed,
      verbose=OPTIONS$verbose
   )
   
   if(OPTIONS$verbose){ message('\n## Making germline mutation profile...') }
   mut_profile$germ <- mkMutProfileSnvIndel(
      input$germ,
      gene.identifier=OPTIONS$gene.identifier,
      ignore.additional.evidence=OPTIONS$ignore.additional.evidence,
      tumor.purity=input$purity,
      mode='germline',
      rm.non.selected.genes=T,
      genes.bed=genes_bed
   )
   
   if(OPTIONS$verbose){ message('\n## Making somatic mutation profile...') }
   mut_profile$som <- mkMutProfileSnvIndel(
      input$som,
      gene.identifier=OPTIONS$gene.identifier,
      ignore.additional.evidence=OPTIONS$ignore.additional.evidence,
      tumor.purity=input$purity,
      mode='somatic',
      rm.non.selected.genes=T,
      genes.bed=genes_bed
   )
   
   #========= Calculate hit scores =========#
   if(OPTIONS$verbose){ message('\n## Making pairwise joins of cnv, germ, and som tables...') }
   biall_mut_profile <- mkBialleleMutProfile(mut_profile)
   
   if(OPTIONS$verbose){ message('\n## Determining hit_scores...') }

   if(OPTIONS$verbose){ message('  cnv + germ...') }
   biall_mut_profile$cnv_germ <- (function(){
      df <- biall_mut_profile$cnv_germ

      hit_scores <- do.call(rbind, Map(
         calcHitScore, mode='cnv_germ',

         full_gene_loss = df$full_gene_loss,
         loh = df$loh,

         germ.max_score = df$max_score,
         germ.alt_exists = df$alt_exists,
         germ.ref_loss = df$ref_loss,

         cn_break_in_gene = df$cn_break_in_gene,

         min.hit.score.filter = CUTOFFS$min.hit.score,

         USE.NAMES=F
      ))

      cbind(df, hit_scores)
   })()

   if(OPTIONS$verbose){ message('  cnv + som...') }
   biall_mut_profile$cnv_som <- (function(){
      df <- biall_mut_profile$cnv_som

      hit_scores <- do.call(rbind, Map(
         calcHitScore, mode='cnv_som',

         full_gene_loss = df$full_gene_loss,
         loh = df$loh,

         som.max_score = df$max_score,
         som.alt_exists = df$alt_exists,
         som.ref_loss = df$ref_loss,

         cn_break_in_gene = df$cn_break_in_gene,

         min.hit.score.filter = CUTOFFS$min.hit.score,

         USE.NAMES=F
      ))

      cbind(df, hit_scores)
   })()

   if(OPTIONS$verbose){ message('  germ + som...') }
   biall_mut_profile$germ_som <- (function(){
      df <- biall_mut_profile$germ_som

      hit_scores <- do.call(rbind, Map(
         calcHitScore, mode='germ_som',

         germ.max_score = df$germ.max_score,
         som.max_score = df$som.max_score,

         germ.alt_exists = df$germ.alt_exists,
         som.alt_exists = df$som.alt_exists,

         germ.ref_loss = df$germ.ref_loss,
         som.ref_loss = df$som.ref_loss,
         cn_break_in_gene = df$cn_break_in_gene,

         min.hit.score.filter = CUTOFFS$min.hit.score,

         USE.NAMES=F
      ))

      cbind(df, hit_scores)
   })()

   if(OPTIONS$verbose){ message('  Exporting tables...') }
   for(i in names(biall_mut_profile)){
      write.table(
         biall_mut_profile[[i]],
         gzfile(paste0(out.dir,'/biall_mut_profile_',i,'.txt.gz')),
         sep='\t', quote=F, row.names=F
      )
   }
   #subset(biall_mut_profile$cnv_germ, loh==5)

   # if(OPTIONS$verbose){ message('Exporting mutation profile...') }
   # write.tsv <- function(...){ write.table(..., sep='\t',quote=F,row.names=F) }
   # write.tsv(mut_profile$cnv, gzfile(mut_profile_paths$cnv))
   # write.tsv(mut_profile$germ, gzfile(mut_profile_paths$germ))
   # write.tsv(mut_profile$som, gzfile(mut_profile_paths$som))

   #========= Gene diplotypes =========#
   if(OPTIONS$verbose){ message('\n## Determining gene diplotypes...') }

   ## Placing cnv_som before cnv_germ prioritizes loh+som events being displayed after overall
   ## getGeneMaxEff step
   gene_diplotypes <- do.call(rbind, lapply(c('cnv_som','cnv_germ','germ_som'), function(i){
      if(OPTIONS$verbose){
         message(sprintf('  %s...',sub('_',' + ',i)))
      }
      getGeneDiplotypes(biall_mut_profile[[i]], i, simplify.snpeff.eff=T)
   }))
   #names(gene_diplotypes) <- c('cnv_som','cnv_germ','germ_som')

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











