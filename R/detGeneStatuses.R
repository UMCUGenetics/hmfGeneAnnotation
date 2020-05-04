#' Determine gene amplifications and biallelic losses
#'
#' @param out.dir Path to output dir
#' @param input.file.paths A list or vector supplying the path output files from the HMF pipeline.
#' @param sample.name Name of the sample
#' This list should be in the form c(germ_vcf='', som_vcf='', gene_cnv='', cnv='')
#' @param bed.file Path to the bed file containing the genes of interest. This bed file should also 
#' contain the column ensembl_gene_id
#' @param java.path Path the the java binary
#' @param snpsift.path Path to the SnpSift.jar
#' @param snpeff.path Path to snpEff.jar
#' @param do.snpeff.ann Annotate SNV/indels variant type with snpeff?
#' @param chrom.arm.split.method Can be 'hmf' or 'universal'. Refer to documentation for 
#' calcChromArmPloidies().
#' @param verbose Show progress messages?
#'
#' @return Writes a diplotypes table to the output dir
#' @export
#'
detGeneStatuses <- function(
   ## Input arguments
   out.dir, input.file.paths=c(germ_vcf=NA, som_vcf=NA, gene_cnv=NA, cnv=NA), sample.name,
   
   ## Package constants
   bed.file=BED_FILE, exons.bed.file=EXONS_BED_FILE,
   java.path=JAVA_PATH, snpsift.path=SNPSIFT_PATH, snpeff.path=SNPEFF_PATH,
   
   ## Defaults for HMF pipeline output
   chrom.arm.split.method='hmf',
   do.filter.vcf=T, do.snpeff.ann=F, do.det.reading.frame=T,
   ignore.chroms=NULL,
   
   verbose=T
){
   
####========= Inputs for debugging =========#

   # if(1==2){
   # 
   #    devtools::load_all('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/')
   # 
   #    ## HMF
   #    vcf_paths <- read.delim('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CUPs_classifier/metadata/vcf_paths.txt', stringsAsFactors=F)
   #    sample.name <- 'CPCT02010422T' ## BRCA2 LOH+som
   #    #sample.name <- 'CPCT02010543T' ## APC som stop gain (5) + som FS (5)
   # 
   #    out.parent.dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CUPs_classifier_data/gene_ann/'
   #    out.dir <- paste0(out.parent.dir,'/',sample.name,'/')
   #    dir.create(out.dir, recursive=T, showWarnings=F)
   # 
   #    input.file.paths <- unlist(vcf_paths[vcf_paths$sample==sample.name,-1])
   #    input.file.paths <- paste0('/Users/lnguyen/',input.file.paths)
   #    names(input.file.paths) <- c('germ_vcf','som_vcf','sv_vcf','gene_cnv','cnv')
   # 
   #    #------
   #    ## PCAWG
   #    vcf_paths <- read.delim('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/PCAWG_2020/scripts/annotate_genes/manifest_samples_with_all_data.txt', stringsAsFactors=F)
   #    sample.name <- '42465bbd-289b-4e96-98fe-76809c5e1520' ## BRCA1 LOH+som
   #    sample.name <- '005e85a3-3571-462d-8dc9-2babfc7ace21' ## Y chrom fail
   #    sample.name <- '07a7c634-bd9a-4fc2-b9fe-87b060ec3d1f' ## empty mut txt
   #    sample.name <- '1bea3a72-3b73-4072-a6bb-96a90119d3ac' ## calcChromArmPloidy fail
   #    sample.name <- '0009b464-b376-4fbc-8a56-da538269a02f' ## test
   # 
   #    out.parent.dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/PCAWG_2020/gene_ann/'
   #    out.dir <- paste0(out.parent.dir,'/',sample.name,'/')
   #    dir.create(out.dir, recursive=T, showWarnings=F)
   # 
   #    input.file.paths <- unlist(vcf_paths[vcf_paths$sample==sample.name,-1])
   #    input.file.paths <- paste0('/Users/lnguyen/',input.file.paths)
   #    names(input.file.paths) <- c('germ_vcf','som_vcf','cnv')
   # 
   #    chrom.arm.split.method='universal'
   #    do.filter.vcf=F
   #    do.snpeff.ann=T
   #    do.det.reading.frame=F
   #    ignore.chroms='Y'
   # 
   #    #------
   #    ## Common args
   #    bed.file=BED_FILE
   #    exons.bed.file=EXONS_BED_FILE
   #    java.path=JAVA_PATH
   #    snpsift.path='/Users/lnguyen/Documents/R_cache/snpEff/SnpSift.jar'
   #    snpeff.path='/Users/lnguyen/Documents/R_cache/snpEff/snpEff.jar'
   #    verbose=T
   # }
   
####========= Sanity checks =========#
   if(!dir.exists(out.dir)){
      stop('out.dir does not exist: ',out.dir)
   }
   
   if(length(names(input.file.paths))==0){
      stop('input.file.paths must be a named vector or list')
   }
   
   if(is.list(input.file.paths)){ ## Coerce list to vector
      input.file.paths <- structure(input.file.paths, names=names(input.file.paths))
   }
   
   # if(anyNA(input.file.paths[c('germ_vcf','som_vcf','cnv')])){
   #    stop('germ_vcf, som_vcf and cnv must be specified in input.file.paths')
   # }
   
   # if(any(input.file.paths[!is.na(input.file.paths)])){
   #    stop('One or more of the input files do not exist')
   # }

####========= Pre-process HMF pipeline output =========#
   preproc_dir <- paste0(out.dir,'/preproc/')
   dir.create(preproc_dir, recursive=T, showWarnings=F)
   
   preproc_files <- c()
   
   #--------- Cnv ---------#
   if(verbose){ message('\n## Calculating chrom arm ploidies...') }
   preproc_files['arm_ploidies'] <- paste0(preproc_dir,'/',sample.name,'.chrom_arm_ploidies.txt')
   
   if( !file.exists(preproc_files['arm_ploidies']) ){
      calcChromArmPloidies(
         input.file.paths['cnv'],
         out.file = preproc_files['arm_ploidies'],
         chrom.arm.split.method = chrom.arm.split.method,
         ignore.chroms = ignore.chroms,
         verbose = verbose
      )
   } else {
      if(verbose){ message('Skipping; output file exists: ',preproc_files['arm_ploidies']) }
   }
   
   #--------- Gene cnv ---------#
   preproc_files['gene_cnv'] <- paste0(preproc_dir,'/',sample.name,'.gene.cnv.txt')
   
   if( !file.exists(preproc_files['gene_cnv']) ){

      if( is.na(input.file.paths['gene_cnv']) ){
         if(verbose){ message('\n## Gene cnv file not specified. Generating from cnv file...') }
         cnvToGeneCnv(
            input.file.paths['cnv'],
            out.file = preproc_files['gene_cnv'],
            exons.bed.file = exons.bed.file,
            verbose = verbose
         )
      } else {
         if(verbose){ message('\n## Subsetting gene cnv...') }
         preProcessGeneCnv(
            gene.cnv.file = input.file.paths['gene_cnv'],
            out.file = preproc_files['gene_cnv'],
            bed.file = bed.file,
            verbose = verbose
         )
      }
   } else {
      if(verbose){ message('Skipping; output file exists: ',preproc_files['gene_cnv']) }
   }
   
   #--------- Somatic/germline vcfs ---------#
   bed_file <- read.delim(bed.file, stringsAsFactors=F, check.names=F)
   colnames(bed_file)[1] <- 'chrom'
   
   #input.file.paths['som_vcf'] <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/PCAWG_2020/scripts/merge_snv_indel_vcfs/test.vcf.gz'
   #counter <- 0
   for(mut_type in c('som','germ')){
      # counter <- counter + 1
      # if(counter==2){ break }
      # mut_type='som'
      
      MUT_VCF <- paste0(mut_type,'_vcf')
      MUT_TXT <- paste0(mut_type,'_txt')
      
      preproc_files[MUT_VCF] <- paste0(preproc_dir,'/',sample.name,'.',mut_type,'.vcf.gz')
      preproc_files[MUT_TXT] <- paste0(preproc_dir,'/',sample.name,'.',mut_type,'.txt.gz')
      
      if(verbose){ message('\n## Processing ',toupper(mut_type),' vcf...') }
      
      if(  !file.exists(preproc_files[MUT_VCF]) ){
         
         if(do.filter.vcf){
            if(verbose){ message('> Filtering vcf on variants in genes in bed file and PASS variants...') }
            filterVcf(
               vcf.file = input.file.paths[MUT_VCF],
               out.file = preproc_files[MUT_VCF],
               mode = mut_type,
               bed.file = bed.file,
               java.path = java.path,
               snpsift.path = snpsift.path
            )
         } else {
            message('> Skipping filtering vcf; copying vcf to out dir...')
            system(paste(
               'cp',
               input.file.paths[MUT_VCF],
               preproc_files[MUT_VCF]
            ))
         }
         
      } else {
         if(verbose){ message('> Skipping filtering vcf; output file exists: ',preproc_files[MUT_VCF]) }
      }
      
      ## Snpeff annotation
      if(do.snpeff.ann){
         
         ann_done <- paste0(preproc_dir,'/ann.',mut_type,'.done')
         
         if(!file.exists(ann_done)){
            if(verbose){ message('> Annotating variants...') }
            annotateVariantType(
               vcf.file = preproc_files[MUT_VCF], 
               out.file = preproc_files[MUT_VCF],
               genome = 'GRCh37.75',
               java.path = java.path, 
               snpeff.path = snpeff.path
            )
            
            file.create(ann_done)
         } else {
            if(verbose){ message("> Skipping annotating variants; done file exists") }
         }
      }
      
      ## Filtering
      if( !file.exists(preproc_files[MUT_TXT]) ){
         if(verbose){ message('> Extracting revelant fields in vcf to txt...') }
         extractVcfFields(
            vcf.file = preproc_files[MUT_VCF],
            out.file = preproc_files[MUT_TXT],
            java.path = java.path,
            snpsift.path = snpsift.path
         )
         
         txt <- read.delim(preproc_files[MUT_TXT], stringsAsFactors=F)
         
         if(!is.null(ignore.chroms)){
            if(verbose){ message('> Removing chromosomes: ',paste(ignore.chroms, collapse=',')) }
            txt <- txt[!(txt$chrom %in% ignore.chroms),]
         }
         
         if(verbose){ message('> Performing extra subsetting for ENSG ids present in bed file...') }
         txt <- txt[txt$ensembl_gene_id %in% bed_file$ensembl_gene_id,]
         
         if(verbose){ message('> Adding HGNC gene ids...') }
         #txt$hgnc_symbol <- ensgToHgncSymbol(txt$ensembl_gene_id)
         txt$hgnc_symbol <- bed_file[match(txt$ensembl_gene_id, bed_file$ensembl_gene_id),'hgnc_symbol']
         
         if(nrow(txt)!=0){
            if(verbose){ message('> Adding ClinVar annotations...') }
            txt$clinvar_sig <- getClinSig(txt, CLINVAR_PATH) ## seqminer::tabix.read returns error if dataframe is empty
         } else {
            if(verbose){ message('> No variants remain after filtering. Skipping adding ClinVar annotations...') }
            txt <- cbind(txt, data.frame(clinsig='')[0,,drop=F])
         }
         
         write.tsv(txt, preproc_files[MUT_TXT])
         rm(txt)
         
      } else {
         if(verbose){ message('> Skipping making txt; output file exists: ',preproc_files[MUT_TXT]) }
      }
   }
   
####========= Make preliminary output =========#
   #--------- Read files ---------#
   if(verbose){ message('\n## Loading input files...') }
   input_table_names <- c('arm_ploidies','gene_cnv','som_txt','germ_txt')
   
   input_tables <- lapply(input_table_names, function(i){ 
         message('Reading file: ', preproc_files[i] )
         read.delim(preproc_files[i], stringsAsFactors=F) 
      }
   )
   names(input_tables) <- input_table_names
   rm(input_table_names)
   
   # input_tables$arm_ploidies <- read.delim(preproc_files['arm_ploidies'], stringsAsFactors=F)
   # input_tables$gene_cnv     <- read.delim(preproc_files['gene_cnv'], stringsAsFactors=F)
   # input_tables$germ         <- read.delim(preproc_files['germ_txt'], stringsAsFactors=F)
   # input_tables$som          <- read.delim(preproc_files['som_txt'], stringsAsFactors=F)
   
   #--------- Make mut profiles ---------#
   mut_profile_dir <- paste0(out.dir,'/mut_profiles/')
   dir.create(mut_profile_dir, recursive=T, showWarnings=F)
   
   mut_profile_paths <- list(
      gene_cnv=paste0(mut_profile_dir,'/mut_profile_gene_cnv.txt.gz'),
      som=paste0(mut_profile_dir,'/mut_profile_som.txt.gz'),
      germ=paste0(mut_profile_dir,'/mut_profile_germ.txt.gz')
   )
   
   if( all(file.exists(unlist(mut_profile_paths))) ){
      mut_profile <- lapply(mut_profile_paths, read.delim, stringsAsFactors=F)
   } else {
      mut_profile <- list()
      
      if(verbose){ message('\n## Annotating gene CNV table...') }
      mut_profile$gene_cnv <- mkMutProfileGeneCnv(
         gene.cnv=input_tables$gene_cnv,
         arm.ploidies=input_tables$arm_ploidies,
         verbose=verbose
      )
      #subset(mut_profile$gene_cnv, hgnc_symbol=='BRCA2')
      
      if(verbose){ message('\n## Annotating germline and somatic mutations...') }
      if(verbose){ message('> som...') }
      mut_profile$som <- mkMutProfileSnvIndel(
         input_tables$som_txt,
         scoring=SCORING_MUT,
         keep.only.first.eff=T,
         filter.no.impact.variants=F,
         verbose=verbose
      )
      #subset(mut_profile$som, hgnc_symbol=='BRCA2')
      
      if(verbose){ message('> germ...') }
      mut_profile$germ <- mkMutProfileSnvIndel(
         input_tables$germ_txt,
         scoring=SCORING_MUT,
         keep.only.first.eff=T,
         filter.no.impact.variants=F,
         verbose=verbose
      )
      
      #--------- Write output ---------#
      if(verbose){ message('\n## Exporting mutation profiles...') }
      for(i in names(mut_profile)){
         write.tsv(mut_profile[[i]],mut_profile_paths[[i]])
      }
   }
   
   ## Force factors to characters
   mut_profile <- lapply(mut_profile, function(i){
      #i=mut_profile[[1]]
      as.data.frame(lapply(i, function(j){
         if(is.factor(j)){ as.character(j) }
         else(j)
      }), stringsAsFactors=F)
   })
   
####========= Make gene diplotypes =========#
   if(verbose){ message('\n## Making gene diplotypes tables...') }
   l_gene_diplotypes <- list()
   
   if(verbose){ message('> cnv_som...') }
   ## cnv_som placed first so that somatic mutations are prioritized over germline mutations
   l_gene_diplotypes$cnv_som <- mkGeneDiplotypesCnvMut(
      mut.profile.cnv = mut_profile$gene_cnv, 
      mut.profile.mut = mut_profile$som,
      mut.origin = 'som',
      verbose = verbose
   ) 
   
   if(verbose){ message('> cnv_germ...') }
   l_gene_diplotypes$cnv_germ <- mkGeneDiplotypesCnvMut(
      mut.profile.cnv = mut_profile$gene_cnv, 
      mut.profile.mut = mut_profile$germ, 
      mut.origin = 'germ',
      verbose=verbose
   )
   
   if(verbose){ message('> som_som...') }
   l_gene_diplotypes$som_som <- mkGeneDiplotypesMutMut(
      mut.profile.mut1 = mut_profile$som, 
      mut.profile.mut2 = mut_profile$som, 
      diplotype.origin = 'som_som',
      verbose = verbose
   )
   
   if(verbose){ message('> germ_som...') }
   l_gene_diplotypes$germ_som <- mkGeneDiplotypesMutMut(
      mut.profile.mut1 = mut_profile$germ, 
      mut.profile.mut2 = mut_profile$som, 
      diplotype.origin = 'germ_som',
      verbose=verbose
   )

   # subset(l_gene_diplotypes$cnv_som, hgnc_symbol %in% c('BRCA1','BRCA2'))
   # View(subset(l_gene_diplotypes$cnv_germ, hgnc_symbol %in% c('BRCA1','BRCA2')))
   # lapply(l_gene_diplotypes, function(i){ which(apply(i,1,anyNA)) })
   
####========= Calculate hit scores and determine gene max effect =========#
   if(verbose){ message('\n## Merging diplotype origins into one table...') }
   gene_diplotypes <- do.call(rbind, l_gene_diplotypes)
   rownames(gene_diplotypes) <- NULL
   
   ## Fix to remove duplicate rows for deep_deletion and trunc (due to being reported in both cnv_germ and cnv_som)
   gene_diplotypes <- unique(gene_diplotypes)
   
   if(verbose){ message('\n## Calculating hit_scores...') }
   gene_diplotypes <- insColAfter(
      gene_diplotypes,
      calcHitScores(gene_diplotypes),
      after='hgnc_symbol'
   )
   
   if(verbose){ message('\n## Determining biallelic hit type...') }
   gene_diplotypes <- (function(){
      df <- data.frame(
         hit_type=detBiallHitType(gene_diplotypes),
         biallelic_status=detBiallStatus(gene_diplotypes)
      )
      insColAfter(gene_diplotypes, df, after='diplotype_origin')
   })()
   
   if(verbose){ message('\n## Determining most pathogenic diplotype per gene...') }
   gene_diplotypes_max <- (function(){
      df <- getGeneDiplotypeMaxEff(gene_diplotypes, colname='hit_score')
      df <- df[order(df$hgnc_symbol),]
      #df$hit_score_boosted <- NULL
      return(df)
   })()
   
   if(do.det.reading.frame){
      if(verbose){ message('\n## Determining resulting reading frames due to multi-frameshifts...') }
      gene_diplotypes_max <- cbind(
         gene_diplotypes_max,
         detReadingFrame(gene_diplotypes_max, mut_profile$germ, mut_profile$som)
      )
   }

   # subset(gene_diplotypes_max, diplotype_origin %in% c('germ_som','som_som'))
   # subset(gene_diplotypes, is.na(hit_score_boosted))
   # subset(gene_diplotypes_max, hgnc_symbol %in% c('BRCA1','BRCA2'))
   # subset(gene_diplotypes, hgnc_symbol %in% c('BRCA1','BRCA2'))
   # table(gene_diplotypes_max$a1)
   
   if(verbose){ message('\n## Appending gene amplification info...') }
   gene_diplotypes_max <- (function(){
      df <- mut_profile$gene_cnv[
         match(gene_diplotypes_max$ensembl_gene_id, mut_profile$gene_cnv$ensembl_gene_id),
         c('amp_ratio_arm','amp_ratio_local','amp_ratio')
      ]
      insColAfter(gene_diplotypes_max, df, after='hgnc_symbol')
   })()
   
   #========= Export diplotype tables =========#
   if(verbose){ message('\n## Exporting gene diplotype tables...') }
   write.tsv(gene_diplotypes,paste0(out.dir,'/gene_diplotypes_full.txt.gz'))
   write.tsv(gene_diplotypes_max,paste0(out.dir,'/gene_diplotypes_with_amps.txt.gz'))
   
}
