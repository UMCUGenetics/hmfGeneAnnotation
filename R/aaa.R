#========= Load default paths on package load =========#
.onLoad <- function(libname, pkgname){
   
   ## Java
   assign('JAVA_PATH', system('which java',intern=T), envir=parent.env(environment()))
   
   ## Paths within the package
   pkg_constant_paths <- list(
      c('SNPSIFT_PATH','../dep/snpEff/SnpSift.jar'),
      c('SNPEFF_PATH','../dep/snpEff/snpEff.jar'), ## snpEff has requires databases that are downloaded to ./ the first time
      c('BED_FILE','misc/cosmic_cancer_gene_census_20200225.bed'),
      
      c('CENTROMERE_POSITIONS','db/centromere_positions_hg19.txt'),
      c('GENES_ENST2ENSG','db/human_genes_enst2ensg.txt.gz'),
      c('GENES_HGNC','db/hgnc_gene_names.txt.gz'),
      c('CLINVAR_PATH','db/clinvar.txt.bgz')
   )
   
   for(i in pkg_constant_paths){
      assign( i[1], system.file(i[2], package='hmfGeneAnnotation'), envir=parent.env(environment()) )
   }
   
   ## Load clinvar and snpeff scoring tables
   scoring_txt_files <- list(clinvar='scoring/clinvar_scoring.txt',snpeff='scoring/snpeff_scoring.txt')
   scoring_txt_files <- lapply(scoring_txt_files, function(i){ system.file(i,package='hmfGeneAnnotation') })
   
   SCORING_MUT <- lapply(scoring_txt_files, function(i){
      df <- read.delim(i, stringsAsFactors=F)
      structure(df$score, names=df$annotation)
   })
   
   assign('SCORING_MUT', SCORING_MUT, envir=parent.env(environment()))
}

ONE_ARMED_CHROMS <- c(13,14,15,21,22)
