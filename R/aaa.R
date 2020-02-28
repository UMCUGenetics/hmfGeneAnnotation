#========= Load default paths on package load =========#
.onLoad <- function(libname, pkgname){
   
   ##
   assign(
      'JAVA_PATH', system('which java',intern=T), 
      envir=parent.env(environment())
   )
   
   ##
   pkg_constant_paths <- list(
      c('SNPSIFT_PATH','dep/SnpSift.jar'),
      c('BED_FILE','misc/cosmic_cancer_gene_census_20200225.bed'),
      c('GENES_ENST2ENSG','db/human_genes_enst2ensg.txt.gz'),
      c('GENES_HGNC','db/hgnc_gene_names.txt.gz'),
      c('CLINVAR_PATH','db/clinvar.txt.bgz')
   )
   
   for(i in pkg_constant_paths){
      assign( i[1], system.file(i[2], package='hmfGeneAnnotation'), envir=parent.env(environment()) )
   }
   
}

#========= Constants =========#
SCORING_MUT <- list()

SCORING_MUT$clinvar <- c(
   'Pathogenic'=5,
   'Pathogenic/Likely_pathogenic'=5,
   'Likely_pathogenic'=4,
   'Uncertain_significance'=3,
   'Conflicting_interpretations_of_pathogenicity'=3,
   'Likely_benign'=2,
   'Benign/Likely_benign'=2,
   'Benign'=1
)

SCORING_MUT$snpeff <- c(
   'frameshift_variant'= 5,
   'rare_amino_acid_variant'=5,
   'transcript_ablation'=5,
   'chromosome_number_variation'=4,
   'exon_loss_variant'=4,
   'start_lost'=4,
   'stop_gained'=4,
   'stop_lost'=4,
   'splice_acceptor_variant'=3,
   'splice_donor_variant'=3,
   'missense_variant'=3,
   'coding_sequence_variant'=3,
   'disruptive_inframe_deletion'=3,
   'disruptive_inframe_insertion'=3,
   'splice_region_variant'=2,
   '3_prime_UTR_truncation'=2,
   '5_prime_UTR_truncation'=2,
   'conservative_inframe_insertion'=2,
   'conservative_inframe_deletion'=2,
   'inframe_insertion'=2,
   'inframe_deletion'=2,
   'regulatory_region_ablation'=2,
   'TFBS_ablation'=2,
   '5_prime_UTR_premature_start_codon_gain_variant'=1,
   '3_prime_UTR_variant'=1,
   '5_prime_UTR_variant'=1,
   'initiator_codon_variant'=1,
   'TF_binding_site_variant'=1
   ## All other variant types have a score of 0
)

# names(SCORING_MUT$snpeff)

# SCORING_TABLES_DIR <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/scoring'
# 
# mkScoringVector <- function(scoring.table.path){
#    df <- read.delim(scoring.table.path,stringsAsFactors=F)
#    v <- df$score
#    names(v) <- df$ann
#    return(v)
# }
# SCORING_MUT <- list(
#    ## databases
#    snpeff = mkScoringVector(paste0(SCORING_TABLES_DIR,'/snpeff_scoring.txt')),
#    clinvar = mkScoringVector(paste0(SCORING_TABLES_DIR,'/clinvar_scoring.txt'))
# )
# #SNPEFF_SIMPLE_ANN_LOOKUP <- read.delim(paste0(SCORING_TABLES_DIR,'/snpeff_scoring.txt'))[c('ann','ann_s2')]

