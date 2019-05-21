## Parameters for detGeneStatuses.R

options(stringsAsFactors = F)

#--------- Paths ---------#
ROOT_DIR='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'
if(dir.exists('/Users/lnguyen/')){
   ROOT_DIR <- paste0('/Users/lnguyen/', ROOT_DIR)
}

GENES_HGNC <- read.delim(paste0(ROOT_DIR,'/data/gene_selection/hgnc_gene_names_20190311.txt'))

#GENES_ENST2ENSG <- paste0(ROOT_DIR,'/data/gene_selection/human_genes_enst2ensg.txt.gz')

#--------- Options ---------#
OPTIONS <- list(
   overwrite.mut.profile=F,
   overwrite.biall.mut.profile=T,
   keep.only.first.eff=T,
   gene.identifier='ensembl_gene_id',
   ignore.additional.evidence=T,
   simplify.snpeff.eff=F,
   verbose=T
)

#--------- Scoring ---------#
SCORING_TABLES_DIR <- paste0(ROOT_DIR,'/data/variant_significance/')

mkScoringVector <- function(scoring.table.path){
   df <- read.delim(scoring.table.path,stringsAsFactors=F)
   v <- df$score
   names(v) <- df$ann
   return(v)
}

SCORING <- list(
   ## databases
   snpeff = mkScoringVector(paste0(SCORING_TABLES_DIR,'/snpeff/snpeff_scoring.txt')),
   clinvar = mkScoringVector(paste0(SCORING_TABLES_DIR,'/clinvar/clinvar_scoring.txt')),
   enigma = mkScoringVector(paste0(SCORING_TABLES_DIR,'/enigma/enigma_scoring.txt')),

   ## Use integer values for main evidance
   full_gene_loss=10,
   loh=5,
   
   ## Only used by getGeneMaxEff() for prioritizing as shown below
   score_boost.full_gene_loss=300,
   score_boost.loh_som=200, ## LOH + som VUS likely has more impact than LOH + germ VUS
   score_boost.loh_germ=100,

   ## Additional evidence
   germ.is_hotspot_mut=0.1,
   som.is_hotspot_mut=0.1,
   
   cn_break_in_gene=0.01,
   
   germ.ref_loss=0.01, ## Check if adjusted REF AD is high --> no ref loss
   som.ref_loss=0.01,
   
   germ.alt_exists=0.001, ## Give bonus points when ALT AD is good
   som.alt_exists=0.001
)

SNPEFF_SIMPLE_ANN_LOOKUP <- read.delim(paste0(SCORING_TABLES_DIR,'/snpeff/snpeff_scoring.txt'))[c('ann','ann_s2')]

IS_DEF_MIN_HIT_SCORE <- c(
   'full_gene_loss' = SCORING$full_gene_loss,

   'loh+som' = SCORING$loh + 5 + 
      SCORING$som.is_hotspot_mut +
      SCORING$som.alt_exists +
      min(SCORING$som.ref_loss, SCORING$cn_break_in_gene),

   'loh+germ' = SCORING$loh + 5 +
      SCORING$germ.is_hotspot_mut +
      SCORING$germ.alt_exists + 
      min(SCORING$germ.ref_loss, SCORING$cn_break_in_gene),

   'germ+som' = 5 + 5 +
      SCORING$germ.is_hotspot_mut + SCORING$germ.alt_exists + min(SCORING$germ.ref_loss, SCORING$cn_break_in_gene) +
      SCORING$som.is_hotspot_mut + SCORING$som.alt_exists + min(SCORING$som.ref_loss, SCORING$cn_break_in_gene)
)

#--------- Cutoffs ---------#
CUTOFFS <- list(
   max.max.copy.number=0.2, ## full_gene_loss
   min.min.copy.number=5, ## amp
   max.min.minor.allele.ploidy=0.2, ## loh
   min.cn.diff.in.gene=0.9, ## cn_break_in_gene
   min.adj.tumor.ad.alt=10, ## germ.alt_exists, som.alt_exists
   min.ad.diff.score=1.5, ## germ.ref_loss, som.ref_loss
   min.hit.score.filter=4

   # min.cadd.phred=20,
   # min.mcap.score=0.88,
   # scap.cutoffs = (function(){
   #    df <- read.delim(paste0(SCORING_TABLES_DIR,'/SCAP/scap_cutoffs.txt'),stringsAsFactors=F)
   #    v <- df$final
   #    names(v) <- df$group
   #    return(v)
   # })()
)




