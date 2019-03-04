## Parameters for detGeneStatuses.R

#--------- Paths ---------#
ROOT_DIR='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/'
if(dir.exists('/Users/lnguyen/')){
   ROOT_DIR <- paste0('/Users/lnguyen/', ROOT_DIR)
}

GENES_ENST2ENSG <- paste0(ROOT_DIR,'/data/gene_selection/human_genes_enst2ensg.txt.gz')
GENES_BED <- paste0(ROOT_DIR,'/data/gene_selection/genes.bed')

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
   full_gene_loss=30,
   loh=20,
   
   ## Additional evidence
   cn_break_in_gene_loh=0.2, ## When CN break happens with LOH, gene is more likely to be deficient
   cn_break_in_gene_germ_som=0.1,
   germ.ref_loss=0.1,
   
   ## Give bonus points when ALT AD is good
   germ.alt_exists=0.01,
   som.alt_exists=0.01
)

#--------- Cutoffs ---------#
CUTOFFS <- list(
   max.max.copy.number=0.2, ## full_gene_loss
   min.min.copy.number=5, ## amp
   max.min.minor.allele.ploidy=0.2, ## loh
   min.cn.diff.in.gene=0.9, ## cn_break_in_gene
   min.adj.tumor.ad.alt=10, ## germ.alt_exists, som.alt_exists
   min.germ.ad.diff.score=1.5, ## germ.ref_loss
   min.hit.score=6, ## hit_type
   
   min.cadd.phred=20,
   min.mcap.score=0.88,
   scap.cutoffs = (function(){
      df <- read.delim(paste0(SCORING_TABLES_DIR,'/SCAP/scap_cutoffs.txt'),stringsAsFactors=F)
      v <- df$final
      names(v) <- df$group
      return(v)
   })()
)

#--------- Options ---------#
OPTIONS <- list(
   verbose=T,
   keep.only.first.eff=T,
   overwrite.mut.profiles=T
)


