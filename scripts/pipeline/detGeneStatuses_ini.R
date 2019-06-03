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
   keep.only.first.eff=T,
   gene.identifier='ensembl_gene_id',
   verbose=T
)

#--------- Scoring ---------#
SCORING_TABLES_DIR <- paste0(ROOT_DIR,'/data/variant_significance/')

SCORING_CNV <- list(
   full_gene_loss=c(5,5),
   trunc=c(5,5),
   loh=c(5,0)
)

mkScoringVector <- function(scoring.table.path){
   df <- read.delim(scoring.table.path,stringsAsFactors=F)
   v <- df$score
   names(v) <- df$ann
   return(v)
}

SCORING_MUT <- list(
   ## databases
   snpeff = mkScoringVector(paste0(SCORING_TABLES_DIR,'/snpeff/snpeff_scoring.txt')),
   clinvar = mkScoringVector(paste0(SCORING_TABLES_DIR,'/clinvar/clinvar_scoring.txt')),
   enigma = mkScoringVector(paste0(SCORING_TABLES_DIR,'/enigma/enigma_scoring.txt'))
)

DIPLOTYPE_ORIGIN_RANK <- c('cnv_cnv','cnv_som','cnv_germ','germ_som')

SNPEFF_SIMPLE_ANN_LOOKUP <- read.delim(paste0(SCORING_TABLES_DIR,'/snpeff/snpeff_scoring.txt'))[c('ann','ann_s1')]

#--------- Cutoffs ---------#
CUTOFFS <- list(
   full.gene.loss.max.max.copy.number=0.3, ## full_gene_loss
   trunc.max.min.copy.number=0.3,
   loh.max.min.minor.allele.ploidy=0.2, ## loh
   min.hit.score.filter=4
)




