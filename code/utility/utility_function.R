# -------------------------------------------------------------------------------------------------- #
# Description: Required packages and utility function
# -------------------------------------------------------------------------------------------------- #

requiredPackages = c(
  'BiocManager',
  "msigdbr",
  'MASS',
  "readr",
  "tibble",
  "purrr",
  "data.table",
  "survival",
  "plyr",
  'dplyr'
)
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

requiredPackages_Bioc = c(
  'affy',
  'annotate',
  'AnnotationDbi',
  "clusterProfiler"
  "enrichplot",
  'genefilter', 
  'dplyr', 
  'GEOquery', 
  'hgu133a.db', 
  'hgu133plus2frmavecs', 
  'hgu133a2frmavecs', 
  'hgu133afrmavecs', 
  'hgu133plus2.db',
  'frma',
  "meta"
)

for(p in requiredPackages_Bioc){
  if(!require(p,character.only = TRUE)) BiocManager::install(p, force = T)
  library(p, character.only = TRUE)
}

#### TNBC filter function

mix.obj <- function(p, x){
  e <- p[1]*dnorm((x-p[2])/p[3])/p[3]+(1-p[1])*dnorm((x-p[4])/p[5])/p[5]
  if(any(e<=0)) Inf else -sum(log(e))
}

#### tpm convert function

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
