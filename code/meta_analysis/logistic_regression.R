# -------------------------------------------------------------------------------------------------- #
# Description: Logistic regression 
# -------------------------------------------------------------------------------------------------- #
# Load packages and create directory
# -------------------------------------------------------------------------------------------------- #

setwd(".")
dir.base <- getwd()
source(file.path(dir.base, "code/utility/utility_function.R"))
source(file.path(dir.base, "code/utility/meta_function.R"))
dir.exp <- file.path(dir.base, "/Data/Expr/")
dir.clinical <- file.path(dir.base, "/Data/clinical/")
dir.results <- file.path(dir.base, "Results/")
dir.glm <- file.path(dir.results, "/GLM/")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

# -------------------------------------------------------------------------------------------------- #
# Load expression data and cliincal data
# -------------------------------------------------------------------------------------------------- #

#### Datasets selected in meta analysis

datasets_selected_in_meta <- c(
  "GSE25055", 
  "GSE25065",
  "GSE41998",
  "GSE16446",
  "GSE18864",
  "GSE18728",
  "GSE20194",
  "GSE20271",
  "GSE22093",
  "GSE32646",
  "VUMC",
  "GSE22226",
  "GSE192341",
  "GSE163882",
  "GSE164458"
)

#### Load expression datasets

exp_files <- list.files(
  path = file.path(dir.exp),
  pattern = "*.csv",
  full.names = T
)
exp_files.select <- exp_files[stringr::str_extract(exp_files, "GSE[0-9]+|VUMC") %in% datasets_selected_in_meta]

list.exp <- plyr::alply(exp_files.select, .margin = 1, function(f){
    exp <- read_csv(f)
  }
)
names(list.exp) <- stringr::str_extract(exp_files.select, "GSE[0-9]+|VUMC")

#### Load clinical datasets

clinical_files <- list.files(
  path = file.path(dir.clinical),
  pattern = "*.csv",
  full.names = T
)
clinical_files.select <- clinical_files[stringr::str_extract(clinical_files, "GSE[0-9]+|VUMC") %in% datasets_selected_in_meta]

list.clinical <- plyr::alply(clinical_files.select, .margin = 1, function(f){
    clinical <- read_csv(f)
  }
)

names(list.clinical) <- stringr::str_extract(clinical_files.select, "GSE[0-9]+|VUMC")
list.clinical <- list.clinical[names(list.exp)]
list.clinical$GSE154524$Age <- as.factor(list.clinical$GSE154524$Age)

# -------------------------------------------------------------------------------------------------- #
# Logistic regression
# -------------------------------------------------------------------------------------------------- #

doParallel::registerDoParallel(cores = 6)

for(i in c(1:length(list.exp))){
  message("Begin fitting glm functions for ", names(list.exp)[i])
  
  results <- get_pCR_coef(
    expr = list.exp[[i]],
    clin = list.clinical[[i]], 
    adjust = "Age"
  )
  
  write_csv(
    results, 
    file.path(dir.glm, paste0("GLM_", names(list.exp)[i], ".csv"))
  )
}


