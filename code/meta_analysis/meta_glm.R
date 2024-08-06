# -------------------------------------------------------------------------------------------------- #
# Description: Logistic regression meta analysis
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
dir.meta <- file.path(dir.results, "/meta_analysis/")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

# -------------------------------------------------------------------------------------------------- #
# All cohorts
# -------------------------------------------------------------------------------------------------- #
# Load GLM results
# -------------------------------------------------------------------------------------------------- #

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

glm_files <- list.files(
  path = file.path(dir.glm),
  pattern = "*.csv",
  full.names = T
)
glm_files.select <- glm_files[stringr::str_extract(glm_files, "GSE[0-9]+|VUMC") %in% datasets_selected_in_meta]

list.glm <- plyr::alply(glm_files.select, .margin = 1, function(f){
    data <- read_csv(f, col_types = readr::cols())
    dataset <- paste0(
      stringr::str_extract(basename(f),paste(datasets_selected_in_meta, collapse = "|")),
      "_",
      stringr::str_extract(toupper(basename(f)),"GLM")
    )
    data <- data %>% dplyr::select(contains(c("Symbol","Std. Error","z value","Pr(>|z|)","Estimate")))
    data <- data %>% rename_with(
      .fn = function(x) {
        paste0(dataset,"_",x)
      },
      contains(c("Std. Error","z value","Pr(>|z|)","Estimate"))
    )
    data 
  }
)
names(list.glm) <- stringr::str_extract(glm_files.select, "GSE[0-9]+|VUMC")

multi_cohorts <- Reduce(
  function(x,y, ...) full_join(
    x,
    y,
    ..., 
    by = "Symbol"
  ),
  list.glm
) %>% unique() 

#### Select genes that contain more than 4 cohorts
#multi_cohorts_k <- aaply(
#  multi_cohorts, .margin = 1, function(f){
#    length(na.omit(f[,grep("Estimate",colnames(f),value = T)]  %>% as.numeric()))
#  }, .expand = F
#) %>% as.numeric()
#multi_cohorts_k_5 <- multi_cohorts[multi_cohorts_k >= 5,]
# -------------------------------------------------------------------------------------------------- #
# meta analysis
# -------------------------------------------------------------------------------------------------- #

doParallel::registerDoParallel(6)
meta_results <- calculate_meta_analysis(multi_cohorts, method.tau = "DL")
meta_results <- meta_results %>% dplyr::mutate(
  p.adj = p.adjust(pVal.random, method = "BH"),
  .after = pVal.random
) %>% dplyr::mutate(odds.ratio = exp(estimate), .after = Symbol)


write_csv(
  meta_results %>% arrange(pVal.random),
  file.path(dir.meta, "GLM_meta_analysis_all_cohorts_results.csv")
)

write_csv(
  meta_results %>% arrange(pVal.random) %>% filter(k >= 5),
  file.path(dir.meta, "GLM_meta_analysis_all_cohorts_results_filter_k_5.csv")
)

# -------------------------------------------------------------------------------------------------- #
# Exclude GSE22226
# -------------------------------------------------------------------------------------------------- #
# Load GLM results
# -------------------------------------------------------------------------------------------------- #

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
  "GSE192341",
  "GSE163882",
  "GSE154524"
)

glm_files <- list.files(
  path = file.path(dir.glm),
  pattern = "*.csv",
  full.names = T
)
glm_files.select <- glm_files[stringr::str_extract(glm_files, "GSE[0-9]+|VUMC") %in% datasets_selected_in_meta]

list.glm <- plyr::alply(glm_files.select, .margin = 1, function(f){
  data <- read_csv(f, col_types = readr::cols())
  dataset <- paste0(
    stringr::str_extract(basename(f),paste(datasets_selected_in_meta, collapse = "|")),
    "_",
    stringr::str_extract(toupper(basename(f)),"GLM")
  )
  data <- data %>% dplyr::select(contains(c("Symbol","Std. Error","z value","Pr(>|z|)","Estimate")))
  data <- data %>% rename_with(
    .fn = function(x) {
      paste0(dataset,"_",x)
    },
    contains(c("Std. Error","z value","Pr(>|z|)","Estimate"))
  )
  data 
}
)
names(list.glm) <- stringr::str_extract(glm_files.select, "GSE[0-9]+|VUMC")

multi_cohorts <- Reduce(
  function(x,y, ...) full_join(
    x,
    y,
    ..., 
    by = "Symbol"
  ),
  list.glm
) %>% unique() 

# -------------------------------------------------------------------------------------------------- #
# meta analysis
# -------------------------------------------------------------------------------------------------- #

doParallel::registerDoParallel(6)
meta_results <- calculate_meta_analysis(multi_cohorts, method.tau = "DL")
meta_results <- meta_results %>% dplyr::mutate(
  p.adj = p.adjust(pVal.random, method = "BH"),
  .after = pVal.random
) %>% dplyr::mutate(odds.ratio = exp(estimate), .after = Symbol)


write_csv(
  meta_results %>% arrange(pVal.random),
  file.path(dir.meta, "GLM_meta_analysis_exclude_GSE22226_cohorts_results.csv")
)

write_csv(
  meta_results %>% arrange(pVal.random) %>% filter(k >= 5),
  file.path(dir.meta, "GLM_meta_analysis_exclude_GSE22226_cohorts_results_filter_k_5.csv")
)


