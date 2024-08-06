# -------------------------------------------------------------------------------------------------- #
# Description: Cox regression meta analysis
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
dir.x <- file.path(dir.results, "/GLM/")
dir.meta <- file.path(dir.results, "/meta_analysis/")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

# -------------------------------------------------------------------------------------------------- #
# All cohorts
# -------------------------------------------------------------------------------------------------- #
# Load Cox results
# -------------------------------------------------------------------------------------------------- #

cox_files <- list.files(
  path = file.path(dir.cox),
  pattern = "*.csv",
  full.names = T
)
cox_files.select <- cox_files[stringr::str_extract(cox_files, "GSE[0-9]+|VUMC") %in% datasets_selected_in_meta]

list.cox <- plyr::alply(cox_files.select, .margin = 1, function(f){
  data <- read_csv(f, col_types = readr::cols())
  dataset <- paste0(
    stringr::str_extract(basename(f),paste(datasets_selected_in_meta, collapse = "|")),
    "_",
    stringr::str_extract(toupper(basename(f)),"COX")
  )
  data <- data %>% dplyr::select(contains(c("Symbol","coef", "exp(coef)", "se(coef)","z","Pr(>|z|)")))
  data <- data %>% rename_with(
    .fn = function(x) {
      paste0(dataset,"_",x)
    },
    contains(c("coef", "exp(coef)", "se(coef)","z","Pr(>|z|)"))
  )
  data 
}
)
names(list.cox) <- stringr::str_extract(cox_files.select, "GSE[0-9]+|VUMC")

multi_cohorts <- Reduce(
  function(x,y, ...) full_join(
    x,
    y,
    ..., 
    by = "Symbol"
  ),
  list.cox
) %>% unique() 

# -------------------------------------------------------------------------------------------------- #
# meta analysis
# -------------------------------------------------------------------------------------------------- #

doParallel::registerDoParallel(6)
meta_results <- calculate_meta_analysis(multi_cohorts, sm = "HR", method.tau = "DL")
meta_results <- meta_results %>% dplyr::mutate(
  p.adj = p.adjust(pVal.random, method = "BH"),
  .after = pVal.random
) %>% dplyr::mutate(HR = exp(estimate), .after = Symbol)


write_csv(
  meta_results %>% arrange(pVal.random),
  file.path(dir.meta, "Cox_meta_analysis_all_cohorts_results.csv")
)

write_csv(
  meta_results %>% arrange(pVal.random) %>% filter(k >= 5),
  file.path(dir.meta, "Cox_meta_analysis_all_cohorts_results_filter_k_5.csv")
)

# -------------------------------------------------------------------------------------------------- #
# Exclude GSE22226
# -------------------------------------------------------------------------------------------------- #
# meta analysis
# -------------------------------------------------------------------------------------------------- #

list.cox <- list.cox[!names(list.cox) %in% "GSE22226"]
multi_cohorts <- Reduce(
  function(x,y, ...) full_join(
    x,
    y,
    ..., 
    by = "Symbol"
  ),
  list.cox
) %>% unique() 

doParallel::registerDoParallel(6)
meta_results <- calculate_meta_analysis(multi_cohorts, sm = "HR", method.tau = "DL")
meta_results <- meta_results %>% dplyr::mutate(
  p.adj = p.adjust(pVal.random, method = "BH"),
  .after = pVal.random
) %>% dplyr::mutate(HR = exp(estimate), .after = Symbol)

write_csv(
  meta_results %>% arrange(pVal.random),
  file.path(dir.meta, "Cox_meta_analysis_exclude_GSE22226_cohorts_results.csv")
)

write_csv(
  meta_results %>% arrange(pVal.random) %>% filter(k >= 5),
  file.path(dir.meta, "Cox_meta_analysis_exclude_GSE22226_cohorts_results_filter_k_5.csv")
)
