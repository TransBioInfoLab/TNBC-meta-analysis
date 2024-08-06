# -------------------------------------------------------------------------------------------------- #
# Description: Functions related to meta analysis
# -------------------------------------------------------------------------------------------------- #

#### Logistic regression

get_pCR_coef <- function(expr, clin, adjust = "Age"){
  
  features <- colnames(clin)
  expr_mat <- expr[rowSums(expr[,-1]) != 0, ] %>% as.data.frame()
  if(any(toupper(features) %in% toupper(adjust))){
    adj_var <- features[toupper(features) %in% toupper(adjust)]
    fn <- formula(paste("pCR ~ gene + ", paste( adj_var, collapse = "+")))
  } else {
    fn <- formula(paste("pCR ~ gene"))
  }
  print(fn)
  
  results <- plyr::adply(
    .data = expr_mat,
    .margins = 1, 
    .fun = function(e){
      Symbol <- e[1]$Symbol
      dat <- data.frame(clin, gene = t(e[-1])[,1])
      dat <- dat[!is.na(dat$pCR),]
      dat$pCR <- factor(dat$pCR)
      dat$pCR <- relevel(dat$pCR, ref = "RD")
      tryCatch({
        suppressMessages({
          g <- glm(fn, data = dat, family = 'binomial')
          if(g$converged){
            sum_results <- summary(g)
            gene_results <- t(coef(sum_results) [grep("gene",rownames(coef(sum_results)),value = TRUE),]) %>% as.data.frame()
            gene_results <- gene_results %>% dplyr::mutate(
              Symbol = Symbol,
              .before = 1
            )
          } else gene_results <- NULL
        })
        return(gene_results)
      }, error = function(e) {message(e);return(NULL)})
    }, .progress = "time",.inform = F, .parallel = T, .expand = F ,.id = "Symbol"
  )
  
  message("Finish fitting glm functions. Filter abnormal standard error..")
    
  results <- results %>% filter(
    `Std. Error` < sqrt(nrow(clin) - 1)
  ) %>% na.omit()
  
  return(results)
}

#### Cox model

get_Surv_coef <- function(expr, clin, adjust = "Age"){
  
  features <- colnames(clin)
  expr_mat <- expr[rowSums(expr[,-1]) != 0, ] %>% as.data.frame()
  if(any(toupper(features) %in% toupper(adjust))){
    adj_var <- features[toupper(features) %in% toupper(adjust)]
    adj_var_noAge <- adj_var[adj_var != "Age"]
    clin <- clin %>% mutate_at(adj_var_noAge, factor)
    fn <- formula(paste("Surv(Time, Event) ~ gene + ", paste( adj_var, collapse = "+")))
  } else {
    fn <- formula(paste("Surv(Time, Event) ~ gene"))
  }
  print(fn)
  
  results <- plyr::adply(
    .data = expr_mat,
    .margins = 1, 
    .fun = function(e){
      Symbol <- e[1]$Symbol
      dat <- data.frame(clin, gene = t(e[-1])[,1])
      tryCatch({
        suppressMessages({
          c <- coxph(fn, data = dat)
          sum_results <- summary(c)
          gene_results <- t(coef(sum_results) [grep("gene",rownames(coef(sum_results)),value = TRUE),]) %>% as.data.frame()
          gene_results <- gene_results %>% dplyr::mutate(
            Symbol = Symbol,
            .before = 1
          )
        })
        return(gene_results)
      }, error = function(e) {message(e);return(NULL)})
    }, .parallel = T, .progress = "time",.inform = F, .expand = F ,.id = "Symbol"
  )
  
  message("Finish fitting cox functions. Filter abnormal standard error..")
  
  results <- results %>% filter(
    `se(coef)` < sqrt(nrow(clin) - 1) & `Pr(>|z|)` > 1e-15
  ) %>% na.omit()
  
  return(results)
}


#### meta-analysis

calculate_meta_analysis <- function(multi_cohorts, sm = "OR", ...){
  plyr::adply(
    .data = multi_cohorts, 
    .margins = 1, 
    .fun =  function(row){
      
      est <- row[grep("Estimate|_coef",colnames(row))] %>% as.numeric
      
      direction <-  paste(
        ifelse(
          is.na(est), ".",
          ifelse(est > 0, "+", "-")
        ), collapse = "")
      
      se <- row[grep("Std. Error|se",colnames(row))] %>% as.numeric
      cohort <- gsub("_Std. Error|_se","",grep("Std. Error|se",colnames(row),value = T))
      
      df <- data.frame(
        cohort = cohort,
        est = est,
        se = se,
        stringsAsFactors = FALSE
      )
      tryCatch({
        suppressMessages({
          f <- metagen(
            TE = est,
            seTE = se,
            data = df,
            sm = sm,
            random = T,
            fixed = F,
            ...
          )
        
        est <- na.omit(est)  
        tibble::tibble(
          Symbol = row$Symbol,
          estimate = f$TE.random,
          se = f$seTE.random,
          zval.random = f$zval.random,
          pVal.random = f$pval.random,
          pValQ = f$pval.Q,
          direction = direction,
          num.coherence = ifelse(f$TE.random >= 0, length(est[est>=0]), length(est[est<0])),
          k = f$k
        ) })
      }, error = function(e) {message(e);return(NULL)})
    }  , .progress = "time",
    .parallel = T,
    .id = "Symbol",
    .expand = F
  )
}
