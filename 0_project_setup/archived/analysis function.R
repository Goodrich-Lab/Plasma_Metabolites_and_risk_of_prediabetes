#' @description
#' Multiple exposures and multiple outcomes analysis.
#' This script works for continuous or dichotomous outcomes.

#' @param exposures Independent variable
#' @param outcomes Dependent variable
#' @param covars Covariates adjusted for when fitting regression models
#' @param outcome_family Options: "gaussian", "binomial"; Default: "gaussian"
#' "gaussian" if outcomes are continous, "binomial" if outcomes are dichotomous
#' @param data Input data: a dataframe with all exposures, outcomes and covariates for analysis 
#' 
#' @usage 
#' data_model <- replace_with_name_of_your_data
#' ind_var <- c('pfhxs', 'pfos', 'pfda', 'pfna' ,'pfuda')
#' dep_var <- c("nafld_out_1", "nafld_out_2")
#' covar <- c("age", "sex")
#' 
#' final <- model_output(exposures = ind_var, outcomes = dep_var, covars = covar,
#' outcome_family = "binomial", data = data_model)
#' 
#' @return
#' dataframe contains analysis result ("outcome", "exposures", "term",
#' "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high",
#' "odds.ratio", "conf.low(odds.ratio)","conf.high(odds.ratio)").
#'  Notes: "odds.ratio", "conf.low (odds.ratio)", "conf.high (odds.ratio)" are not in linear regression result.

model_output <- function(exposures,
                         outcomes,
                         covars = NULL,
                         outcome_family = "gaussian",
                         data){
                             
  # check if outcomes, exposures, covars are in the dataset
  if(!all(outcomes %in% colnames(data)))
    stop("Not all outcome(s) is/are in the given data")
  if(!all(exposures %in% colnames(data)))
    stop("Not all exposure(s) is/are in the given data")
  if(!is.null(covars) & !all(covars %in% colnames(data)))
    stop("Not all covar(s) is/are in the given data")
    
  # check outcomes and exposures type
  if(length(unique(unlist(lapply(data[outcomes],class))))!=1)
    stop("Outcomes all have to be the same variable type (numeric, factor, etc)")
  if(length(unique(unlist(lapply(data[exposures],class))))!=1 & !"matrix" %in% (unique(unlist(lapply(data[exposures],class)))) )
    stop("Exposures all have to be the same variable type (numeric, factor, etc)")
    
  # Pivot longer on both outcomes and exposures
  data_l <- data %>% 
    tidylog::pivot_longer(cols = all_of(outcomes), 
                          names_to = "outcome", 
                          values_to = "outcome_value") %>% 
    tidylog::pivot_longer(cols = all_of(exposures),
                          names_to = "exposures", 
                          values_to = "exposure_concentration") %>% 
    filter(!is.na(outcome_value), 
           !is.na(exposure_concentration)) %>% 
    droplevels()
  
  
  # check number of unique values for each covar
  if(!is.null(covars)){
    
    temp <- data_l %>%
      group_by(outcome, exposures) %>%
      dplyr::summarise(across(all_of(covars),
                              ~length(unique(.)))) %>%
      ungroup()
    if(!is_empty(temp %>% dplyr::select(where(~ any(.== 1))))){
      stop("Modeling can be applied only to factors with 2 or more levels, please check covariate(s): ", 
           str_c(colnames(temp %>% dplyr::select(where(~ any(.== 1)))), collapse =", " ))
    }
    
  }
    
    # Run models with/without covars
    if(!is.null(covars)){
    # Run models with covars
    mod_output <- data_l %>% 
      group_by(outcome, exposures) %>% 
      nest() %>% 
      mutate(
        output = map(data, 
                     ~glm(as.formula(paste0("outcome_value ~ exposure_concentration + ",
                                            str_c(covars, collapse = "+"))),
                          data = .,
                          family = outcome_family) %>% 
                       tidy(.,conf.int = TRUE)))
  }else{
    # Run models without covars
    mod_output <- data_l %>% 
      group_by(outcome, exposures) %>% 
      nest() %>% 
      mutate(
        output = map(data, 
                     ~glm(as.formula(outcome_value ~ exposure_concentration),
                          data = .,
                          family = outcome_family) %>%
                       tidy(., conf.int = TRUE))) 
  }
  
  # Create final dataframe of all results
  mod_output_df <- mod_output %>% dplyr::select(-data) %>%
    unnest(output)
  
  if(outcome_family == "binomial"){
    mod_output_df <- mod_output_df %>%
      mutate(odds.ratio = exp(estimate),
             `conf.low(odds.ratio)` = exp(conf.low),
             `conf.high(odds.ratio)` = exp(conf.high))
  }
  return(mod_output_df)
}


