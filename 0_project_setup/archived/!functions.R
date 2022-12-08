# Functions


# Read Data from HPC -----------------------------------------------------
read_data_from_hpc <- function(file_path, n_col_in_df){
  
  # Read in data, without headers.
  df_for_colnames <- read.table(file_path,
                                sep = ",", 
                                na.strings = "",
                                as.is = TRUE, 
                                fill = TRUE,
                                header = FALSE)
  
  # Get temp col names. Col names are not indexed correctly though- they 
  # Need to be shifted to the right by 1, and we need to 
  row1 <- df_for_colnames[1,] %>% as.character(.)
  # Create column names for readining in new data
  col_names <- c("model_term",
                 row1[-length(row1)], 
                 paste0(row1[1:n_col_in_df], ".999999"))
  
  # Read in final data with headers
  suppressWarnings(
    df_original <- read.table(file_path,
                              sep = ",",col.names = col_names,
                              na.strings = "",
                              as.is = TRUE, fill = TRUE, header = TRUE) %>% 
      dplyr::select(-model_term))
  
  # Clean Col names
  df_clean_names <- df_original %>% 
    janitor::clean_names() %>% 
    rename_all(~str_replace(., "x2_5", "lcl") %>% 
                 str_replace(., "x97_5", "ucl") %>% 
                 str_replace(., "p_val", "p") %>% 
                 str_replace(., "var_names", "var"))
  
  # Get new column names in a dataframe
  new_colnames <- tibble(cnms = colnames(df_clean_names)) %>% 
    mutate(variable = str_split_fixed(cnms, "_", n = 2)[,1], 
           group = str_split_fixed(cnms, "_", n = 2)[,2] %>% 
             if_else(. == "", "0", .) %>% 
             as.numeric)
  
  # create a list of colnames by column group
  cnms_bygroup <- split(new_colnames, new_colnames$group)
  
  # create a list of sol result by column groups
  results_list <- cnms_bygroup %>% 
    modify(~dplyr::select(df_clean_names, all_of(.$cnms)))
  
  # rename all cols to match names, then cbind
  results_df <- results_list %>% 
    modify(~rename_all(., ~str_split_fixed(., "_", 2)[,1])) %>% 
    bind_rows() %>% 
    filter(!is.na(metabolite))
  
  # Separate beta and pips
  results_final <- results_df %>% 
    mutate(term = case_when(str_detect(var, "beta") ~ "beta", 
                            str_detect(var, "gamma") ~ "pip",
                            str_detect(var, "psi") ~ "beta", 
                            str_detect(var, "eta") ~ "eta", 
                            TRUE ~ var), 
           var = str_remove(var, ".beta") %>% 
             str_remove(".gamma") %>% 
             str_remove("eta.") %>%
             str_replace("psi", "mixture")) %>% 
    rename(exposure = var, 
           estimate = mean, 
           p_value = p) %>% 
    dplyr::select(exposure, term, everything())
  
  
  return(results_final)
}

# Create Functions for summarizing Data ---------------
transpose_ft <- function(ft) {
  dataout <- ft %>%
    mutate(name = str_c(mz, time, sep = "_")) %>%
    select(name, everything(), -mz, -time) %>%
    gather(file_name, val, 2:ncol(.)) %>%
    spread(name, val) %>% 
    mutate(file_name = str_remove(file_name, "_mz_xml"))
  
  return(dataout)
}

# # Rename Compounds -----------------------
rename_pfas <- function(pfas_names, include_asterisk = FALSE, 
                        arrange_by_class = FALSE){
  x <- tibble(pfas = pfas_names)
  
  suppressWarnings(
    pfas2 <-  x %>%
      mutate(pfas2 = case_when(
        pfas == "pfhxs" ~ "PFHxS",
        pfas == "pfhps" ~ "PFHpS",
        pfas == "pfpes" ~ "PFPeS",
        pfas == "pfhpa" ~ "PFHpA",
        pfas == "nmefosaab" ~ "N-MeFOSAA-b†", 
        pfas == "pfuda" ~ "PFUnDA†",
        pfas == "pfds" ~ "PFDS†",
        pfas == "netfosaa" ~ "N-EtFOSAA†",
        pfas == "pfns" ~ "PFNS†",
        pfas == "pfbs" ~ "PFBS†",
        pfas == "x82fts" ~ "8:2 FTS†", 
        pfas == "pfhxa" ~ "PFHxA†", 
        pfas == "pfdoa" ~ "PFDoDA†",
        pfas == "Mixture effect" ~ "Mixture effect",
        TRUE ~ toupper(pfas)) %>% 
          as.factor() %>% 
          fct_relevel(., 
                      "PFOS", "PFOA", "PFHxS", "PFNA", "PFHpS","PFDA", "PFPeS", 
                      "PFHpA","N-MeFOSAA-b†","N-EtFOSAA†","PFDS†","PFBS†", 
                      "8:2 FTS†", "PFDoDA†", "PFUnDA†","PFNS†","PFHxA†",
                      "Mixture effect")) 
    )
    
    if(include_asterisk == TRUE){ 
      suppressWarnings(
        
       pfas2 <-  x %>%
        mutate(pfas2 = case_when(
          pfas == "pfhxs" ~ "PFHxS",
          pfas == "pfhps" ~ "PFHpS",
          pfas == "pfpes" ~ "PFPeS",
          pfas == "pfhpa" ~ "PFHpA",
          pfas == "nmefosaab" ~ "N-MeFOSAA-b*", 
          pfas == "pfuda" ~ "PFUnDA*",
          pfas == "pfds" ~ "PFDS*",
          pfas == "netfosaa" ~ "N-EtFOSAA*",
          pfas == "pfns" ~ "PFNS*",
          pfas == "pfbs" ~ "PFBS*",
          pfas == "x82fts" ~ "8:2 FTS*", 
          pfas == "pfhxa" ~ "PFHxA*", 
          pfas == "pfdoa" ~ "PFDoDA*", 
          pfas == "Mixture effect" ~ "Mixture effect",
          TRUE ~ toupper(pfas)) %>% 
            as.factor() %>% 
            fct_relevel(., 
                        "PFOS", "PFOA", "PFHxS", "PFNA", "PFHpS","PFDA", "PFPeS",
                        "PFHpA","N-MeFOSAA-b*","N-EtFOSAA*","PFDS*","PFBS*", 
                        "8:2 FTS*", "PFDoDA*", "PFUnDA*","PFNS*","PFHxA*", 
                        "Mixture effect")) )
    }
    
    if(arrange_by_class == TRUE){ 
      suppressWarnings(
        
      pfas2 <-  pfas2 %>% 
        # left_join(lod, by = "pfas") %>% 
        mutate(pfas2 = fct_relevel(pfas2, 
                                   "PFOS", 
                                   "PFHxS",
                                   "PFHpS",
                                   "PFPeS",
                                   "PFOA", 
                                   "PFNA", 
                                   "PFDA", 
                                   "PFHpA")))
    }
    
    return(pfas2$pfas2)
}


# mxtune_v11 - force the covariates be in the model


mxtune<-function(X,Y,U = NULL, Z = as.matrix(rep(1, ncol(X))),
                 c = 1/2,
                 alpha.est.init = rep(0,ncol(Z)+1),
                 epsilon = 6,
                 max_s = 20,
                 margin_s = 0.00001,
                 maxstep = 100,
                 margin = 0.001,
                 maxstep_inner = 100,
                 tol.inner = 0.001,
                 compute.likelihood = FALSE,
                 verbosity = 1,
                 standardize=T,intercept = T){
  
  
  ## Initialize alpha
  n = nrow(X);p=ncol(X);k = length(unique(Y)); q = ncol(Z)-1
  alpha.est.old = alpha.est.init
  s = 1
  likelihood.score = c()
  #beta.star = rep(0,p)
  #scale.func <- function(x) sqrt(sum((x-mean(x))^2)/length(x))
  #sx <- scale(X,scale=apply(X, 2, scale.func))
  #sy <- as.vector(scale(Y, scale=scale.func(Y)))
  
  ## calculate alpha.max
  a.max = NULL
  for (i in 1:k) {
    y.temp = ifelse(Y==unique(Y)[i],1,0)
    a.max = cbind(a.max, max( abs(t(y.temp - mean(y.temp)*(1-mean(y.temp))) %*% X ) )/ (c * n))
  } 
  
  alpha.max = min(max(abs(a.max))/abs(ifelse(Z==0,0.01,Z)))
  alpha.min = alpha.max * ifelse(n<p, 0.01, 0.0001)
  
  ## reparameterize Z
  Z.original = Z
  Z = sweep(Z,2,colMeans(Z))
  Z = cbind(rep(1,p), Z)
  
  
  while (s < max_s){ ## iteratively compute alpha and beta mode
    ## Given alpha, update beta mode, Y hat etc
    lambda = exp(Z%*%alpha.est.old)
    gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c)) 
    A = (lambda^2*c^2 + 2*lambda*(1-c))/2 # V
    
    ## use glmnet to update beta.mode
    mode.est=coef(glmnet(X,Y,alpha=0, family = 'multinomial',penalty.factor = A,lambda = sum(A)/p/n, standardize = F, intercept = FALSE))
    
    # cat.prob = NULL
    # for (i in 1:k) {
    #   cat.prob = cbind(cat.prob,exp(X%*%mode.est[[i]][-1]))
    # }
    cat.prob = sapply(mode.est, function(g) exp(X%*%g[-1]))
    t.mode = cat.prob/rowSums(cat.prob) # softmax function
    
    B_inv = apply(t.mode, 2, function(g) as.vector(1/(g*(1-g))))
    
    Y_HAT = NULL
    for (i in 1:k) {
      # B_inv = cBind(B_inv, as.vector(1/(t.mode[,i]*(1-t.mode[,i]))))
      y = ifelse(Y==names(mode.est)[i],1,0)
      cat.y = X%*%mode.est[[i]][-1] + B_inv[,i]*(y - t.mode[,i])
      Y_HAT = cbind(Y_HAT,cat.y)
    } 
    
    
    if(compute.likelihood == TRUE){
      V = diag(as.vector(gamma_inv))
      likelihood.score = c(Likelihood.multinomial(mode.est,X,Y,V),likelihood.score)
    }
    
    ## Given beta mode, Y hat etc, update alpha
    cat(yellow$italic$bold("Start estimating alpha:\n"))
    alpha.est.new <- estimate.alpha(X,Y_HAT,Z,c = c,B_inv,alpha.init = alpha.est.old, alpha.max, alpha.min,
                                    epsilon,
                                    maxstep = maxstep, 
                                    margin = margin,
                                    maxstep_inner = maxstep_inner,
                                    tol.inner = tol.inner,
                                    compute.likelihood = F,verbosity = 1)
    print(alpha.est.new)
    
    # Check convergence 
    if(sum(abs(alpha.est.new - alpha.est.old)) < margin_s ){
      cat(red$bold("Done!\n"))
      break
    }
    cat("Difference between alpha_old and alpha_new:",sum(abs(alpha.est.new - alpha.est.old)),"\n")
    
    alpha.est.old <- alpha.est.new
    cat(green$italic("#---"),green$bold("Outer loop Iteration",s,"Done"),green$italic("---#\n"),sep = "")
    s <- s+1
  }
  
  ## restore alpha and Z
  #alpha.est.old = alpha.est.old[-1]
  #Z = Z[,-1]
  
  tauEst = exp(Z%*%alpha.est.old)
  pen_vec= tauEst/n
  pen_vec_cov = c(pen_vec, rep(0,ncol(U)))
  
  pen_vec_cov[pen_vec_cov > 1e3] <- 1e3
  C = sum(pen_vec_cov)/p
  
  obj = glmnet(cbind(X, U),Y,family = "multinomial",alpha = c, lambda = C, penalty.factor = pen_vec_cov, type.measure="mae")
  cus.coef = coef(obj,x=cbind(X, U),y=Y,family = "multinomial",alpha = c, exact=TRUE,s= C, penalty.factor = pen_vec_cov,standardize=standardize,intercept = intercept)
  
  return(list(model = obj, C = C,pen_vec = pen_vec_cov, cus.coef = cus.coef,alpha.hat = alpha.est.old,n_iter = s-1,likelihood.score = likelihood.score))
  
}

estimate.alpha <- function(X,Y,Z,c,B_inv,
                           alpha.init, alpha.max, alpha.min,
                           epsilon,
                           maxstep,
                           margin,
                           maxstep_inner,
                           tol.inner,
                           compute.likelihood,
                           verbosity){
  n = nrow(X);p=ncol(X);q = ncol(Z);k = ncol(B_inv)
  ## Initialize
  alpha.old = alpha.init
  itr = 1
  
  while(itr < maxstep){
    # Given alpha, update theta
    lambda = exp(Z%*%alpha.old)
    gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
    
    # Sigma_y = NULL; theta = NULL
    # for (i in 1:k) {
    #   Sigma_y[[i]] = diag(B_inv[,i]) + (t(t(X)*c(gamma)))%*%t(X)
    #   #print(Sigma_y[[i]])
    #   #print(solve(Sigma_y[[i]],X))
    #   theta = rbind(theta, colSums(X*solve(Sigma_y[[i]],X)))  
    # }
    
    Sigma_y = apply(B_inv, 2 , function(g) list(diag(g) + (t(t(X)*c(gamma)))%*%t(X)))
    theta = t(sapply(Sigma_y, function(g) colSums(X*solve(g[[1]],X))))
    
    
    # Given theta, update alpha
    update.result <-update_alpha(X,Y,Z,c = c,B_inv = B_inv, alpha.old = alpha.old,alpha.max, alpha.min, theta = theta, epsilon, maxstep_inner = maxstep_inner,tol.inner = tol.inner)
    alpha.new <- update.result$alpha.est
    
    # Check convergence 
    if(sum(abs(alpha.new - alpha.old)) < margin ){
      #cat(green$bold("Innerloop Done!\n"))
      break
    }
    alpha.old <- alpha.new
    cat(blue$italic("#-----------------"),magenta("Inner loop Iteration",itr,"Done"),blue$italic("-----------------#\n"),sep = "")
    itr <- itr+1
    #print(k)
    
  }
  return(alpha.new)
}

update_alpha<-function(X,Y,Z,c,B_inv,alpha.old,alpha.max, alpha.min, theta, epsilon, maxstep_inner,tol.inner){
  ## initial
  alpha.iner.old = alpha.old
  i_inner = 1
  n=nrow(X)
  p=ncol(X)
  k = ncol(B_inv)
  B = 1/B_inv
  while (i_inner < maxstep_inner){
    # given alpha update delta
    lambda = exp(Z%*%alpha.iner.old)
    gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
    
    X_B_sqrt = NULL; Y_B_sqrt = NULL; sd_y = NULL; C = NULL;delta.est = NULL
    for (i in 1:k) {
      X_B_sqrt[[i]] = X*sqrt(B[,i])
      Y_B_sqrt[[i]] = Y[,i]*sqrt(B[,i])
      sd_y[[i]] <- sqrt(var(Y_B_sqrt[[i]])*(n-1)/n)
      C[[i]]=sum(1/gamma)/p* sd_y[[i]]*1/n
      
      delta.est[[i]]=coef(glmnet(X_B_sqrt[[i]],Y_B_sqrt[[i]],alpha=0, penalty.factor = 1/gamma,lambda = C[[i]], standardize = F, intercept = FALSE))[-1]
    }
    
    ## given delta update alpha with constrain
    #print(alpha.old[,i],likelihood.ungrpalpha.theta,likelihood.ungrpalpha.theta.gradient,theta = theta[i,],delta=delta.est[[i]])
    alpha.iner.new <- optim(alpha.old,likelihood.alpha.theta,likelihood.alpha.theta.gradient,c =c,Z=Z,theta = theta,delta=delta.est,
                            k = k,method = "L-BFGS-B", upper = c(alpha.max*epsilon, rep(Inf, length(alpha.old)-1)))$par
    if (sum(abs(alpha.iner.new - alpha.iner.old)) < tol.inner){
      break
    }
    i_inner = i_inner + 1
    alpha.iner.old <- alpha.iner.new
  }
  return(list(alpha.est=alpha.iner.old,inner_iter = i_inner))
}

likelihood.alpha.theta<-function(Z,c,alpha,theta,delta,k){
  lambda = exp(Z %*% alpha) 
  gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
  
  L.theta.alpha = NULL
  for (i in 1:k) {
    L.theta.alpha[i] = t(theta[i,])%*%gamma + delta[[i]]^2%*%(1/gamma)
  }
  
  return(as.numeric(sum(L.theta.alpha)))
  print(as.numeric(L.theta.alpha))
}

likelihood.alpha.theta.gradient<-function(Z,c,alpha,theta,delta,k){
  
  lambda = exp(Z %*% alpha) 
  gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
  
  dev_gamma_alpha = as.vector((-2*(2*(1-c) + 2*c^2*lambda))/(2*lambda*(1-c) + (c*lambda)^2)^2 *lambda) *Z
  dev_gamma = NULL;L.prime = NULL
  for (i in 1:k) {
    dev_gamma[[i]] = theta[i,] - delta[[i]]^2/(gamma^2)
    # L.prime = rbind(L.prime, crossprod(dev_gamma[[i]], dev_gamma_alpha))
  }
  L.prime = t(sapply(dev_gamma, function(g) crossprod(g, dev_gamma_alpha)))
  
  return(colSums(L.prime))
  print(L.prime)
}
