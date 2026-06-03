select_KL_edpm_vb <- function(data,
                              num_vars,
                              cat_vars,
                              outcome,
                              eventvar,
                              event_only,
                              K_values,
                              L_values,
                              max_iter,
                              tol,
                              n_starts,
                              seed,
                              eps,
                              verbose) {
  
  
  if (missing(data)) stop("Input 'data' must be supplied.")
  if (missing(num_vars)) stop("Input 'num_vars' must be supplied.")
  if (missing(cat_vars)) stop("Input 'cat_vars' must be supplied.")
  if (missing(outcome)) stop("Input 'outcome' must be supplied.")
  if (missing(eventvar)) stop("Input 'eventvar' must be supplied.")
  if (missing(event_only)) stop("Input 'event_only' must be supplied.")
  if (missing(K_values)) stop("Input 'K_values' must be supplied.")
  if (missing(L_values)) stop("Input 'L_values' must be supplied.")
  if (missing(max_iter)) stop("Input 'max_iter' must be supplied.")
  if (missing(tol)) stop("Input 'tol' must be supplied.")
  if (missing(n_starts)) stop("Input 'n_starts' must be supplied.")
  if (missing(seed)) stop("Input 'seed' must be supplied.")
  if (missing(eps)) stop("Input 'eps' must be supplied.")
  if (missing(verbose)) stop("Input 'verbose' must be supplied.")
  
  required_vars <- c(num_vars, cat_vars, outcome, eventvar)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(
      "The following variables are missing from the data: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  set.seed(seed)
  
  # -----------------------------
  # 1) Build analysis data
  # -----------------------------
  df_model <- data[, required_vars, drop = FALSE]
  df_model <- stats::na.omit(df_model)
  
  for (v in cat_vars) {
    df_model[[v]] <- as.factor(df_model[[v]])
  }
  
  df_model[[outcome]] <- as.numeric(df_model[[outcome]])
  df_model[[eventvar]] <- as.integer(df_model[[eventvar]])
  
  df_model <- df_model[df_model[[outcome]] > 0, , drop = FALSE]
  
  if (event_only) {
    df_model <- df_model[df_model[[eventvar]] == 1, , drop = FALSE]
  }
  
  N <- nrow(df_model)
  if (N == 0) stop("No observations remain after preprocessing.")
  
  if (verbose) {
    cat("Number of observations used:", N, "\n")
  }
  
  # -----------------------------
  # 2) Standardize numeric predictors
  # -----------------------------
  X_num_raw <- as.matrix(df_model[, num_vars, drop = FALSE])
  
  scale_center <- colMeans(X_num_raw)
  scale_scale <- apply(X_num_raw, 2, stats::sd)
  
  if (any(scale_scale == 0)) {
    stop(
      "The following numeric predictors have zero standard deviation: ",
      paste(names(scale_scale)[scale_scale == 0], collapse = ", ")
    )
  }
  
  X_num <- sweep(X_num_raw, 2, scale_center, "-")
  X_num <- sweep(X_num, 2, scale_scale, "/")
  
  # -----------------------------
  # 3) Dummy-code categorical predictors for regression
  # -----------------------------
  cat_formula <- stats::as.formula(
    paste("~", paste(cat_vars, collapse = " + "))
  )
  
  dummies <- stats::model.matrix(cat_formula, data = df_model)[, -1, drop = FALSE]
  
  X_reg <- as.matrix(cbind(X_num, dummies))
  
  # Add intercept column
  X_design <- cbind(Intercept = 1, X_reg)
  P <- ncol(X_design)
  
  # -----------------------------
  # 4) Categorical variables for enrichment part
  # -----------------------------
  X_cat_df <- df_model[, cat_vars, drop = FALSE]
  cat_levels <- lapply(X_cat_df, levels)
  cat_dims <- sapply(cat_levels, length)
  
  for (j in seq_along(cat_vars)) {
    X_cat_df[[j]] <- as.integer(X_cat_df[[j]])
  }
  
  X_cat <- as.matrix(X_cat_df)
  P_cat <- ncol(X_cat)
  max_cat <- max(cat_dims)
  
  # -----------------------------
  # 5) Log-response
  # -----------------------------
  y_log <- log(df_model[[outcome]] + eps)
  
  # -----------------------------
  # Helper: log-sum-exp
  # -----------------------------
  log_sum_exp <- function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  }
  
  # -----------------------------
  # Helper: Dirichlet ELBO contribution
  # -----------------------------
  dirichlet_elbo_term <- function(alpha_prior, alpha_post) {
    
    E_log_theta <- digamma(alpha_post) - digamma(sum(alpha_post))
    
    log_p_theta <- lgamma(sum(alpha_prior)) -
      sum(lgamma(alpha_prior)) +
      sum((alpha_prior - 1) * E_log_theta)
    
    log_q_theta <- lgamma(sum(alpha_post)) -
      sum(lgamma(alpha_post)) +
      sum((alpha_post - 1) * E_log_theta)
    
    log_p_theta - log_q_theta
  }
  
  # -----------------------------
  # Helper: fit one K,L pair
  # -----------------------------
  fit_one_KL <- function(K, L, start_id) {
    
    set.seed(seed + 1000 * K + 100 * L + start_id)
    
    # q(z_i = k, s_i = l)
    r <- array(stats::runif(N * K * L), dim = c(N, K, L))
    
    for (i in seq_len(N)) {
      r[i, , ] <- r[i, , ] / sum(r[i, , ])
    }
    
    beta <- matrix(0, nrow = K, ncol = P)
    sigma2 <- rep(stats::var(y_log), K)
    
    # Priors from EDPM code:
    # intercept[k] ~ N(0, 10^2)
    # beta[k,p]   ~ N(0, 5^2)
    # log_sigma_y[k] ~ N(0, 1.5^2)
    beta_prior_precision <- diag(c(1 / 10^2, rep(1 / 5^2, P - 1)))
    
    alpha_pi_prior <- rep(1, K)
    alpha_omega_prior <- rep(1, L)
    
    alpha_psi_prior <- array(0, dim = c(P_cat, max_cat))
    
    for (j in seq_len(P_cat)) {
      alpha_psi_prior[j, 1:cat_dims[j]] <- 1
    }
    
    objective_trace <- numeric(max_iter)
    
    for (iter in seq_len(max_iter)) {
      
      # --------------------------------------------------------
      # Expected counts
      # --------------------------------------------------------
      N_kl <- matrix(0, nrow = K, ncol = L)
      
      for (k in seq_len(K)) {
        for (l in seq_len(L)) {
          N_kl[k, l] <- sum(r[, k, l])
        }
      }
      
      N_k <- rowSums(N_kl)
      
      # --------------------------------------------------------
      # q(pi)
      # --------------------------------------------------------
      alpha_pi <- alpha_pi_prior + N_k
      E_log_pi <- digamma(alpha_pi) - digamma(sum(alpha_pi))
      
      # --------------------------------------------------------
      # q(omega)
      # --------------------------------------------------------
      alpha_omega <- matrix(0, nrow = K, ncol = L)
      E_log_omega <- matrix(0, nrow = K, ncol = L)
      
      for (k in seq_len(K)) {
        alpha_omega[k, ] <- alpha_omega_prior + N_kl[k, ]
        E_log_omega[k, ] <- digamma(alpha_omega[k, ]) -
          digamma(sum(alpha_omega[k, ]))
      }
      
      # --------------------------------------------------------
      # q(psi)
      # --------------------------------------------------------
      alpha_psi <- array(0, dim = c(K, L, P_cat, max_cat))
      E_log_psi <- array(0, dim = c(K, L, P_cat, max_cat))
      
      for (k in seq_len(K)) {
        for (l in seq_len(L)) {
          for (j in seq_len(P_cat)) {
            
            counts_j <- rep(0, cat_dims[j])
            
            for (c in seq_len(cat_dims[j])) {
              counts_j[c] <- sum(r[X_cat[, j] == c, k, l])
            }
            
            alpha_psi[k, l, j, 1:cat_dims[j]] <-
              alpha_psi_prior[j, 1:cat_dims[j]] + counts_j
            
            E_log_psi[k, l, j, 1:cat_dims[j]] <-
              digamma(alpha_psi[k, l, j, 1:cat_dims[j]]) -
              digamma(sum(alpha_psi[k, l, j, 1:cat_dims[j]]))
          }
        }
      }
      
      # --------------------------------------------------------
      # Update regression parameters for each Y-cluster
      # This uses the actual Normal priors from the EDPM model.
      # No ridge penalty is added.
      # --------------------------------------------------------
      for (k in seq_len(K)) {
        
        w_k <- apply(r[, k, , drop = FALSE], 1, sum)
        W_sum <- sum(w_k)
        
        if (W_sum < eps) {
          
          beta[k, ] <- rep(0, P)
          sigma2[k] <- stats::var(y_log)
          
        } else {
          
          # Weighted Gaussian regression update with Normal prior
          XtW <- t(X_design) * (w_k / sigma2[k])
          
          A <- XtW %*% X_design + beta_prior_precision
          b <- XtW %*% y_log
          
          beta[k, ] <- as.numeric(solve(A, b))
          
          resid_k <- as.numeric(y_log - X_design %*% beta[k, ])
          sigma2[k] <- sum(w_k * resid_k^2) / W_sum
          sigma2[k] <- max(sigma2[k], eps)
        }
      }
      
      # --------------------------------------------------------
      # Update q(z_i, s_i)
      # --------------------------------------------------------
      for (i in seq_len(N)) {
        
        log_r_i <- matrix(0, nrow = K, ncol = L)
        
        for (k in seq_len(K)) {
          
          mu_ik <- sum(X_design[i, ] * beta[k, ])
          
          log_y <- -0.5 * log(2 * pi * sigma2[k]) -
            0.5 * ((y_log[i] - mu_ik)^2 / sigma2[k])
          
          for (l in seq_len(L)) {
            
            log_x <- 0
            
            for (j in seq_len(P_cat)) {
              xij <- X_cat[i, j]
              log_x <- log_x + E_log_psi[k, l, j, xij]
            }
            
            log_r_i[k, l] <- E_log_pi[k] +
              E_log_omega[k, l] +
              log_y +
              log_x
          }
        }
        
        normalizer <- log_sum_exp(as.vector(log_r_i))
        r[i, , ] <- exp(log_r_i - normalizer)
      }
      
      # --------------------------------------------------------
      # Compute VB objective
      # --------------------------------------------------------
      objective <- 0
      
      for (i in seq_len(N)) {
        for (k in seq_len(K)) {
          
          mu_ik <- sum(X_design[i, ] * beta[k, ])
          
          log_y <- -0.5 * log(2 * pi * sigma2[k]) -
            0.5 * ((y_log[i] - mu_ik)^2 / sigma2[k])
          
          for (l in seq_len(L)) {
            
            if (r[i, k, l] > 0) {
              
              log_x <- 0
              
              for (j in seq_len(P_cat)) {
                xij <- X_cat[i, j]
                log_x <- log_x + E_log_psi[k, l, j, xij]
              }
              
              objective <- objective +
                r[i, k, l] * (
                  E_log_pi[k] +
                    E_log_omega[k, l] +
                    log_y +
                    log_x -
                    log(r[i, k, l])
                )
            }
          }
        }
      }
      
      # Dirichlet prior/posterior terms
      objective <- objective + dirichlet_elbo_term(alpha_pi_prior, alpha_pi)
      
      for (k in seq_len(K)) {
        objective <- objective +
          dirichlet_elbo_term(alpha_omega_prior, alpha_omega[k, ])
      }
      
      for (k in seq_len(K)) {
        for (l in seq_len(L)) {
          for (j in seq_len(P_cat)) {
            objective <- objective +
              dirichlet_elbo_term(
                alpha_psi_prior[j, 1:cat_dims[j]],
                alpha_psi[k, l, j, 1:cat_dims[j]]
              )
          }
        }
      }
      
      # Actual Normal prior contribution from EDPM model
      for (k in seq_len(K)) {
        
        objective <- objective +
          stats::dnorm(beta[k, 1], mean = 0, sd = 10, log = TRUE)
        
        if (P > 1) {
          objective <- objective +
            sum(stats::dnorm(beta[k, -1], mean = 0, sd = 5, log = TRUE))
        }
        
        log_sigma_k <- log(sqrt(sigma2[k]))
        
        objective <- objective +
          stats::dnorm(log_sigma_k, mean = 0, sd = 1.5, log = TRUE)
      }
      
      objective_trace[iter] <- objective
      
      if (iter > 1) {
        diff_obj <- abs(objective_trace[iter] - objective_trace[iter - 1])
        
        if (diff_obj < tol) {
          objective_trace <- objective_trace[seq_len(iter)]
          break
        }
      }
    }
    
    # Posterior modal cluster assignments
    z_prob <- matrix(0, nrow = N, ncol = K)
    
    for (k in seq_len(K)) {
      z_prob[, k] <- apply(r[, k, , drop = FALSE], 1, sum)
    }
    
    z_hat <- apply(z_prob, 1, which.max)
    
    s_hat <- rep(NA_integer_, N)
    
    for (i in seq_len(N)) {
      k_i <- z_hat[i]
      s_hat[i] <- which.max(r[i, k_i, ])
    }
    
    list(
      K = K,
      L = L,
      start_id = start_id,
      final_objective = tail(objective_trace, 1),
      objective_trace = objective_trace,
      beta = beta,
      sigma2 = sigma2,
      r = r,
      z_hat = z_hat,
      s_hat = s_hat,
      alpha_pi = alpha_pi,
      alpha_omega = alpha_omega,
      alpha_psi = alpha_psi
    )
  }
  
  # -----------------------------
  # Grid search over K and L
  # -----------------------------
  all_fits <- list()
  result_rows <- list()
  counter <- 1
  
  for (K in K_values) {
    for (L in L_values) {
      
      if (verbose) {
        cat("\nFitting EDPM VB for K =", K, ", L =", L, "\n")
      }
      
      best_fit_this_KL <- NULL
      best_objective_this_KL <- -Inf
      
      for (start_id in seq_len(n_starts)) {
        
        if (verbose) {
          cat("  Start", start_id, "of", n_starts, "\n")
        }
        
        fit <- fit_one_KL(K, L, start_id)
        
        if (fit$final_objective > best_objective_this_KL) {
          best_objective_this_KL <- fit$final_objective
          best_fit_this_KL <- fit
        }
      }
      
      all_fits[[counter]] <- best_fit_this_KL
      
      result_rows[[counter]] <- data.frame(
        K = K,
        L = L,
        VB_objective = best_objective_this_KL,
        n_iter = length(best_fit_this_KL$objective_trace)
      )
      
      counter <- counter + 1
    }
  }
  
  results_table <- do.call(rbind, result_rows)
  results_table <- results_table[order(-results_table$VB_objective), ]
  rownames(results_table) <- NULL
  
  best_K <- results_table$K[1]
  best_L <- results_table$L[1]
  
  best_index <- which(
    sapply(all_fits, function(fit) {
      fit$K == best_K &&
        fit$L == best_L &&
        fit$final_objective == results_table$VB_objective[1]
    })
  )[1]
  
  best_fit <- all_fits[[best_index]]
  
  if (verbose) {
    cat("\nSelected K =", best_K, "\n")
    cat("Selected L =", best_L, "\n")
    cat("Best VB objective =", results_table$VB_objective[1], "\n")
  }
  
  return(
    list(
      best_K = best_K,
      best_L = best_L,
      results_table = results_table,
      best_fit = best_fit,
      processed_data = df_model,
      X_design = X_design,
      X_reg = X_reg,
      X_cat = X_cat,
      y_log = y_log,
      scale_center = scale_center,
      scale_scale = scale_scale,
      cat_levels = cat_levels
    )
  )
}