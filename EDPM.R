edpm_predict_log <- function(data,
                             num_vars,
                             cat_vars,
                             outcome,
                             eventvar,
                             event_only,
                             K,
                             L,
                             niter,
                             nburnin,
                             thin,
                             nchains,
                             seed,
                             eps,
                             return_details,
                             return_samples,
                             verbose) {
  
  # -----------------------------
  # 0) Package and input checks
  # -----------------------------
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("Package 'nimble' is required.")
  }
  
  if (missing(data)) stop("Input 'data' must be supplied.")
  if (missing(num_vars)) stop("Input 'num_vars' must be supplied.")
  if (missing(cat_vars)) stop("Input 'cat_vars' must be supplied.")
  if (missing(outcome)) stop("Input 'outcome' must be supplied.")
  if (missing(eventvar)) stop("Input 'eventvar' must be supplied.")
  if (missing(event_only)) stop("Input 'event_only' must be supplied.")
  if (missing(K)) stop("Input 'K' must be supplied.")
  if (missing(L)) stop("Input 'L' must be supplied.")
  if (missing(niter)) stop("Input 'niter' must be supplied.")
  if (missing(nburnin)) stop("Input 'nburnin' must be supplied.")
  if (missing(thin)) stop("Input 'thin' must be supplied.")
  if (missing(nchains)) stop("Input 'nchains' must be supplied.")
  if (missing(seed)) stop("Input 'seed' must be supplied.")
  if (missing(eps)) stop("Input 'eps' must be supplied.")
  if (missing(return_details)) stop("Input 'return_details' must be supplied.")
  if (missing(return_samples)) stop("Input 'return_samples' must be supplied.")
  if (missing(verbose)) stop("Input 'verbose' must be supplied.")
  
  required_vars <- c(num_vars, cat_vars, outcome, eventvar)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(
      "The following variables are missing from the data: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  if (K < 1) stop("'K' must be at least 1.")
  if (L < 1) stop("'L' must be at least 1.")
  if (niter <= nburnin) stop("'niter' must be larger than 'nburnin'.")
  if (thin < 1) stop("'thin' must be at least 1.")
  if (nchains < 1) stop("'nchains' must be at least 1.")
  
  set.seed(seed)
  
  # -----------------------------
  # 1) Build model data
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
  P_reg <- ncol(X_reg)
  
  # -----------------------------
  # 4) Categorical variables as integers for enrichment part
  # -----------------------------
  X_cat_df <- df_model[, cat_vars, drop = FALSE]
  cat_levels <- lapply(X_cat_df, levels)
  cat_dims <- sapply(cat_levels, length)
  
  for (j in seq_along(cat_vars)) {
    X_cat_df[[j]] <- as.integer(X_cat_df[[j]])
  }
  
  X_cat <- as.matrix(X_cat_df)
  P_cat <- ncol(X_cat)
  
  # -----------------------------
  # 5) Response on log scale
  # -----------------------------
  y_log <- log(df_model[[outcome]] + eps)
  
  # -----------------------------
  # 6) Dirichlet hyperparameters
  # -----------------------------
  max_cat <- max(cat_dims)
  
  alpha_phi <- array(0, dim = c(P_cat, max_cat))
  for (j in seq_len(P_cat)) {
    alpha_phi[j, 1:cat_dims[j]] <- 1.0
  }
  
  data_list <- list(
    y_log = as.numeric(y_log),
    X_reg = X_reg,
    X_cat = X_cat,
    alpha = rep(1, K),
    alpha_x = rep(1, L),
    alpha_phi = alpha_phi
  )
  
  constants <- list(
    N = N,
    P_reg = P_reg,
    P_cat = P_cat,
    K = K,
    L = L,
    cat_dims = as.integer(cat_dims)
  )
  
  # -----------------------------
  # 7) EDPM nimble model
  # -----------------------------
  edpm_code <- nimble::nimbleCode({
    for (i in 1:N) {
      
      z[i] ~ dcat(pi[1:K])
      s[i] ~ dcat(omega[z[i], 1:L])
      
      mu_y[i] <- intercept[z[i]] +
        inprod(beta[z[i], 1:P_reg], X_reg[i, 1:P_reg])
      
      y_log[i] ~ dnorm(mu_y[i], sd = sigma_y[z[i]])
      
      for (j in 1:P_cat) {
        X_cat[i, j] ~ dcat(psi[z[i], s[i], j, 1:cat_dims[j]])
      }
    }
    
    pi[1:K] ~ ddirch(alpha[1:K])
    
    for (k in 1:K) {
      
      omega[k, 1:L] ~ ddirch(alpha_x[1:L])
      
      intercept[k] ~ dnorm(0, sd = 10)
      
      for (p in 1:P_reg) {
        beta[k, p] ~ dnorm(0, sd = 5)
      }
      
      log_sigma_y[k] ~ dnorm(0, sd = 1.5)
      sigma_y[k] <- exp(log_sigma_y[k])
      
      for (l in 1:L) {
        for (j in 1:P_cat) {
          psi[k, l, j, 1:cat_dims[j]] ~
            ddirch(alpha_phi[j, 1:cat_dims[j]])
        }
      }
    }
  })
  
  # -----------------------------
  # 8) Initial values
  # -----------------------------
  psi_init <- array(0, dim = c(K, L, P_cat, max_cat))
  
  for (k in 1:K) {
    for (l in 1:L) {
      for (j in 1:P_cat) {
        psi_init[k, l, j, 1:cat_dims[j]] <- rep(1 / cat_dims[j], cat_dims[j])
      }
    }
  }
  
  inits <- list(
    z = sample(1:K, N, replace = TRUE),
    s = sample(1:L, N, replace = TRUE),
    pi = rep(1 / K, K),
    omega = matrix(1 / L, nrow = K, ncol = L),
    intercept = stats::rnorm(K, 0, 0.5),
    beta = matrix(stats::rnorm(K * P_reg, 0, 0.3),
                  nrow = K,
                  ncol = P_reg),
    log_sigma_y = stats::rnorm(K, 0, 0.2),
    psi = psi_init
  )
  
  # -----------------------------
  # 9) Build and run MCMC
  # -----------------------------
  if (verbose) {
    cat("Building NIMBLE model...\n")
  }
  
  model <- nimble::nimbleModel(
    code = edpm_code,
    data = data_list,
    constants = constants,
    inits = inits
  )
  
  cmodel <- nimble::compileNimble(model)
  
  monitors <- c("pi", "omega", "intercept", "beta", "log_sigma_y", "z", "s")
  
  conf <- nimble::configureMCMC(model, monitors = monitors)
  mcmc <- nimble::buildMCMC(conf)
  cmcmc <- nimble::compileNimble(mcmc, project = model)
  
  if (verbose) {
    cat("Running MCMC...\n")
  }
  
  samples <- nimble::runMCMC(
    cmcmc,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = nchains,
    samplesAsCodaMCMC = TRUE
  )
  
  # -----------------------------
  # 10) Combine chains into one matrix
  # -----------------------------
  if (inherits(samples, "mcmc.list")) {
    sample_mat <- do.call(rbind, lapply(samples, as.matrix))
  } else {
    sample_mat <- as.matrix(samples)
  }
  
  # -----------------------------
  # 11) Helper functions for posterior summaries
  # -----------------------------
  order_by_first_index <- function(cols) {
    ind <- sapply(cols, function(x) {
      as.integer(regmatches(x, gregexpr("[0-9]+", x))[[1]][1])
    })
    cols[order(ind)]
  }
  
  mode_assign <- function(mat) {
    apply(mat, 2, function(v) {
      tab <- table(v)
      as.integer(names(tab)[which.max(tab)])
    })
  }
  
  extract_vector_mean <- function(sample_mat, varname, len) {
    cols <- grep(paste0("^", varname, "\\["), colnames(sample_mat), value = TRUE)
    
    out <- numeric(len)
    
    for (cc in cols) {
      ind <- as.integer(regmatches(cc, gregexpr("[0-9]+", cc))[[1]][1])
      out[ind] <- mean(sample_mat[, cc])
    }
    
    out
  }
  
  extract_matrix_mean <- function(sample_mat, varname, nrow_out, ncol_out) {
    cols <- grep(paste0("^", varname, "\\["), colnames(sample_mat), value = TRUE)
    
    out <- matrix(NA_real_, nrow = nrow_out, ncol = ncol_out)
    
    for (cc in cols) {
      inds <- as.integer(regmatches(cc, gregexpr("[0-9]+", cc))[[1]])
      out[inds[1], inds[2]] <- mean(sample_mat[, cc])
    }
    
    out
  }
  
  # -----------------------------
  # 12) Modal cluster assignments
  # -----------------------------
  z_cols <- grep("^z\\[", colnames(sample_mat), value = TRUE)
  s_cols <- grep("^s\\[", colnames(sample_mat), value = TRUE)
  
  z_cols <- order_by_first_index(z_cols)
  s_cols <- order_by_first_index(s_cols)
  
  z_all <- sample_mat[, z_cols, drop = FALSE]
  s_all <- sample_mat[, s_cols, drop = FALSE]
  
  z_mode <- mode_assign(z_all)
  s_mode <- mode_assign(s_all)
  
  # -----------------------------
  # 13) Posterior mean regression parameters
  # -----------------------------
  intercept_mean <- extract_vector_mean(sample_mat, "intercept", K)
  beta_mean <- extract_matrix_mean(sample_mat, "beta", K, P_reg)
  
  # -----------------------------
  # 14) EDPM prediction on log scale
  # -----------------------------
  log_prediction <- numeric(N)
  
  for (i in 1:N) {
    k <- z_mode[i]
    log_prediction[i] <- intercept_mean[k] + sum(beta_mean[k, ] * X_reg[i, ])
  }
  
  names(log_prediction) <- rownames(df_model)
  
  if (verbose) {
    cat("Finished. Returning EDPM predictions on log scale.\n")
  }
  
  # -----------------------------
  # 15) Output
  # -----------------------------
  if (!return_details) {
    return(log_prediction)
  }
  
  output <- list(
    log_prediction = log_prediction,
    processed_data = df_model,
    X_reg = X_reg,
    X_cat = X_cat,
    z_mode = z_mode,
    s_mode = s_mode,
    posterior_mean_intercept = intercept_mean,
    posterior_mean_beta = beta_mean,
    scale_center = scale_center,
    scale_scale = scale_scale,
    cat_levels = cat_levels,
    K = K,
    L = L,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = nchains
  )
  
  if (return_samples) {
    output$samples <- samples
  }
  
  return(output)
}