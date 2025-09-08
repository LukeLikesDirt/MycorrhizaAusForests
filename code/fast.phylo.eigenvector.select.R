library(foreach)
library(doParallel)
library(Matrix)

# Fixed statistics with without tree

# Fast Pagel's λ test via eigen-decomposition of the VCV #######################
eigenPagelPrep <- function(vcv) {
  n <- nrow(vcv)
  eig <- eigen(vcv, symmetric = TRUE)
  h <- vcv[1,1]  # Assuming all diagonals are equal
  ones <- rep(1, n)
  z <- crossprod(eig$vectors, ones)
  list(U = eig$vectors, d = eig$values, z = z, h = h)
}

lambdaTest_fast <- function(x, vcv, prep = NULL) {
  # Pre-compute eigen decomposition once
  if (is.null(prep)) {
    prep <- eigenPagelPrep(vcv)
  }
  
  U <- prep$U; d <- prep$d; z <- prep$z; h <- prep$h
  y <- crossprod(U, x)        # y = U^T %*% x
  n <- length(x)
  lambda.max <- max(vcv)/max(vcv[lower.tri(vcv)])
  
  # Profile log-likelihood as a function of λ
  logLik_lambda <- function(lambda) {
    v <- lambda * d + (1 - lambda) * h
    if (any(v <= 0)) return(-Inf)
    inv_v <- 1 / v
    sum_y2_v <- sum(y^2 * inv_v)
    sum_z2_v <- sum(z^2 * inv_v)
    sum_zy_v <- sum(y * z * inv_v)
    Q <- sum_y2_v - (sum_zy_v)^2 / sum_z2_v
    if (Q <= 0) return(-Inf)
    sigma2 <- Q / n
    logdetV <- sum(log(v))
    -0.5 * (n * log(sigma2) + logdetV)
  }
  
  # Optimize over [0, lambda.max]
  opt <- suppressWarnings( optimize(logLik_lambda, interval = c(0, lambda.max), maximum = TRUE) )
  Lambda <- opt$maximum
  L1 <- opt$objective
  L0 <- logLik_lambda(0)
  pvalue <- pchisq(2 * (L1 - L0), df = 1, lower.tail = FALSE)
  list(Lambda = Lambda, pvalue = pvalue, prep = prep)
}

# Function for building oriAbouheif PEMs #######################################

# # Example:
# # Generate eigenvector matrices using the oriAbouheif proximity matrix:
# pem <- build_oriAbouheif_PEMs(phylo_tree, prefix = "PEM")
# 
# # Verify the attributes are set correctly 
# print(head(attr(pem, "values")))
# print(paste("Number of positive PEMs:", ncol(pem)))

build_oriAbouheif_PEMs <- function(tree, threshold = 0, prefix = "PEM") {
  # Compute oriAbouheif proximities
  prox <- proxTips(tree, method = "oriAbouheif", normalize = "none")
  
  # Apply threshold if specified (like MPSEM's 'a' parameter)
  if (threshold > 0) {
    prox[prox < threshold] <- 0
  }
  
  # Ensure symmetry
  prox <- (prox + t(prox)) / 2
  
  # Eigendecomposition
  eig_result <- eigen(prox, symmetric = TRUE)
  
  # Keep positive eigenvalues (standard for PEMs)
  positive_idx <- which(eig_result$values > .Machine$double.eps^0.5)
  
  if (length(positive_idx) == 0) {
    warning("No positive eigenvalues found for tree")
    return(NULL)
  }
  
  eigenvectors <- eig_result$vectors[, positive_idx, drop = FALSE]
  eigenvalues <- eig_result$values[positive_idx]
  
  # Convert to data.frame format (matching your original approach)
  pem_df <- as.data.frame(eigenvectors)
  
  # Set proper column names (matching your naming convention)
  colnames(pem_df) <- paste0(prefix, seq_len(ncol(pem_df)))
  
  # Set row names from tree tips
  rownames(pem_df) <- rownames(prox)
  
  # Set attributes (matching your original structure)
  attr(pem_df, "values") <- eigenvalues
  names(attr(pem_df, "values")) <- colnames(pem_df)
  attr(pem_df, "method") <- "oriAbouheif"
  attr(pem_df, "proximity_matrix") <- prox
  
  return(pem_df)
}
# Method-specific test statistic computations ##################################
compute_test_statistic <- function(x, prox_sparse, method = "Abouheif") {
  n <- length(x)
  x_centered <- x - mean(x)
  
  # Method-specific computations
  if (method == "Abouheif") {
    # Standard Abouheif test with frequency matrix transformation
    total_sum <- sum(prox_sparse)
    if (total_sum == 0) return(0)
    
    A <- prox_sparse / total_sum
    d_i <- Matrix::rowSums(A)
    
    if (any(d_i == 0)) {
      warning("Some rows have zero sum in proximity matrix")
      d_i[d_i == 0] <- 1e-10
    }
    
    D_inv_sqrt <- 1 / sqrt(d_i)
    x_d_centered <- x_centered * D_inv_sqrt
    
    numerator <- as.numeric(crossprod(x_d_centered, A %*% x_d_centered))
    denominator <- sum(x_d_centered^2)
    
  } else if (method == "oriAbouheif") {
    # Original Abouheif test statistic (no frequency transformation)
    # Direct application of proximity weights
    numerator <- as.numeric(crossprod(x_centered, prox_sparse %*% x_centered))
    denominator <- sum(x_centered^2)
    
  } else if (method %in% c("patristic", "nNodes", "sumDD")) {
    # For distance-based methods, use Moran's I formulation
    # Convert distances to weights (inverse relationship)
    W <- prox_sparse
    if (method %in% c("patristic", "nNodes", "sumDD")) {
      # For distance matrices, convert to similarity weights
      max_dist <- max(W)
      if (max_dist > 0) {
        W <- max_dist - W
        W[W < 0] <- 0  # Ensure non-negative weights
      }
    }
    
    # Row-standardize the weight matrix
    row_sums <- Matrix::rowSums(W)
    if (any(row_sums == 0)) {
      warning("Some rows have zero sum in weight matrix")
      row_sums[row_sums == 0] <- 1e-10
    }
    
    # Create row-standardized matrix
    W_std <- W / row_sums
    
    # Moran's I computation
    numerator <- as.numeric(crossprod(x_centered, W_std %*% x_centered))
    denominator <- sum(x_centered^2)
    
  } else {
    stop("Unknown method: ", method)
  }
  
  if (denominator == 0) return(0)
  return(numerator / denominator)
}

# Optimised phylo.moran.test with method-specific computations ################
phylo.moran.test_fast <- function(x, prox, method = "Abouheif", nperm = 999, one.tailed = TRUE) {
  if (!is.numeric(x)) x <- as.numeric(x)
  
  # Validate method
  valid_methods <- c("oriAbouheif", "patristic", "nNodes", "Abouheif", "sumDD")
  if (!method %in% valid_methods) {
    stop("Method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Pre-compute sparse matrix if not already sparse
  if (!inherits(prox, "sparseMatrix")) {
    prox <- Matrix::Matrix(prox, sparse = TRUE)
  }
  
  # Observed statistic using method-specific computation
  obs_stat <- compute_test_statistic(x, prox, method = method)
  
  # Pre-allocate result vector
  sim_stats <- numeric(nperm)
  
  # Vectorized permutation test with optimized sampling
  n <- length(x)
  
  # Use faster sampling and vectorized operations
  for (i in seq_len(nperm)) {
    # Efficient random permutation
    perm_idx <- sample.int(n, n, replace = FALSE)
    x_perm <- x[perm_idx]
    sim_stats[i] <- compute_test_statistic(x_perm, prox, method = method)
  }
  
  # Vectorized p-value calculation
  p_value <- if (one.tailed) {
    (sum(sim_stats >= obs_stat) + 1) / (nperm + 1)
  } else {
    (sum(abs(sim_stats) >= abs(obs_stat)) + 1) / (nperm + 1)
  }
  
  list(obs = obs_stat, sim = sim_stats, pvalue = p_value)
}

# Even faster version using parallel permutations ##############################
phylo.moran.test_parallel <- function(x, prox, method = "Abouheif", nperm = 999, one.tailed = TRUE, ncores = NULL) {
  if (!is.numeric(x)) x <- as.numeric(x)
  
  # Validate method
  valid_methods <- c("oriAbouheif", "patristic", "nNodes", "Abouheif", "sumDD")
  if (!method %in% valid_methods) {
    stop("Method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Pre-compute sparse matrix if not already sparse
  if (!inherits(prox, "sparseMatrix")) {
    prox <- Matrix::Matrix(prox, sparse = TRUE)
  }
  
  # Observed statistic using method-specific computation
  obs_stat <- compute_test_statistic(x, prox, method = method)
  
  # Setup parallel backend
  if (is.null(ncores)) ncores <- min(parallel::detectCores() - 1, 4)
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl))
  
  # Export necessary objects
  parallel::clusterExport(cl, c("compute_test_statistic", "prox", "x", "method"), envir = environment())
  parallel::clusterEvalQ(cl, library(Matrix))
  
  # Parallel permutation test
  sim_stats <- parallel::parSapply(cl, seq_len(nperm), function(i) {
    x_perm <- sample(x)
    compute_test_statistic(x_perm, prox, method = method)
  })
  
  # Calculate p-value
  p_value <- if (one.tailed) {
    (sum(sim_stats >= obs_stat) + 1) / (nperm + 1)
  } else {
    (sum(abs(sim_stats) >= abs(obs_stat)) + 1) / (nperm + 1)
  }
  
  list(obs = obs_stat, sim = sim_stats, pvalue = p_value)
}

# Abouheif preprocessing with method-specific proximity computation ############
precompute_abouheif_matrices <- function(phylo_matrix, tree = NULL, method = "Abouheif") {
  # Try to extract or compute proximity matrix
  if (method == "Abouheif" && !is.null(attr(phylo_matrix, "proximity_matrix"))) {
    # Use pre-computed proximity matrix from phylo_matrix attributes
    prox <- attr(phylo_matrix, "proximity_matrix")
  } else if (method %in% c("oriAbouheif", "nNodes", "patristic", "sumDD")) {
    # Check if phylo_matrix has the required proximity matrix stored
    if (!is.null(attr(phylo_matrix, "proximity_matrix")) && 
        !is.null(attr(phylo_matrix, "method")) && 
        attr(phylo_matrix, "method") == method) {
      # Use pre-computed proximity matrix with matching method
      prox <- attr(phylo_matrix, "proximity_matrix")
    } else if (!is.null(tree)) {
      # Fall back to computing from tree if available
      if (!requireNamespace("adephylo", quietly = TRUE)) {
        stop("Package 'adephylo' is required for method = '", method, "' when proximity matrix is not pre-computed")
      }
      
      # Map method names to proxTips methods
      prox_method <- switch(method,
                            "oriAbouheif" = "oriAbouheif",
                            "nNodes" = "nNodes", 
                            "patristic" = "patristic",
                            "sumDD" = "sumDD"
      )
      
      # Compute proximity using adephylo::proxTips
      prox <- adephylo::proxTips(tree, method = prox_method, normalize = "none")
    } else {
      # Try to infer proximity matrix from phylo_matrix based on method
      prox <- infer_proximity_from_matrix(phylo_matrix, method)
    }
  } else {
    # Use provided matrix directly (assume it's already a proximity matrix)
    prox <- phylo_matrix
  }
  
  # Ensure symmetry and convert to sparse
  prox <- (prox + t(prox)) / 2
  prox_sparse <- Matrix::Matrix(prox, sparse = TRUE)
  
  # Pre-compute any invariant quantities
  list(
    prox_sparse = prox_sparse,
    n = nrow(prox_sparse),
    method = method
  )
}

# Helper function to infer proximity matrix from phylogenetic matrix ##########
infer_proximity_from_matrix <- function(phylo_matrix, method) {
  # This function attempts to derive proximity matrices from the phylogenetic matrix
  # when tree is not available
  
  if (method == "oriAbouheif") {
    # For oriAbouheif, if phylo_matrix is a VCV matrix, we can derive proximities
    # The oriAbouheif proximity is related to the phylogenetic covariance structure
    if (is.matrix(phylo_matrix) && isSymmetric(phylo_matrix)) {
      # Assume phylo_matrix is VCV matrix and derive oriAbouheif proximities
      # oriAbouheif proximities are proportional to shared evolutionary history
      prox <- phylo_matrix
      # Normalize to [0,1] range typically expected for oriAbouheif
      max_val <- max(prox)
      if (max_val > 0) {
        prox <- prox / max_val
      }
    } else {
      stop("Cannot infer oriAbouheif proximities from provided matrix. Please provide tree object or pre-computed proximity matrix.")
    }
  } else if (method == "patristic") {
    # For patristic distances, convert VCV to distances
    if (is.matrix(phylo_matrix) && isSymmetric(phylo_matrix)) {
      # Patristic distance = diag(i) + diag(j) - 2*VCV(i,j)
      n <- nrow(phylo_matrix)
      diag_vals <- diag(phylo_matrix)
      prox <- matrix(0, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          prox[i,j] <- diag_vals[i] + diag_vals[j] - 2 * phylo_matrix[i,j]
        }
      }
      rownames(prox) <- rownames(phylo_matrix)
      colnames(prox) <- colnames(phylo_matrix)
    } else {
      stop("Cannot infer patristic distances from provided matrix. Please provide tree object or pre-computed proximity matrix.")
    }
  } else if (method == "nNodes") {
    # nNodes requires tree structure - cannot be inferred from VCV alone
    stop("Method 'nNodes' requires tree object or pre-computed proximity matrix with node counts.")
  } else if (method == "sumDD") {
    # sumDD requires tree structure - cannot be inferred from VCV alone  
    stop("Method 'sumDD' requires tree object or pre-computed proximity matrix.")
  } else {
    # Default: use phylo_matrix as-is
    prox <- phylo_matrix
  }
  
  return(prox)
}

# Main function ################################################################
fast.phylo.eigenvector.select <- function(
    x, phylo_eigenvectors, phylo_matrix,
    method = c("abouheif", "lambda"),
    abouheif_method = c("Abouheif", "oriAbouheif", "patristic", "nNodes", "sumDD"),
    tree = NULL,
    PE.all = FALSE, nperm = 999, nperm.global = 999,
    alpha = 0.05, verbose = FALSE,
    max_selection_pe = NULL, max_output_pe = NULL,
    use_parallel_moran = TRUE
) {
  method <- match.arg(method)
  abouheif_method <- match.arg(abouheif_method)
  
  x <- x[rownames(phylo_eigenvectors)]
  moran_values <- attr(phylo_eigenvectors, "values")
  positive_idx <- which(moran_values > 0)
  
  # Input validation
  if (any(is.na(x))) stop("NA entries in x")
  if (length(x) != nrow(phylo_eigenvectors)) stop("Length of x must match rows in phylo_eigenvectors")
  if (is.null(names(x))) stop("x must be a named vector")
  if (!all(names(x) %in% rownames(phylo_eigenvectors))) stop("Names in x must match rownames in phylo_eigenvectors")
  if (is.null(max_selection_pe)) max_selection_pe <- length(positive_idx)
  if (is.null(max_output_pe)) max_output_pe <- length(positive_idx)
  if (max_selection_pe < max_output_pe) stop("max_selection_pe must be >= max_output_pe")
  
  # Check if tree or pre-computed proximities are available for certain methods
  if (method == "abouheif" && abouheif_method %in% c("nNodes", "sumDD") && 
      is.null(tree) && 
      (is.null(attr(phylo_matrix, "proximity_matrix")) || 
       is.null(attr(phylo_matrix, "method")) || 
       attr(phylo_matrix, "method") != abouheif_method)) {
    stop("Tree object or pre-computed proximity matrix is required for abouheif_method = '", abouheif_method, "'")
  }
  
  if (length(positive_idx) == 0) {
    if (verbose) message("No eigenvectors with positive phylogenetic autocorrelation found")
    result <- list(global.test = list(obs = NA, pvalue = 1))
    if (PE.all) result$PE.all <- phylo_eigenvectors
    return(result)
  }
  
  # Pre-filter PEMs by R-squared (vectorized)
  PEM_positive_temp <- phylo_eigenvectors[, positive_idx, drop = FALSE]
  r_squared <- apply(PEM_positive_temp, 2, function(me) {
    # Faster R-squared calculation
    fit <- .lm.fit(cbind(1, me), x)
    1 - sum(fit$residuals^2) / sum((x - mean(x))^2)
  })
  
  top_pe_idx <- order(r_squared, decreasing = TRUE)[1:min(max_selection_pe, length(r_squared))]
  positive_idx <- positive_idx[top_pe_idx]
  PEM_positive <- phylo_eigenvectors[, positive_idx, drop = FALSE]
  positive_moran <- moran_values[positive_idx]
  
  # Prepare proximity matrices or VCV for different methods
  if (method == "abouheif") {
    abouheif_prep <- precompute_abouheif_matrices(phylo_matrix, tree, abouheif_method)
    prox_sparse <- abouheif_prep$prox_sparse
  } else {
    vcv <- phylo_matrix
    prep <- eigenPagelPrep(vcv)
  }
  
  # Global test with optimized functions
  if (method == "abouheif") {
    if (use_parallel_moran && nperm.global > 500) {
      global_test <- phylo.moran.test_parallel(x, prox_sparse, method = abouheif_method, nperm = nperm.global)
    } else {
      global_test <- phylo.moran.test_fast(x, prox_sparse, method = abouheif_method, nperm = nperm.global)
    }
  } else {
    tmp <- lambdaTest_fast(x, vcv, prep)
    global_test <- list(obs = tmp$Lambda, pvalue = tmp$pvalue)
  }
  
  res <- if (PE.all) list(global.test = global_test, PE.all = phylo_eigenvectors) else list(global.test = global_test)
  if (global_test$pvalue >= alpha) {
    if (verbose) message("No significant phylogenetic signal")
    return(res)
  }
  
  if (verbose) message(sprintf(
    "Global test significant (p = %.4f for %s%s). Starting PEM selection with %d variables.",
    global_test$pvalue, method, 
    if(method == "abouheif") paste0(" (", abouheif_method, ")") else "",
    ncol(PEM_positive)
  ))
  
  p <- global_test$pvalue
  PE.sel <- min.stat <- p.vector <- integer(0)
  current_x <- x
  
  # Setup parallel and ensure cleanup
  cores <- parallel::detectCores() - 1
  doParallel::registerDoParallel(cores)
  on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  
  # Pre-compute design matrices for faster fitting
  design_matrices <- lapply(seq_len(ncol(PEM_positive)), function(i) {
    cbind(1, PEM_positive[, i])
  })
  
  # Selection loop with optimizations
  while (p < alpha && length(PE.sel) < max_output_pe) {
    stat.vector <- foreach::foreach(
      i = seq_len(ncol(PEM_positive)),
      .combine = c,
      .packages = c("Matrix"),
      .export = c("lambdaTest_fast", "compute_test_statistic", "phylo.moran.test_fast", "abouheif_method")
    ) %dopar% {
      if (i %in% PE.sel) return(NA)
      
      X <- design_matrices[[i]]
      temp_fit <- .lm.fit(X, current_x)  # Faster than lm.fit
      temp_resid <- current_x - X %*% temp_fit$coefficients
      names(temp_resid) <- names(current_x)
      
      if (method == "abouheif") {
        # Use method-specific computation (no permutation needed here)
        compute_test_statistic(temp_resid, prox_sparse, method = abouheif_method)
      } else {
        lambdaTest_fast(temp_resid, vcv, prep)$Lambda
      }
    }
    
    idx.min <- which.min(abs(stat.vector))
    PE.sel <- c(PE.sel, idx.min)
    
    if (verbose) message(sprintf(
      "Testing PE #%d (original index: %d)",
      length(PE.sel), positive_idx[idx.min]
    ))
    
    # Update residuals using faster fitting
    fit <- .lm.fit(design_matrices[[idx.min]], current_x)
    current_x <- current_x - design_matrices[[idx.min]] %*% fit$coefficients
    names(current_x) <- names(x)
    
    # Post-test with optimized functions
    if (method == "abouheif") {
      if (use_parallel_moran && nperm > 500) {
        test_res <- phylo.moran.test_parallel(current_x, prox_sparse, method = abouheif_method, nperm = nperm)
      } else {
        test_res <- phylo.moran.test_fast(current_x, prox_sparse, method = abouheif_method, nperm = nperm)
      }
    } else {
      tmp <- lambdaTest_fast(current_x, vcv, prep)
      test_res <- list(obs = tmp$Lambda, pvalue = tmp$pvalue)
    }
    p <- test_res$pvalue
    min.stat <- c(min.stat, test_res$obs)
    p.vector <- c(p.vector, p)
    
    if (verbose) message(sprintf(" → Test obs = %.4f, p = %.4f", test_res$obs, p))
  }
  
  if (verbose) message(sprintf("Procedure stopped (p = %.4f >= alpha = %.2f)", p, alpha))
  
  doParallel::stopImplicitCluster()
  
  selected_PEM <- PEM_positive[, PE.sel, drop = FALSE]
  original_indices <- positive_idx[PE.sel]
  
  stat_label <- if (method == "abouheif") paste0("Iresid_", abouheif_method) else "lambda"
  summary_df <- data.frame(
    variables = paste0("ME", original_indices),
    order = original_indices,
    moran_eigenvalue = positive_moran[PE.sel],
    stat = min.stat,
    pvalue = p.vector,
    row.names = seq_along(PE.sel)
  )
  colnames(summary_df)[4] <- stat_label
  
  res$PE.select <- selected_PEM
  res$summary <- summary_df
  return(res)
}