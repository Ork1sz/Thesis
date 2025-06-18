#' @importFrom methods new
create.block <- function(feature.names, selected.features, block.columns) {
  if (length(feature.names) > 0) {
    selected.features <- feature.names[selected.features]
  }
  return(new("block", features = selected.features, block.columns = block.columns))
}

get.blocks <- function(threshold.matrix, feature.names) {
  p <- nrow(threshold.matrix)
  k <- ncol(threshold.matrix)
  columns <- 1:k
  blocks <- list()
  
  if (any(colSums(threshold.matrix == 0) == p)) {
    stop("Zero column. Reduce threshold value.\n")
  }
  
  if (any(rowSums(threshold.matrix == 0) == k)) {
    if (k < p) {
      stop("Zero row. Add more loadings, or reduce threshold value.\n")
    } else {
      stop("Zero row. Reduce threshold value.\n")
    }
    
  }
  
  while (length(columns) != 0) {
    block.columns <- columns[1]
    while (!identical(block.columns, find.blocks(threshold.matrix, block.columns))) {
      block.columns <- find.blocks(threshold.matrix, block.columns)
    }
    
    if (length(block.columns) == 1) {
      block.idx <- which(threshold.matrix[, block.columns] != 0)
    } else {
      block.idx <- which(rowSums(threshold.matrix[, block.columns] != 0) > 0)
    }
    
    blocks[[length(blocks) + 1]] <- create.block(
      feature.names = feature.names,
      selected.features = block.idx,
      block.columns = block.columns
    )
    
    columns <- columns[!columns %in% block.columns]
  }
  
  return(blocks)
}

get.threshold.matrix <- function(loadings, threshold) {
  loadings[which(abs(loadings) <= threshold)] <- 0
  loadings[which(abs(loadings) != 0)] <- 1
  
  return(loadings)
}

find.blocks <- function(matrix, column.idx) {
  if (length(column.idx) == 1) {
    rows <- which(matrix[, column.idx] != 0)
  } else {
    rows <- which(rowSums(matrix[, column.idx] != 0) > 0)
  }
  
  if (length(rows) == 1) {
    cols <- which(matrix[rows, ] != 0)
  } else {
    cols <- which(colSums(matrix[rows, ] != 0) > 0)
  }
  
  return(cols)
}


#' @title Block
#'
#' @description Class used within the package to store the structure and
#' information about the detected blocks.
#' @slot features numeric vector that contains the the variables
#' corresponding to this block.
#' @slot block.columns numeric vector that contains the indices of the
#' singular vectors corresponding to this block.
#' @export
setClass("block", slots = c(features = "vector", block.columns = "vector"))



#' @title Block Detection Using Singular Vectors (BD-SVD).
#'
#' @description Performs BD-SVD iteratively to reveal the block structure. Splits the data matrix into one (i.e., no split)
#' or two submatrices, depending on the structure of the first sparse loading \eqn{v} (which is a sparse approximation of the
#' first right singular vector, i.e., a vector with many zero values) that mirrors the shape of the covariance matrix. This
#' procedure is continued iteratively until the block diagonal structure has been revealed.
#' 
#' An extended version of BD-SVD can use multiple sparse vectors (`q > 1`). This enables two key features:
#' 1.  **Automatic q-Selection**: When `q.selection = TRUE`, the algorithm automatically
#'     determines the number of vectors to use at each split. The `mode`
#'     argument specifies the selection method, which includes the Eigenvalue Ratio method
#'     (Ahn & Horenstein, 2013), the Empirical Kaiser Criterion (Braeken & van Assen, 2017),
#'     and Parallel Analysis (Horn, 1965).
#' 2.  **Weighted Sparsity Allocation**: For **any** split that uses multiple vectors
#'     (`q > 1`), the algorithm performs a weighted allocation of sparsity. You **must**
#'     specify the allocation strategy using the `weight_method` argument.
#'
#' The data matrix ordered according to this revealed block diagonal structure can be obtained by \code{\link{bdsvd.structure}}.
#'
#' @param X Data matrix of dimension \eqn{n}x\eqn{p} with possibly \eqn{p >> n}.
#'
#' @param q Initial guess for the number of singular vectors. This is ignored if
#'   `q.selection = TRUE`. If `q.selection = FALSE`, manually setting `q > 1`
#'   is ignored with a warning, and `q` is reset to `1`.
#'   
#' @param q.selection Logical; if `TRUE`, automatically selects the number of
#'   singular vectors (`q`) at each split.
#'   
#' @param mode Character; the method for automatic q-selection ("ERM", "EKC", "PAR").
#'   Required if `q.selection = TRUE`.
#'   
#' @param weight_method Character; the method for allocating sparsity when `q > 1`.
#'   Must be one of "equal",  or "explained_var".
#'   
#' @param dof.lim Interval limits for the number of non-zero components in the sparse loading (degrees of freedom).
#' If \eqn{S} denotes the support of \eqn{v}, then the cardinality of the support, \eqn{|S|},
#' corresponds to the degrees of freedom. Default is \code{dof.lim <- c(0, p-1)} which is highly recommended to check for
#' all levels of sparsity.
#'
#' @param anp Which regularization function should be used for the HBIC. \code{anp = "1"} implements \eqn{a_{np} = 1} which corresponds
#' to the BIC, \code{anp = "2"} implements \eqn{a_{np} = 1/2 log(np)} which corresponds to the regularization used by Bauer (2024), and \code{anp = "3"}
#' implements \eqn{a_{np} = log(log(np))} which corresponds to the regularization used by Wang et al. (2009) and Wang et al. (2013).
#'
#' @param standardize Standardize the data to have unit variance. Default is \code{TRUE}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loading.
#' Default is \code{200}.
#'
#' @param trace Print out progress as iterations are performed. Default is \code{TRUE}.
#'
#' @details
#' The sparse loadings are computed using the method by Shen & Huang (2008), implemented by Baglama, Reichel, and Lewis in \code{\link[irlba]{ssvd}} \{\link[irlba]{irlba}\}.
#'
#' @return
#' A list containing the feature names of the submatrices of \code{X}. The length of the list equals
#' the number of submatrices.
#'
#' @seealso \code{\link{bdsvd.structure}}, \code{\link{bdsvd.ht}}, \code{\link{single.bdsvd}}
#'
#' @references
#' \cite{Ahn, S. C., & Horenstein, A. R. (2013). Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203-1227.}
#' \cite{Bauer, J.O. (2024). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat.}
#' \cite{Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion. Psychological methods, 22(3), 450.}
#' \cite{Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 30, 179-185.}
#' \cite{Wang, H., B. Li, and C. Leng (2009). Shrinkage tuning parameter selection with a diverging number of parameters, J. R. Stat. Soc. B 71 (3), 671–683.}
#' \cite{Wang, L., Y. Kim, and R. Li (2013). Calibrating nonconvex penalized regression in ultra-high dimension, Ann. Stat. 41 (5), 2505–2536.}
#'
#' @examples
#' #Replicate the simulation study (c) from Bauer (2024).
#'
#' \dontrun{
#' p <- 500 #Number of variables
#' n <- 500 #Number of observations
#' b <- 10  #Number of blocks
#' design <- "c" #Simulation design "a", "b", "c", or "d".
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#' colnames(X) <- seq_len(p)
#'
#' bdsvd(X, standardize = FALSE)
#' }
#'
#' @importFrom irlba ssvd
#'
#' @export
bdsvd <- function(X,
                  q = 1, 
                  q.selection = FALSE,
                  mode = NULL,# number of singular vectors
                  dof.lim,
                  anp = "2",
                  standardize = TRUE,
                  max.iter,
                  trace = FALSE,
                  weight_method = NULL) {
  
  ## 1) Input Validation
  
  # Ensure q is a positive integer.
  if (!is.numeric(q) || q < 1) stop("`q` must be a positive integer.")
  
  # Handle case where user provides q > 1 withouth automatic selection
  # This warns the user and defaults to the standard single-vector BD-SVD.
  if (!q.selection && q != 1) {
    warning(
      "Manual setting of q > 1 is ignored when `q.selection` is FALSE.\n",
      "Resetting q to 1 and will run single vector BD-SVD for splits.",
      call. = FALSE
    )
    q <- 1
  }
  
  # Check for missing values in the data matrix.
  if (anyNA(X)) stop("X contains missing values (NA).")
  
  # Get the number of columns.
  p <- ncol(X)
  
  # Check if q exceeds the number of columns.
  if (q > p) stop("`q` cannot exceed the number of columns p.")
  
  
  # If q-selection is enabled, a method ('mode') must be provided.
  if (!is.character(mode) && q.selection) {
    stop("`mode` must be a character string (e.g. \"EKC\" or \"PAR\") when `q.selection = TRUE`.")
  }
  
  # If the data matrix has no column names, assign them.
  if (length(colnames(X)) == 0) {
    colnames(X) <- as.character(seq_len(p))
  }
  
  # Set default values for the degrees of freedom limits if not provided.
  if (missing(dof.lim)) dof.lim <- c(0, p - 1)
  dof.lim <- sort(round(dof.lim))
  dof.lim[1] <- max(dof.lim[1], 0)
  dof.lim[2] <- min(dof.lim[2], p - 1)
  
  # Validate the choice for the HBIC penalty.
  ANP <- c("1","2","3")
  if (!anp %in% ANP) stop("`anp` must be one of ", paste(ANP, collapse = ", "))
  
  # Set a default for the maximum number of iterations.
  if (missing(max.iter)) {
    max.iter <- 500
  }
  
  
  # --- 2) Initialize the Algorithm ---
  
  # Create a list to hold sub-matrices that need to be split.
  sub.matrices <- list(colnames(X))
  results <- list()
  
  # Counter for discovered blocks.
  b <- 1 
  
  # Vector to store the q value used for each split.
  chosen_qs     <- integer(0)
  
  
  # --- 3) Main Loop For Block Detection ---
  # Continue processing as long as there are sub-matrices in the queue.
  while (length(sub.matrices) > 0) {
    # Get the next sub-matrix (by column names) from the front of the queue.
    cols <- sub.matrices[[1]]
    sub.matrices <- sub.matrices[-1]
    
    # if only one column, it's a final block
    if (length(cols) == 1) {
      if (trace) cat("Block", b, ":", cols, "\n")
      results[[b]] <- cols
      b <- b + 1
      next
    }
    
    # Extract the data for the current sub-matrix.
    Xi <- X[, cols, drop = FALSE]
    
    # If q-selection is enabled, call pick.q to determine the number of singular vectors.
    if (q.selection) {
      q <- pick.q(
        X = Xi,
        method = mode
      )
    }
    
    # Ensure q is not larger than the number of columns in the current sub-matrix.
    q <- min(q, ncol(Xi))
    
    # Record the value of q used for this split.
    chosen_qs <- c(chosen_qs, q)
    
    
    if(trace){
      cat(">> Chosen q = ", q , "\n")
    }
    
    # Find the optimal degrees of freedom (dof) for the split.
    dof.opt <- if (q == 1) {
      # For a single vector, use the standard hyperparameter tuning.
      bdsvd.ht(X = Xi, dof.lim = dof.lim, anp = anp,
               standardize = standardize, max.iter = max.iter)$dof
      
    } else {
      # For multiple vectors, use multivariate hyperparameter tunning.
      multi.bdsvd.ht(X = Xi, q = q, dof.lim = dof.lim, anp = anp,
                     standardize = standardize, max.iter = max.iter, weight_method = weight_method)$dof
    }
    
    if(trace){    cat(">> Chosen support for this submatrix (|\n Xi | = ", ncol(Xi), "): ",ncol(Xi) - dof.opt, "\n")
    }
    # Run the appropriate splitting function based on the value of q.
    sub.results <- if (q == 1) {
      single.bdsvd(X = Xi, dof = dof.opt,
                   standardize = standardize, max.iter = max.iter)
    } else {
      # For q > 1, a weight_method must be provided to allocate sparsity.
      multi.single.bdsvd.weighted(X = Xi, q = q, dof = dof.opt,
                                  standardize = standardize, max.iter = max.iter, weight_method = weight_method)
      
    }
    
    # collect or queue up the subblocks
    if (length(sub.results) == 1) {
      if (trace) cat("Block", b, ":", sub.results[[1]], "\n")
      results[[b]] <- sub.results[[1]]
      b <- b + 1
    } else {
      sub.matrices <- c(sub.matrices, sub.results)
    }
  }
  
  # sort final blocks by size 
  results <- results[order(sapply(results, length), decreasing = TRUE)]
  class(results) <- "bdsvd"
  
  # Attach history of chosen q's
  attr(results, "chosen_qs") <- chosen_qs
  
  return(results)
}



#' @title Covariance Matrix Simulation for BD-SVD
#'
#' @description This function generates covariance matrices based on the simulation studies described in Bauer (2024).
#'
#' @param p Number of variables.
#'
#' @param b Number of blocks. Only required for simulation design "c" and "d".
#'
#' @param design Simulation design "a", "b", "c", "d" or the new "e".
#'
#' @return
#' A covariance matrix according to the chosen simulation design.
#'
#' @references \cite{Bauer, J.O. (2024). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat.}
#'
#' @examples
#' #The covariance matrix for simulation design (a) is given by
#' Sigma <- bdsvd.cov.sim(p = 500, b = 500, design = "a")
#' 
#' The covariance matrix for the noisy block design (e)
#' Sigma_e <- bdsvd.cov.sim(p = 500, b = 10, design = "e", noise_level = 0.1)
#'
#' @importFrom stats runif
#' @importFrom stats runif rnorm
#' @importFrom Matrix nearPD
#' @export
bdsvd.cov.sim <- function(p = p,
                          b,
                          design = design,
                          noise_level = 0.05
) {
  
  # --- Input Validation ---
  # Define the valid design options.
  DESIGN <- c("a", "b", "c", "d","e")
  if (!(design %in% DESIGN))
    stop(paste(design), " is an invalid design")
  
  
  #--- Design "a": Identity Matrix ---
  # This represents uncorrelated variables with unit variance.
  if (design == "a") {
    Sigma <- diag(1, p, p)
    return(Sigma)
  }
  
  # --- Design "b": Diagonal Matrix with Random Variances ---
  # This represents uncorrelated variables with heterogeneous variances.
  if (design == "b") {
    Sigma <- diag(stats::runif(p, 1, 5), p, p)
    return(Sigma)
  }
  
  # --- Design "c": Block-Diagonal with Equi-correlation ---
  # Creates a matrix with `b` equally sized blocks, where each block has the same correlation `rho`.
  if (design == "c") {
    if (missing(b))
      stop("Number of blocks b required for simulation design c.")
    
    if (!(p %% b == 0))
      stop("p must be divisible by b so that blocks of equal size can be created.")
    
    rho0 <- 0.2
    epsilon <- 0.1
    
    d <- p / b
    Sigma <- matrix(0, p, p)
    for (i in 1:b) {
      rho <- stats::runif(1, rho0 - epsilon, rho0 + epsilon)
      Sigma[(1 + d * (i - 1)) : (d + d * (i - 1)), (1 + d * (i - 1)) : (d + d * (i - 1))] <-
        ((1 - rho) * diag(1, d, d) + 2 * rho * rep(1, d) %*% t(rep(1, d)))
    }
    return(Sigma)
  }
  
  # --- Design "d": Block-Diagonal with Autoregressive Structure ---
  # Creates a matrix with `b` blocks, where each block has a toeplitz-like structure.
  if (design == "d") {
    if (missing(b))
      stop("Number of blocks b required for simulation design d.")
    if (!(p %% b == 0))
      stop("p must be divisible by b so that blocks of equal size can be created.")
    
    rho0 <- 0.45
    epsilon <- 0.15
    omega <- 0.1
    
    d <- p / b
    Sigma <- matrix(0, p, p)
    for (B in 1:b) {
      rho <- stats::runif(1, rho0 - epsilon, rho0 + epsilon)
      Rii <- matrix(0, d, d)
      for (i in 1:d) {
        for (j in 1:d) {
          Rii[i, j] <- (-1)^(i + j) * rho^(abs(i - j)^(omega))
        }
      }
      diag(Rii) <- 1
      Sigma[(1 + d * (B - 1)) : (d + d * (B - 1)), (1 + d * (B - 1)):(d + d * (B - 1))] <- Rii
    }
    return(Sigma)
  }
  
  # --- Design "e": Leaky-Block Matrix (Block-Diagonal + Noise) ---
  # Creates a noisy block-diagonal matrix that is still positive definite.
  if (design == "e") {
    if (missing(b))
      stop("Number of blocks b required for simulation design e.")
    
    # 1. Start with a base block-diagonal matrix (using the logic from design "c")
    Sigma_block <- bdsvd.cov.sim(p = p, b = b, design = "c")
    
    # 2. Create a symmetric noise matrix
    E <- matrix(rnorm(p * p, mean = 0, sd = noise_level), p, p)
    E <- (E + t(E)) / 2 # Ensure symmetry
    diag(E) <- 0         # No noise on the diagonal itself
    
    # 3. Add the noise to the block matrix
    Sigma_leaky <- Sigma_block + E
    
    # 4. Ensure the final matrix is positive definite and return
    # We set corr=TRUE if we want the diagonal to be exactly 1.
    Sigma_final <- as.matrix(Matrix::nearPD(Sigma_leaky, corr = TRUE)$mat)
    
    return(Sigma_final)
  }
  
}



#' @title Hyperparameter Tuning for BD-SVD
#'
#' @description Finds the number of non-zero elements of the sparse loading according to the high-dimensional
#' Bayesian information criterion (HBIC).
#'
#' @param X Data matrix of dimension \eqn{n x p} with possibly \eqn{p >> n}.
#'
#' @param dof.lim Interval limits for the number of non-zero components in the sparse loading (degrees of freedom).
#' If \eqn{S} denotes the support of \eqn{v}, then the cardinality of the support, \eqn{|S|},
#' corresponds to the degrees of freedom. Default is \code{dof.lim <- c(0, p-1)} which is highly recommended to check for
#' all levels of sparsity.
#'
#' @param anp Which regularization function should be used for the HBIC. \code{anp = "1"} implements \eqn{a_{np} = 1} which corresponds
#' to the BIC, \code{anp = "2"} implements \eqn{a_{np} = 1/2 log(np)} which corresponds to the regularization used by Bauer (2024), and \code{anp = "3"}
#' implements \eqn{a_{np} = log(log(np))} which corresponds to the regularization used by Wang et al. (2009) and Wang et al. (2013).
#'
#' @param standardize Standardize the data to have unit variance. Default is \code{TRUE}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loading.
#' Default is \code{200}.
#'
#' @details
#' The sparse loadings are computed using the method by Shen & Huang (2008), implemented in
#' the \code{irlba} package. The computation of the HBIC is outlined in Bauer (2024).
#'
#' @return
#' \item{dof}{
#'   The optimal number of nonzero components (degrees of freedom) according to the HBIC.
#' }
#' \item{BIC}{
#'   The HBIC for the different numbers of nonzero components.
#' }
#'
#' @seealso \code{\link{bdsvd}}, \code{\link{single.bdsvd}}
#'
#' @references \cite{Bauer, J.O. (2024). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat.}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#' @references \cite{Wang, H., B. Li, and C. Leng (2009). Shrinkage tuning parameter selection with a diverging number of parameters, J. R. Stat. Soc. B 71 (3), 671–683.}
#' @references \cite{Wang, L., Y. Kim, and R. Li (2013). Calibrating nonconvex penalized regression in ultra-high dimension, Ann. Stat. 41 (5), 2505–2536.}
#'
#' @examples
#' #Replicate the illustrative example from Bauer (2024).
#'
#'
#' p <- 300 #Number of variables. In Bauer (2024), p = 3000
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#' colnames(X) <- seq_len(p)
#'
#' ht <- bdsvd.ht(X)
#' plot(0:(p-1), ht$BIC[,1], xlab = "|S|", ylab = "HBIC", main = "", type = "l")
#' single.bdsvd(X, dof = ht$dof, standardize = FALSE)
#'
#' @importFrom irlba ssvd
#'
#' @export
bdsvd.ht <- function(X,
                     dof.lim,
                     standardize = TRUE,
                     anp = "2",
                     max.iter
) {
  
  if (anyNA(X)) {
    stop("X contains missing value indicator (NA)")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (missing(dof.lim)) {
    dof.lim <- c(0, p - 1)
  }
  dof.lim <- sort(round(dof.lim))
  if (dof.lim[1] < 0) {
    dof.lim[1] <- 0
  }
  if (dof.lim[2] > p - 1) {
    dof.lim[2] <- p - 1
  }
  dof.grid <- dof.lim[1]:dof.lim[2]
  
  X <- scale(X, center = TRUE, scale = standardize)
  
  ANP <- c("1", "2", "3")
  if (!(anp %in% ANP))
    stop(paste(anp), " is an invalid option for anp")
  
  if (anp == "2") {
    a_np <- function(n, p) {
      1 / 2 * log(n * p)
    }
  } else if (anp == "1") {
    a_np <- function(n, p) {
      1
    }
  } else if (anp == "3") {
    a_np <- function(n, p) {
      log(log(n * p))
    }
  }
  
  if (missing(max.iter)) {
    max.iter <- 500
  }
  
  BIC <- vector(length = length(dof.grid))
  i <- 1
  for (dof in dof.grid) {
    v <- tryCatch(suppressWarnings(irlba::ssvd(x = X, k = 1, n = p - dof, maxit = max.iter)$v), error = function(e) e)
    if (inherits(v, "error")) {
      v <- matrix(0, nrow = p, ncol = 1)
    }
    
    u <- X %*% v
    BIC[i] <- log(norm(X - u %*% t(v), type = "F")^2 / n / p) + sum(v != 0) * log(n * p) / n / p *  a_np(n, p)
    i <- i + 1
  }
  
  dof <- dof.grid[order(BIC)[1]]
  BIC <- cbind.data.frame(BIC)
  rownames(BIC) <- dof.grid
  
  return(list(dof = dof, BIC = BIC))
}



#' @title Data Matrix Structure According to the Detected Block Structure.
#'
#' @description Either sorts the data matrix \eqn{X} according to the detected block structure \eqn{X_1 , ... , X_b}, ordered by the number
#' of variables that the blocks contain. Or returns the detected submatrices each individually in a list object.
#'
#' @param X Data matrix of dimension \eqn{n x p} with possibly \eqn{p >> n}.
#'
#' @param block.structure Output of \code{bdsvd()} or \code{single.bdsvd()} which identified the block structure.
#'
#' @param output Should the output be the data matrix ordered according to the blocks (\code{"matrix"}), or
#' a list containing the submatrices (\code{"submatrices"}). Default is \code{"matrix"}.
#'
#' @param block.order A vector that contains the order of the blocks detected by \code{bdsvd()} or \code{single.bdsvd()}.
#' The vector must contain the index of each blocks exactly once. Default is \code{1:b} where \code{b} is the total number of blocks.
#'
#' @return
#' Either the data matrix \code{X} with columns sorted according to the detected blocks, or a list containing the detected
#' submatrices.
#'
#' @seealso \code{\link{bdsvd}}, \code{\link{single.bdsvd}}
#'
#' @references \cite{Bauer, J.O. (2024). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat.}
#'
#' @examples
#' #Toying with the illustrative example from Bauer (2024).
#'
#'
#' p <- 150 #Number of variables. In Bauer (2024), p = 3000.
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#' colnames(X) <- seq_len(p)
#'
#' #Compute iterative BD-SVD
#' bdsvd.obj <- bdsvd(X, standardize = FALSE)
#'
#' #Obtain the data matrix X, sorted by the detected blocks
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "matrix") )
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "matrix", block.order = c(2,1,3)) )
#'
#' #Obtain the detected submatrices X_1, X_2, and X_3
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "submatrices")[[1]] )
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "submatrices")[[2]] )
#' colnames(bdsvd.structure(X, bdsvd.obj, output = "submatrices")[[3]] )
#'
#' @export
bdsvd.structure <- function(X,
                            block.structure,
                            output = "matrix",
                            block.order
) {
  
  if (!inherits(block.structure, "bdsvd"))
    stop("block.structure must be the outcome of bdsvd() or single.bdsvd().")
  
  OUTPUT <- c("matrix", "submatrices")
  if (!(output %in% OUTPUT))
    stop(paste(output), " is an invalid argument for output")
  
  p <- ncol(X)
  if (length(colnames(X)) == 0) {
    colnames(X) <- as.character(1:p)
  }
  
  b <- length(block.structure)
  if (b == p)
    return(result)
  
  ifelse(missing(block.order),
         block.order <- seq_along(block.structure),
         block.order <- as.integer(block.order))
  if (!identical(sort(block.order), 1:b)) {
    stop("block.order must contain the index of each block exactly once.")
  }
  
  if (output == "matrix") {
    result <- X[, colnames(X) %in% block.structure[[block.order[1]]], drop = FALSE]
    for (i in block.order[-1]) {
      result <- cbind.data.frame(result,  X[, colnames(X) %in% block.structure[[i]], drop = FALSE])
    }
  } else {
    result <- list()
    for (i in block.order) {
      result <- c(result, list(X[, colnames(X) %in% block.structure[[block.order[i]]], drop = FALSE]))
    }
  }
  
  return(result)
}



#' @title Block Detection
#'
#' @description This function returns the block structure of a matrix.
#'
#' @param V Numeric matrix which either contains the loadings or is a covariance matrix.
#'
#' @param threshold All absolute values of \code{V} below the threshold are set to zero.
#'
#' @return
#' An object of class \code{Block} containing the features and columns indices corresponding to each detected block.
#'
#' @seealso \code{\link{bdsvd}}, \code{\link{single.bdsvd}}
#'
#' @references \cite{Bauer, J.O. (2024). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat.}
#'
#' @examples
#' #In the first example, we replicate the simulation study for the ad hoc procedure
#' #Est_0.1 from Bauer (2024). In the second example, we manually compute the first step
#' #of BD-SVD, which can be done using the bdsvd() and/or single.bdsvd(), for constructed
#' #sparse loadings
#'
#' #Example 1: Replicate the simulation study (a) from Bauer (2024) for the ad hoc
#' #procedure Est_0.1.
#'
#'\dontrun{
#' p <- 500 #Number of variables
#' n <- 125 #Number of observations
#' b <- 500 #Number of blocks
#' design <- "a"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=Sigma)
#' colnames(X) <- 1:p
#'
#' #Perform the ad hoc procedure
#' detect.blocks(cvCovEst::scadEst(dat = X, lambda = 0.2), threshold = 0)
#' }
#'
#' #Example 2: Manually compute the first step of BD-SVD
#' #for some loadings V that mirror the two blocks
#' #("A", "B") and c("C", "D").
#'
#' V <- matrix(c(1,0,
#'               1,0,
#'               0,1,
#'               0,1), 4, 2, byrow = TRUE)
#'
#' rownames(V) <- c("A", "B", "C", "D")
#' detected.blocks <- detect.blocks(V)
#'
#' #Variables in block one with corresponding column index:
#' detected.blocks[[1]]@features
#' detected.blocks[[1]]@block.columns
#'
#' #Variables in block two with corresponding column index:
#' detected.blocks[[2]]@features
#' detected.blocks[[2]]@block.columns
#'
#' @export
detect.blocks <- function(V,
                          threshold = 0
) {
  
  if (missing(V)) {
    stop("V is required.")
  }
  
  if (length(rownames(V)) == 0) {
    rownames(V) == as.character(seq_len(nrow(V)))
  }
  
  threshold.matrx <- get.threshold.matrix(V, threshold)
  
  return(get.blocks(threshold.matrx, rownames(V)))
}



#' @title Single Iteration of Block Detection Using Singular Vectors (BD-SVD).
#'
#' @description Performs a single iteration of BD-SVD: splits the data matrix into one (i.e., no split)
#' or two submatrices, depending on the structure of the first sparse loading \eqn{v} (which is a sparse
#' approximation of the first right singular vector, i.e., a vector with many zero values) that mirrors the
#' shape of the covariance matrix.
#'
#' @param X Data matrix of dimension \eqn{n x p} with possibly \eqn{p >> n}.
#'
#' @param dof Number of non-zero components in the sparse loading (degrees of freedom). If
#' \eqn{S} denotes the support of \eqn{v}, then the cardinality of the support, \eqn{|S|},
#' corresponds to the degrees of freedom.
#'
#' @param standardize Standardize the data to have unit variance. Default is \code{TRUE}.
#'
#' @param max.iter How many iterations should be performed for computing the sparse loading.
#' Default is \code{200}.
#'
#' @details
#' The sparse loadings are computed using the method by Shen & Huang (2008), implemented in
#' the \code{irlba} package.
#'
#' @return
#' A list containing the feature names of the submatrices of \code{X}. It is either of length one (no
#' split) or length two (split into two submatrices).
#'
#' @seealso \code{\link{bdsvd}}, \code{\link{bdsvd.ht}}
#'
#' @references \cite{Bauer, J.O. (2024). High-dimensional block diagonal covariance structure detection using singular vectors, J. Comput. Graph. Stat.}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#'
#' @examples
#' #Replicate the illustrative example from Bauer (2024).
#'
#' \dontrun{
#'
#' p <- 300 #Number of variables. In Bauer (2024), p = 3000.
#' n <- 500 #Number of observations
#' b <- 3   #Number of blocks
#' design <- "c"
#'
#' #Simulate data matrix X
#' set.seed(1)
#' Sigma <- bdsvd.cov.sim(p = p, b = b, design = design)
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#' colnames(X) <- 1:p
#'
#' ht <- bdsvd.ht(X)
#' plot(0:(p-1), ht$BIC[,1], xlab = "|S|", ylab = "HBIC", main = "", type = "l")
#' single.bdsvd(X, dof = ht$dof, standardize = FALSE)
#'
#' }
#'
#' @importFrom irlba ssvd
#'
#' @export
single.bdsvd <- function(X,
                         dof,
                         standardize = TRUE,
                         max.iter
) {
  
  if (anyNA(X)) {
    stop("X contains missing value indicator (NA)")
  }
  
  p <- ncol(X)
  
  if (missing(dof)) {
    stop("Enter the degrees of freedom")
  }
  
  if (missing(max.iter)) {
    max.iter <- 500
  }
  
  X <- scale(X, center = TRUE, scale = standardize)
  if (length(colnames(X)) == 0) {
    colnames(X) <- as.character(seq_len(p))
  }
  
  if (length(unique(colnames(X))) != p) {
    stop("Variable names are not unique")
  }
  
  feature.names <- colnames(X)
  eigen <- list()
  
  v <- tryCatch(suppressWarnings(irlba::ssvd(x = X, k = 1, n = p - dof, maxit = max.iter)$v), error = function(e) e)
  if (inherits(v, "error")) {
    stop("dof = ", dof, " do not fit the structure of the singular vectors. You can use bdsvd.ht to find a suitable value for dof.")
  }
  
  if (length(which(v == 0)) == 0) {
    return(list(feature.names))
  }
  
  eigen$vectors <- cbind(v, 0)
  eigen$vectors[which(v == 0), 2] <- 1
  rownames(eigen$vectors) <- feature.names
  
  
  detected.blocks <- detect.blocks(eigen$vectors, 0)
  result <- list()
  result[[1]] <- detected.blocks[[1]]@features
  result[[2]] <- detected.blocks[[2]]@features
  class(result) <- "bdsvd"
  
  return(result)
}


### Block detection for BDSVD q >1 

#' @title Threshold Matrix Conversion
#' @description
#' Turn a p×q loading matrix V into a 0/1 mask using a cutoff.
#'
#' @param loadings Numeric p×q matrix of loadings.
#' @param threshold Numeric cutoff: entries ≤ threshold become 0.
#'
#' @return Integer p×q mask of 0s and 1s.
#' @export
multi.get.threshold.matrix <- function(loadings, threshold) {
  
  # --- Input Validation ---
  # Ensure the 'loadings' input is a matrix.
  if (!is.matrix(loadings)) stop("`loadings` must be a matrix.")
  # Ensure the 'threshold' is a single numeric value.
  if (!is.numeric(threshold) || length(threshold) != 1) stop("`threshold` must be one number.")
  
  # --- Create Binary Mask ---
  # Create a logical matrix where TRUE corresponds to elements whose
  # absolute value is greater than the threshold.
  mask <- abs(loadings) > threshold
  # Convert the logical matrix (TRUE/FALSE) into an integer matrix (1/0).
  res <- matrix(as.integer(mask), nrow = nrow(loadings), ncol = ncol(loadings))
  # Preserve the original row and column names from the loadings matrix.
  dimnames(res) <- dimnames(loadings)
  # Return final binary mask
  return(res)
}

#' @title Build Multivariate Blocks
#' @description
#' Assigns each feature to the first loading that selects it (a "greedy" approach),
#' then creates one block for each loading's claimed features, plus a final
#' block for any features that remain unselected.
#'
#' @param mask An integer p x q matrix of 0s and 1s, where mask[i, j] = 1
#'   indicates that feature i is selected by loading j.
#' @param feature.names A character vector of feature names.
#'
#' @return List of "block" objects for each loading and the remainder\.
#' @export
multi.get.blocks <- function(mask, feature.names) {
  # Get the dimensions of the feature mask.
  p <- nrow(mask)
  q <- ncol(mask)
  
  # Create a vector to track which loading has "claimed" each feature.
  # It's initialized to NA, meaning no feature is assigned yet.
  assigned <- rep(NA_integer_, p)
  
  # Greedy assignment: first loading to pick each feature keeps it
  for (j in seq_len(q)) {
    rows_j <- which(mask[, j] == 1)
    new_rows <- rows_j[is.na(assigned[rows_j])]
    assigned[new_rows] <- j
  }
  
  # Initialize an empty list to store the final block objects.
  blocks <- list()
  
  # Build a block for each loading
  for (j in seq_len(q)) {
    # Get the indices of all features that were assigned to loading 'j'.
    idx <- which(assigned == j)
    
    # If this loading claimed any features, create a block for them.
    if (length(idx) > 0) {
      blocks[[length(blocks) + 1]] <- create.block(
        feature.names = feature.names,
        selected.features = idx,
        block.columns = j
      )
    }
  }
  
  # Find any features that were not selected by ANY of the loadings.
  rem_idx <- which(is.na(assigned))
  
  # If there are unassigned features, group them into one final "remainder" block.
  if (length(rem_idx) > 0) {
    blocks[[length(blocks) + 1]] <- create.block(
      feature.names = feature.names,
      selected.features = rem_idx,
      block.columns = NA # This block doesn't correspond to a specific loading.
    )
  }
  
  # Return the complete list of block objects.
  return(blocks)
}



#' @title Hyperparameter Tuning for Multivariate BD-SVD (Weighted)
#'
#' @description
#' Finds the optimal sparsity (degrees of freedom) for `q` sparse loadings by
#' minimizing a multivariate HBIC. This version incorporates weighted sparsity
#' allocation directly into the tuning process for more consistent results.
#'
#' @param X Data matrix (n x p).
#' @param q Number of singular vectors.
#' @param weight_method Method for allocating sparsity ('equal' or'explained_var').
#' @param dof.lim Numeric vector of length 2 for min and max non-zeros per loading.
#' @param standardize Logical; scale columns if TRUE.
#' @param anp Penalty type: "1", "2", "3".
#' @param max.iter Max IRLBA iterations.
#'
#' @return List with `dof` (optimal degrees of freedom) and `BIC` (a data.frame of all tested values).
#' @export
multi.bdsvd.ht <- function(X, q, weight_method ,dof.lim, standardize = TRUE, anp = "2", max.iter = 500) {
  # --- 1. Input Validation and Setup ---
  if (anyNA(X)) stop("X contains NA values.")
  n <- nrow(X); p <- ncol(X)
  if (!is.numeric(q) || q < 1 || q > min(n,p)) stop("`q` must be between 1 and min(n,p)")
  
  # Define the grid of possible degrees of freedom (dof) to test.
  if (missing(dof.lim)) dof.lim <- c(0, p-1)
  dof.lim <- sort(round(dof.lim))
  dof.lim[1] <- max(dof.lim[1], 0)
  dof.lim[2] <- min(dof.lim[2], p-1)
  
  # Ensure that the number of non-zero elements is at least q.
  max_dof <- p - q
  dof.lim[2] <- min(dof.lim[2], max_dof)
  if (dof.lim[2] < dof.lim[1]) stop("No valid dof: must satisfy p - dof >= q.")
  dof.grid <- seq(dof.lim[1], dof.lim[2])
  
  # Standardize the data matrix if requested.
  X <- scale(X, center = TRUE, scale = standardize)
  
  # Define the penalty term for the HBIC based on the 'anp' parameter.
  ANP <- c("1", "2", "3")
  if (!(anp %in% ANP))
    stop(paste(anp), " is an invalid option for anp")
  
  if (anp == "2") {
    a_np <- function(n, p, q) {
      1 / 2 * log(n * p *q)
    }
  } else if (anp == "1") {
    a_np <- function(n, p, q) {
      1
    }
  } else if (anp == "3") {
    a_np <- function(n, p, q) {
      log(log(n * p * q))
    }
  }
  
  if (missing(max.iter)) {
    max.iter <- 500
  }
  
  # --- 2. HBIC Calculation Loop ---
  # Initialize a vector to store the HBIC value for each dof in the grid.
  HBIC_vals <- numeric(length(dof.grid))
  
  # Iterate over each candidate value for the total degrees of freedom.
  for (i in seq_along(dof.grid)) {
    dof <- dof.grid[i]
    
    # Calculate the total number of non-zero elements to allocate across q vectors.
    total_nonzeros <- p-dof
    
    n_vec <- calculate_allocation_vector(X=X,q=q,total_nonzeros = total_nonzeros,weight_method = weight_method)
    
    
    # Compute the sparse SVD using the calculated allocation vector `n_vec`.
    svd_try <- tryCatch(
      suppressWarnings(irlba::ssvd(x = X, k = q, n = n_vec, maxit = max.iter)),
      error = function(e) e
    )
    
    # Calculate the sum of squared residuals (SSR) for this level of sparsity.
    V <- if (inherits(svd_try, "error")) matrix(0, p, q) else svd_try$v
    U <- X %*% V
    SSR <- norm(X - U %*% t(V), "F")^2
    
    # For the penalty, count the number of non-zero elements, which is more robust.
    df_total <- sum(V != 0)
    
    # Calculate the HBIC value according to the formula in the paper.
    HBIC_vals[i] <- log(SSR / (n * p * q)) + (df_total / (n * p * q)) * a_np(n,p,q) * log(n * p * q)
  }
  # --- 3. Find and Return Optimum ---
  # Find the 'dof' value that resulted in the minimum HBIC score.
  best <- which.min(HBIC_vals)
  BIC_df <- data.frame(HBIC = HBIC_vals, row.names = dof.grid)
  # Return the optimal dof and the full table of BIC values.
  return(list(dof = dof.grid[best], BIC = BIC_df))
}



#' @title Pick Number of Components Using a Specified Method
#'
#' @description
#' A wrapper function that selects the number of singular vectors, q, by calling
#' one of several established statistical methods.
#'
#' @param X Numeric matrix (n x p), with no missing values.
#' @param method Character; the method to use. Must be one of "EKC", "PAR", or "ERM".
#' @return Integer: the estimated number of components, q.
#' @export
pick.q <- function(X,
                   method) 
{
  # --- 1. Input Validation ---
  # Ensure both the data matrix and the method are provided.
  if (missing(X) || missing(method)) {
    stop("`X` and `method` are required.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("`X` must be a matrix or data.frame.")
  }
  if (!is.character(method) || length(method) != 1) {
    stop("`method` must be a single string: \"EKC\", \"ERM\" or \"ERM\".")
  }
  
  # --- 2. Data Preparation and Validation ---
  # Convert data.frame to a matrix if necessary.
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix or data.frame.")
  }
  
  if (anyNA(X)) {
    stop("`X` contains missing values; please provide a complete matrix.")
  }
  
  # Check for valid matrix dimensions.
  n <- nrow(X)
  p <- ncol(X)
  if (n < 2 || p < 1) {
    stop("`X` requires at least 2 rows and 1 column.")
  }
  
  
  # --- 3. Dispatch to Selected Method ---
  # Convert the method name to uppercase to make the matching case-insensitive.
  method <- toupper(method)
  raw_q <- switch(
    method,
    # Call the Empirical Kaiser Criterion function.
    EKC = EKC(X),
    # Call the Parallel Analysis function.
    PAR = parallel(X),
    # Call the Eigenvalue Ratio Method function.
    ERM = ERM(X),
    # If the method is not one of the above, stop with an error.
    stop("Unknown method: ", method)
  )
  
  # --- 4. Enforce Valid Bounds ---
  # Ensure the estimated q is at least 1 and no larger than the smaller
  # dimension of the data matrix (min(n, p)).
  raw_q <- min(max(raw_q, 1L), min(n, p))
  
  # Return the final number of components.
  return(raw_q)
}

#' Empirical Kaiser Criterion (EKC) wrapper
#'
#' @param X A numeric matrix or data.frame of raw data.
#' @return Integer: estimated number of factors.
#' @importFrom EFA.dimensions EMPKC
#' @importFrom stats cor
#' @export
EKC <- function(X, corkind = "pearson", verbose = FALSE) {
  # Ensure the required package is available
  if (!requireNamespace("EFA.dimensions", quietly = TRUE)) {
    stop("Package 'EFA.dimensions' is required. Please install it.")
  }
  
  # 1) Input validation
  if (!is.matrix(X) && !is.data.frame(X)) stop("`X` must be a matrix or data.frame.")
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("`X` must be numeric.")
  if (anyNA(X)) stop("`X` contains missing values.")
  
  # --- Suppress Console Output ---
  os_redirect <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
  sink(os_redirect)
  
  result <- tryCatch({
    # --- Conditional logic to handle the square matrix bug ---
    if (nrow(X) == ncol(X)) {
      
      # WORKAROUND PATH: If data is square (n == p), the package's auto-detection fails.
      # We manually compute the correlation matrix.
      # Note: This workaround forces a Pearson correlation in this specific case.
      R <- cor(X)
      N <- nrow(X)
      
      EFA.dimensions::EMPKC(data = R, Ncases = N, verbose = verbose)
      
    } else {
      
      # STANDARD PATH: If data is not square, pass it directly.
      # We explicitly provide Ncases to be robust against any other edge cases.
      EFA.dimensions::EMPKC(
        data    = X,
        corkind = corkind,
        verbose = verbose
      )
    }
    
  }, finally = {
    # CRUCIAL: Always restore console output
    sink()
  })
  # --- End Suppression ---
  
  # Extract and return the number of factors
  return(result$NfactorsEMPKC)
}

#' @title Parallel Analysis via psych::fa.parallel
#' @description
#' Uses Horn's Parallel Analysis (permutation-based) from the psych package
#' to estimate the number of principal components (signals) in X.
#'
#' @param X Numeric data matrix (n × p), rows = observations, columns = variables.
#' @param n.iter Number of PA iterations (permutations or simulations). Default = 500.
#' @param alpha Significance level for the percentile cutoff. Default = 0.05.
#' @return Integer: estimated number of components (signals) \(\hat{q}\).
#' @examples
#' library(psych)
#' set.seed(123)
#' X <- matrix(rnorm(250*500), nrow=250)
#' q_pa <- parallel(X, n.iter=200, alpha=0.05)
#' @import psych
#' @export
parallel <- function(X, n.iter = 500, alpha = 0.05) {
  # 1) Input checks
  if (!is.matrix(X) && !is.data.frame(X))
    stop("`X` must be a matrix or data.frame.")
  X <- as.matrix(X)
  if (!is.numeric(X))
    stop("`X` must be numeric.")
  if (anyNA(X))
    stop("`X` contains missing values.")
  
  # 2) Horn's Parallel Analysis via psych
  #    fa.parallel returns $ncomp (number of components)
  pa_res <- psych::fa.parallel(
    x          = X,
    n.obs      = nrow(X),
    fa         = "pc",
    n.iter     = n.iter,
    quant      = 1 - alpha,
    show.legend= FALSE,
    plot       = FALSE
  )
  
  # 3) Extract and return the PC count
  return(pa_res$ncomp)
}


#' @title Eigenvalue Ratio Method (ER) for Selecting Number of Signals
#' @description Estimates the number of factors by calling the PCA_FN function from the HDRFA package.
#' @param X A numeric matrix with observations in rows and variables in columns.
#' @param rmax The user-supplied maximum number of factors to test. If NULL (the default),
#' a reasonable value is calculated based on the dimensions of X.
#' @return The estimated number of factors (signals) as a single integer.
#' @references Ahn, S.C. and Horenstein, A.R. (2013). Eigenvalue Ratio Test for the Number of Factors.
#' @importFrom HDRFA PCA_FN
#' @export
ERM <- function(X, rmax = NULL) {
  # --- Input Checks ---
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix.")
  if (anyNA(X)) stop("X contains missing values.")
  
  # --- Check for HDRFA Package ---
  # This ensures the required package is installed before proceeding.
  if (!requireNamespace("HDRFA", quietly = TRUE)) {
    stop("Package 'HDRFA' is required. Please install it using: install.packages('HDRFA')")
  }
  
  # --- Set rmax if not provided ---
  # The PCA_FN function requires rmax. If the user doesn't provide one,
  # we calculate a common default value.
  if (is.null(rmax)) {
    n <- nrow(X)
    p <- ncol(X)
    rmax <- floor(min(n, p) / 2)
  }
  
  # --- Call the Package Function ---
  # Apply the PCA_FN function from the HDRFA package to the data.
  num_factors <- HDRFA::PCA_FN(X = X, rmax = rmax)
  
  return(num_factors)
}


#' @title Single Iteration of Block Detection Using q Weighted Sparse Singular Vectors
#'
#' @description
#' Splits the data matrix using q sparse loadings with weighted sparsity allocation.
#' The total number of non-zero elements (determined by `dof`) is distributed
#' among the q vectors according to the specified `weight_method`.
#' 
#' @param X Data matrix (n x p) to block-split.
#' @param q Number of singular vectors to use.
#' @param dof The total number of zero entries (degrees of freedom) across all q vectors.
#' @param standardize Whether to scale columns to unit variance (default TRUE).
#' @param max.iter Maximum IRLBA iterations.
#' @param weight_method The method for allocating sparsity ('equal' or 'explained_var').
#'
#' @return A list of character vectors, with each vector representing a detected block.
#' @export
multi.single.bdsvd.weighted <- function(X, q, dof, standardize = TRUE, max.iter, weight_method = NULL) {
  # --- 1. Input Validation and Setup ---
  if (anyNA(X)) stop("X contains NA values")
  p <- ncol(X)
  if (missing(dof)) stop("`dof` (degrees of freedom) is required")
  if (!is.numeric(q) || q < 1 || q > p) stop("`q` must be between 1 and p")
  
  # Standardize data matrix
  X <- scale(X, center = TRUE, scale = standardize)
  if (length(colnames(X)) == 0) {
    colnames(X) <- as.character(seq_len(p))
  }
  feature.names <- colnames(X)
  
  # Inform the user which allocation method is being used.
  if (q > 1 && !is.null(weight_method) && standardize) {
    cat(">> Using '", weight_method, "' allocation method.\n")
  }
  
  # --- 2. Sparsity Allocation ---
  # Calculate the total number of non-zero elements to be allocated.
  total_nonzeros <- p-dof
  
  n_vec <- calculate_allocation_vector(X=X,q=q,total_nonzeros = total_nonzeros,weight_method = weight_method)
  
  
  # --- 3. Compute Sparse Loadings ---
  # Use tryCatch to handle potential errors from the irlba package.
  svd_try <- tryCatch(
    suppressWarnings(irlba::ssvd(x = X, k = q, n = n_vec, maxit = max.iter)),
    error = function(e) e
  )
  if (inherits(svd_try, "error")) {
    stop("dof = ", dof, " does not fit the structure of the singular vectors. Use bdsvd.ht to find a suitable value for dof.")
  }
  V <- svd_try$v  # p × q matrix
  
  
  
  # --- 4. Extract Blocks ---
  # If the first sparse loading is fully dense, it means no separable structure was found.
  # Return the entire set of features as a single block.
  if (all(V[,1] != 0)) {
    result <- list(feature.names)
    class(result) <- "bdsvd"
    (result)
  }
  
  # Build binary mask for support of all q vectors
  mask <- multi.get.threshold.matrix(V, threshold = 0)
  
  # Extract blocks via greedy assignment
  blk_objs <- multi.get.blocks(mask, feature.names)
  
  # Gather feature names from each block object
  result <- lapply(blk_objs, slot, "features")
  class(result) <- "bdsvd"
  return(result)
}


#' Calculate Weighted Sparsity Allocation Vector
#'
#' @param X The data sub-matrix.
#' @param q The number of singular vectors.
#' @param total_nonzeros The total number of non-zero elements to allocate.
#' @param weight_method The allocation method ('equal' or 'explained_var').
#' @return A numeric vector of length `q` specifying the number of non-zeros for each loading.
#'
calculate_allocation_vector <- function(X, q, total_nonzeros, weight_method) {
  switch(
    weight_method,
    "equal" = {
      base_n <- floor(total_nonzeros / q)
      remainder <- total_nonzeros %% q
      vec <- rep(base_n, q)
      if (remainder > 0) vec[1:remainder] <- vec[1:remainder] + 1
      return(vec)
    },
    "explained_var" = {
      s <- svd(X, nu = 0, nv = q)
      explained_variance_weights <- s$d[1:q]^2
      weights <- if (all(explained_variance_weights == 0)) rep(1/q, q) else explained_variance_weights / sum(explained_variance_weights)
      n_proportional <- weights * total_nonzeros
      vec <- floor(n_proportional)
      remainder <- total_nonzeros - sum(vec)
      if (remainder > 0) {
        fractional_parts <- n_proportional - vec
        order_of_remainders <- order(fractional_parts, decreasing = TRUE)
        vec[order_of_remainders[1:remainder]] <- vec[order_of_remainders[1:remainder]] + 1
      }
      return(vec)
    },
    stop("Invalid `weight_method` provided.")
  )
}



