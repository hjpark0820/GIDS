
## EM_GMM_counts
## -------------------------------------------------------------------------
## Purpose
##   Estimate a two–component mixture on |R| (absolute correlations) using
##   an EM algorithm on histogram counts, after Fisher z-transform. The
##   components are truncated normals on [0, inf) with means (0, mu1) and sds
##   (sigma0, sigma1). From the fitted model, select two cutoffs:
##     epsilon1: threshold used for greedy peeling (structure),
##     epsilon2: higher threshold used to define “edges” in KL tests (strength).
##
## Inputs
##   counts   : numeric vector, histogram counts of |R| over given bin centers.
##   bins     : numeric vector in (0,1), bin centers for |R| (default seq(0.01,0.99,0.02)).
##   max_iter : integer, EM maximum iterations (default 100).
##   tol      : numeric, stop criterion on change of mixing weight (default 1e-6).
##   verbose  : logical, print EM diagnostics if TRUE.
##
## Output
##   list(epsilon1, epsilon2)  — both returned on the original correlation scale.
##
## Details
##   - Applies Fisher z: z = 0.5*log((1+|r|)/(1-|r|)), fits truncated normals on z≥0.
##   - E-step uses current mixture weights and truncated-normal pdfs at the bin grid.
##   - M-step updates lambda (mixing), mu1, sigma0, sigma1 with mild constraints:
##       sigma1 ≤ sigma00 (global sd bound) and sigma1 ≥ sigma0.
##   - epsilon1 maximizes a signal-vs-noise separation objective (obj) 
##   - epsilon2 maximizes an F1-like criterion (obj0) trading precision & sensitivity for the signal.
##   - Returns epsilon1, epsilon2 transformed back to the correlation domain via inverse Fisher z.
EM_GMM_counts <- function(counts, bins=seq(0.01,0.99,by=0.02), max_iter = 100, tol = 1e-6, verbose=FALSE) {
  
  
  library(truncnorm)
  
  # Fisher Z-transformation of grid
  bins = 0.5*log((1+bins)/(1-bins))
  # Maximum values of variance of noise and signal distributions
  sigma00 = sqrt(sum(counts * bins^2)/sum(counts))
  # Initial value of variance of signal distribution
  sigma1 = sqrt(sum(counts * bins^2)/sum(counts) - (sum(counts * bins)/sum(counts))^2)
  # Initial value of variance of noise distribution
  sigma0 = sigma1
  # Initial value of mean of the signal distribution
  mu1 <- 0.2       
  # Initial values of Mixing proportions
  lambda_mix <- c(0.99,0.01)
  # Temporary variable for the proportion of signal distribution
  lambda_mix_temp = 0
  
  # distribution of (transformed) absolute correlation variables (noise + signal) 
  dist = counts/sum(counts)
  # Initial distribution of noise
  dtrunc_x_mu0 <- dtruncnorm(bins, a = 0,  mean = 0, sd = sigma0)
  # Initial distribution of signal
  dtrunc_x_mu1 <- dtruncnorm(bins, a = 0,  mean = mu1, sd = sigma1)
  # grid for epsilon
  epsilon = seq(0,1,by=0.01)
  
  # estimated mean of conditional noise distribution: E[|R|||R|>cutoff,Z=0]
  mean0 = function(cutoff,mu0,sigma0){
    sum(epsilon[epsilon>cutoff]*dtruncnorm(epsilon[epsilon>cutoff],a=0, mean = 0, sd = sigma0) )*(epsilon[2])/(1-ptruncnorm(cutoff,a=0,b=1, mean = 0, sd = sigma0))
  }
  
  # estimated variance of conditional noise distribution: Var[|R|||R|>cutoff,Z=0]
  var0 = function(cutoff,mu0,sigma0){
    sum(epsilon[epsilon>cutoff]^2*dtruncnorm(epsilon[epsilon>cutoff],a=0, mean = 0, sd = sigma0))*(epsilon[2])/(1-ptruncnorm(cutoff,a=0,b=1, mean = 0, sd = sigma0)) - mean0(cutoff,0,sigma0)^2
  }
  
  # estimated mean of conditional signal distribution: E[|R|||R|>cutoff,Z=1]
  mean1 = function(cutoff,mu1,sigma1){
    sum(epsilon[epsilon>cutoff]*dtruncnorm(epsilon[epsilon>cutoff],a=0, mean = mu1, sd = sigma1))*(epsilon[2])/(1-ptruncnorm(cutoff,a=0,b=1, mean = mu1, sd = sigma1))
  }
  
  # estimated variance of conditional signal distribution: Var[|R|||R|>cutoff,Z=1]
  var1 = function(cutoff,mu1,sigma1){
    sum(epsilon[epsilon>cutoff]^2*dtruncnorm(epsilon[epsilon>cutoff],a=0,mean = mu1, sd = sigma1))*(epsilon[2])/(1-ptruncnorm(cutoff,a=0,b=1, mean = mu1, sd = sigma1)) - mean1(cutoff,mu1,sigma1)^2
  }
  
  # estimated CDF of conditional noise distribution: S0(cutoff)
  s0  = function(cutoff,mu0,sigma0){
    1-ptruncnorm(cutoff,a=0,b=1, mean = 0, sd = sigma0)
  }
  
  # estimated CDF of conditional signal distribution: S1(cutoff)
  s1  = function(cutoff,mu1,sigma1){
    1-ptruncnorm(cutoff,a=0,b=1, mean = mu1, sd = sigma1)
  }
  
  ## objective function to choose epsilon1
  obj=function(cutoff,pi1,mu0,sigma0,mu1,sigma1){
    num = (pi1-delta2)*(1-delta1-delta2)*(s1(cutoff,mu1,sigma1)*mean1(cutoff,mu1,sigma1) - s0(cutoff,0,sigma0)*mean0(cutoff,0,sigma0) )^2
    den1 = (sqrt((pi1-delta2)*(1-delta1-delta2)) + 2*delta2 )*var1(cutoff,mu1,sigma1)*s1(cutoff,mu1,sigma1) 
    den2 = (2*(1-delta2) +  sqrt((pi1-delta2)*(1-delta1-delta2)) )*var0(cutoff,0,sigma0)*s0(cutoff,0,sigma0)   
    num/(den1+den2)
  }
  
  ## objective function to choose epsilon2
  obj0=function(cutoff,pi1,mu0,sigma0,mu1,sigma1){
    pre = pi1*s1(cutoff,mu1,sigma1)/( (1-pi1)*s0(cutoff,0,sigma0)+pi1*s1(cutoff,mu1,sigma1))
    sen = s1(cutoff,mu1,sigma1)
    2*pre*sen/(pre+sen)
  }
  
  ## Start of EM algorithm ##
  for (iter in 1:max_iter) {
    
    # E-step: Compute responsibilities
    gamma1 <- lambda_mix[1] * dtrunc_x_mu0
    gamma2 <- lambda_mix[2] * dtrunc_x_mu1
    gamma_sum <- gamma1 + gamma2
    
    gamma1 <- gamma1 / gamma_sum
    gamma2 <- gamma2 / gamma_sum
    
    gamma <- cbind(gamma1, gamma2)
    gamma[is.na(gamma[,1]),1]=0
    gamma[is.na(gamma[,2]),2]=1
    
    # M-step: Update parameters
    Nk <- colSums(counts %*% gamma)
    lambda_mix <- Nk / sum(Nk)
    mu1 <- sum(gamma[, 2] * bins * counts) / Nk[2]
    
    
    # update the variance of noise distribution
    sigma0 = sqrt(sum(bins^2 * counts * gamma[,1])/sum(counts*gamma[,1]))
    # update the variance of signal distribution
    sigma1 = sqrt(sum((bins-mu1)^2 * counts * gamma[,2])/sum(counts *gamma[,2]))
    # if the variance of signal distribution is too big or too small, correct it. (estimation under some constraints)
    if( sigma1 > sigma00 ){
      sigma1 =  sigma00
    } else if (sigma0 > sigma1 ) {
      sigma1 = sigma0
    }
    
    # Update truncated normal densities (only for mu1 and sigma1)
    dtrunc_x_mu0 <- dtruncnorm(bins, a = 0,  mean = 0, sd = sigma0)
    dtrunc_x_mu1 <- dtruncnorm(bins, a = 0,  mean = mu1, sd = sigma1)
    
    # Compute log-likelihood and check convergence
    #new_loglik <- log_likelihood()
    
    if (abs(lambda_mix[1]-lambda_mix_temp) < tol) break
    #prev_loglik <- new_loglik
    lambda_mix_temp = lambda_mix[1]
    
    # update the proportion of signal distribution
    pi1 = lambda_mix[2]
    # Assumption about error rates
    delta1 = 0.1
    delta2 = pi1/2
    # calculate the values of objective function for epsilon1 and epsilon2 for each value in the grid of epsilon
    values1 = unlist(lapply(epsilon[-length(epsilon)],function(x) obj(x,pi1,0,sigma0,mu1,sigma1)))
    values2 = unlist(lapply(epsilon[-length(epsilon)],function(x) obj0(x,pi1,0,sigma0,mu1,sigma1)))
    if(verbose==TRUE) cat("iter: ",iter,"~",(exp(2*epsilon[which.max(values)])-1)/(exp(2*epsilon[which.max(values)])+1),lambda_mix[1],"sigma00: ",sigma00,"sigma0: ",sigma0,"mu1:~",mu1,"sigma1: ",sigma1,'\n')
  }
  
  # choose the cutoffs maximizing the objective functions 
  epsilon1 = epsilon[which.max(values1)]
  epsilon2 = epsilon[which.max(values2)]
  
  # inverse function of Z-transformation
  trans = function(x) {
    (exp(2*x)-1)/(exp(2*x)+1)
  }
  
  # return non-transformed epsilon values
  list(epsilon1=trans(epsilon1),epsilon2=trans(epsilon2))
}


## blockwise_cor_aggregate (used only for large data sets)
## -------------------------------------------------------------------------
## Purpose
##   Compute row/column sums of a hard-thresholded |cor(X, Y)| matrix without
##   materializing the full p×q matrix, by splitting X into blocks and
##   aggregating in parallel.
##
## Inputs
##   X, Y     : numeric matrices with the same number of rows (samples).
##   n_core   : integer, number of worker processes (default 16; one core is left free).
##   threshold: numeric in [0,1], hard threshold (use epsilon1 from EM_GMM_counts).
##
## Output
##   list(rowsum, colsum) where
##     rowsum[i] = sum_{j=1}^q 1{|cor(X_i, Y_j)| ≥ threshold}·|cor(X_i, Y_j)|
##     colsum[j] = sum_{i=1}^p 1{|cor(X_i, Y_j)| ≥ threshold}·|cor(X_i, Y_j)|
##
## Notes
##   - Splits columns of X into contiguous blocks; loops Y’s columns inside each task.
##   - Requires parallel backend (e.g., parallel::makeCluster, parLapply).
##   - Memory-friendly for large p×q; complexity dominated by repeated cor() calls.
##   - Ensure n_core ≤ detectCores()-1 to keep the machine responsive.
blockwise_cor_aggregate <- function(X, Y, n_core = 16, threshold = 0.1) {
  
  # set p and q
  p <- ncol(X)
  q <- ncol(Y)
  
  # set the size of each block
  block_size_x = ceiling(p/n_core)
  
  # Initialize accumulators
  row_sums <- numeric(p)
  col_sums <- numeric(q)
  
  # Create index blocks
  x_blocks <- split(1:p, ceiling(seq_along(1:p) / block_size_x))
  
  # setup the working threads
  if (n_core > detectCores() - 1) num_cores = detectCores() - 1 
  cl <- makeCluster(n_core)
  
  # Export necessary variables to workers
  clusterExport(cl, varlist = c("X", "Y", "x_blocks", "threshold"), envir = environment())
  
  # Compute all block pairs
  results <- parLapply(cl, seq_along(x_blocks), function(i_block) {
    
    ## Assign the row and column sum vectors with the length p and q, respectively
    row_sum <- numeric(ncol(X))
    col_sum <- numeric(ncol(Y))
    
    # set the x index of assigned block
    x_idx <- x_blocks[[i_block]]
    
    # Calculate row and column sums of the assigned block
    for(j in 1:ncol(Y)){
      vec = abs(cor(X[, x_idx, drop = FALSE], Y[, j, drop = FALSE]))
      vec[vec<threshold] = 0
      row_sum[x_idx] <- row_sum[x_idx] + vec
      col_sum[j] <- sum(vec)
    }
    
    list(row_sum = row_sum, col_sum = col_sum)
  })
  
  stopCluster(cl)
  
  # Aggregate across all blocks
  for (res in results) {
    row_sums <- row_sums + res$row_sum
    col_sums <- col_sums + res$col_sum
  }
  
  # return row and column sums
  list(rowsum = row_sums,colsum = col_sums)
}


## Purpose
##   Score a collection of extracted subgraphs by a normalized KL divergence
##   between edge rates inside vs. outside the union of subgraphs, using
##   a binary adjacency W built at epsilon2. Returns KLdiv / (entropy normalizer).
##
## Inputs
##   W     : 0/1 adjacency matrix derived from |cor| at epsilon2 (same dimensions as W1).
##   index : list of subgraphs; each element has $x (rows), $y (cols).
##
## Output
##   numeric scalar: normalized KL (nKL). Returns -Inf if no valid subgraph.
##
## Details
##   - p  = overall edge rate; p1 = edge rate inside union of subgraphs;
##     p0 = edge rate outside.
##   - Normalizer = sum of binary entropies weighted by sizes inside/outside,
##     to keep nKL dimensionless and comparable across selections.
##   - Edge cases handled: empty subgraphs, p1 = 1, or degenerate coverage.
kld_cross00 = function(W,index){
  
  ## Set the number of extracted subgraphs
  n = length(index)
  
  # Calculate the normalized KL-divergence (nKL-div)
  if(length(index[[1]]$x)==0||length(index[[1]]$y)==0) {
    KLdiv = -Inf # if no subgraph is extracted, nKL-div = -infinite
    normalizer = 1 # normalizer: entropy
  } else {
    
    p = mean(W) # p: overall mean of edge connections in the entire graph
    sum1= c() # vector for the numbers of edges in subgraphs
    size = c() # vector for the sizes of subgraphs
    for(i in 1:n){
      sum1[i] = sum(W[index[[i]]$x,index[[i]]$y]) 
      size[i] = length(index[[i]]$x)*length(index[[i]]$y)
    }
    
    p1 = sum(sum1)/sum(size) # p1: the average of edge connection in the entire extracted subgraphs
    W = as.matrix(W)
    p0 = (sum(W) - sum(sum1) ) / (dim(W)[1]*dim(W)[2] - sum(size) ) # p0: the average of edge connection outside the extracted subgraphs
    
    # normalizer (entropy) calculation
    if (p1 < 1){
      normalizer = - sum(size)*(p1*log(p1) + (1-p1)*log(1-p1)) - (length(W)-sum(size))*(p0*log(p0) + (1-p0)*log(1-p0)) 
    } else {
      normalizer = - (length(W)-sum(size))*(p0*log(p0) + (1-p0)*log(1-p0))
    }
    
    # KL-divergence calculation
    KLdiv = 0
    
    if (length(index[[1]]$x)==0 || length(index[[1]]$y)==0) {
      KLdiv = - Inf
    } else  if (p1 == 1){
      KLdiv = sum(size)*(p1*log(p1/p) ) + (length(W)-sum(size))*(p0*log(p0/p) +  (1-p0)*log((1-p0)/(1-p))) 
    } else if (length(index$x)+length(index$y) < dim(W)[1]+dim(W)[2] ){
      KLdiv =  sum(size)*(p1*log(p1/p) +  (1-p1)*log((1-p1)/(1-p))) + (length(W)-sum(size))*(p0*log(p0/p) +  (1-p0)*log((1-p0)/(1-p))) 
    } else {
      KLdiv = -Inf
    }
    
  }
  
  # Return normalized KL-divergence
  KLdiv/normalizer
}


## subg22_s1  — GIDS Phase 1 (screening)
## -------------------------------------------------------------------------
## Purpose
##   Dimensionality reduction (pre-screening) for X and Y by iteratively
##   removing k1 weakest rows or k2 weakest columns (by thresholded |cor|)
##   until target sizes (ps1, qs1) are reached.
##
## Inputs
##   X, Y  : numeric data matrices (same number of rows/samples).
##   rsum  : length-p vector, row sums of hard-thresholded |cor(X, Y)|.
##   csum  : length-q vector, column sums of hard-thresholded |cor(X, Y)|.
##   cut   : numeric, threshold epsilon1 used to compute rsum/csum and updates.
##   ps1   : target number of X features to keep.
##   qs1   : target number of Y features to keep.
##   k1    : rows removed per iteration (granularity on X side).
##   k2    : cols removed per iteration (granularity on Y side).
##   verbose: logical, print progress if TRUE.
##
## Output
##   list(x, y) — screened indices of X and Y to pass to Phase 2.
##
## Details
##   - At each step, compute average weight of k1 lightest rows vs. k2 lightest
##     columns (w.r.t. current thresholded correlations). Remove the side with
##     the smaller average, and update the opposite side’s sums by subtracting
##     contributions of the removed block (recomputing only the needed pieces).
##   - Stops when remaining counts hit targets (≤ps1 and/or ≤qs1).
subg22_s1 <- function(X, Y, rsum, csum, cut, ps1, qs1, k1, k2,verbose=FALSE) {
  
  nr = length(rsum) ## number of rows
  nc = length(csum) ## number of columns
  row <- nr ## initial number of remaining rows
  col <- nc ## initial number of remaining columns
  R <- rsum ## initial row sums for active rows and columns
  C <- csum ## initial col sums for active rows and columns
  r_which_0 = integer(0) ## initial vector for excluded rows
  c_which_0 = integer(0) ## initial vector for excluded columns
  
  # Start of row and column exclusion
  # stop criterion: remaining rows and columns are equal or smallter than target dimensions
  while(length(r_which_0) < nr-ps1 || length(c_which_0) < nc - qs1 ){
    
    row <- nr - length(r_which_0) ## update the number of active rows
    col <- nc - length(c_which_0) ## update the number of active columns.
    
    ## IndR: indice of rows with the k1 minimum sums
    if (length(r_which_0)>0){
      IndR <- ((1:nr)[-r_which_0])[order(R[-r_which_0])[1:k1]]
    } else {
      IndR <- order(R)[1:k1]
    }
    
    ## calculating the average row sums of chosen rows
    resultR <- R[IndR]
    dR <- sum(resultR)/(col*k1) ## average row sum of chosen rows
    
    ## Choosing the columns with the k2 minimum sums
    ## IndC: indice of columns with the k2 minimum sums
    if (length(c_which_0)>0){
      IndC <- ((1:nc)[-c_which_0])[order(C[-c_which_0])[1:k2]]
    } else {
      IndC <- order(C)[1:k2]
    }
    
    ## calculating the average column sums of chosen columns
    resultC <- C[IndC]
    dC <- sum(resultC)/(row*k2) ## average column sum of chosen columns
    
    ## update on the active rows and columns
    if (row <= ps1) { ## case that the row dimension is already sufficiently reduced
      if(length(c_which_0)>0){
        c_which_0 = c(c_which_0,((1:nc)[-c_which_0])[order(C[-c_which_0])[1:(col-qs1)]])
      } else {
        c_which_0 = order(C)[1:qs1]
      }
    } else if (col <= qs1) {   ## case that the column dimension is already sufficiently reduced 
      if(length(r_which_0)>0) {
        r_which_0 = c(r_which_0,((1:nr)[-r_which_0])[order(R[-r_which_0])[1:(row-ps1)]])  
      } else {
        r_which_0 = order(R)[1:ps1]
      }
      
    } else if ( dR <= dC && row > ps1) { ## case that rows are excluded
      if(length(c_which_0)>0){
        cor_temp = abs(cor(X[,IndR],Y[,-c_which_0]))
        cor_temp[cor_temp<cut]=0
        C[-c_which_0] <- C[-c_which_0] - colSums(cor_temp)
      } else { 
        cor_temp = abs(cor(X[,IndR],Y))
        cor_temp[cor_temp<cut]=0
        C <- C - colSums(cor_temp)
      }
      r_which_0 = c(r_which_0,IndR)
      row <- nr - length(r_which_0)
      R[IndR] <- 0
    } else { ## case that columns are excluded
      if(length(r_which_0)>0){
        cor_temp = abs(cor(X[,-r_which_0],Y[,IndC]))
        cor_temp[cor_temp<cut]=0
        R[-r_which_0] <- R[-r_which_0] - rowSums(cor_temp)
      } else {
        cor_temp = abs(cor(X,Y[,IndC]))
        cor_temp[cor_temp<cut]=0
        R <- R - rowSums(cor_temp)
      }
      c_which_0 = c(c_which_0, IndC)
      col <- nc - length(c_which_0)
      C[IndC] <- 0
    }
    
    if(verbose==TRUE) cat("stage1:", "r removed: ",length(r_which_0), "c removed: ",length(c_which_0)  ,"\n")
  }  
  
  result=list()
  result$x=(1:nr)[-sort(r_which_0)]
  result$y=(1:nc)[-sort(c_which_0)]
  
  return(result)
  
}

## GIDS Phase 2: Extract subgraphs from a bipartite weight matrix by greedy peeling
## ------------------------------------------------------------------------------
## Inputs
##   W1     : numeric matrix (|X| x |Y|). Absolute correlation (or general weights).
##   lambda : numeric vector of lambda exponents for density scoring.
##   cut1   : numeric. First hard-threshold (epsilon1); entries < cut1 are zeroed before peeling.
##   cut2   : numeric. Second threshold (epsilon2); forwarded to subg22_2 for refinement.
##   index0 : list with entries $x, $y giving candidate row/col indices to consider.
##            If empty, use all rows/cols of W1.
##   gamma0 : numeric (passed to subg22_2). Typically controls additional criteria (e.g., size/penalty).
##   verbose: logical. If TRUE, prints progress.
##
## Output
##   A list of length length(lambda). For each lambda = lambda[i], you get a list of one
##   or more subgraphs (components). Each subgraph j has fields:
##       result[[i]][[j]]$x : row indices kept for subgraph j (in W1’s coordinate system)
##       result[[i]][[j]]$y : col indices kept for subgraph j (in W1’s coordinate system)
##
## Notes
##   - The greedy “peeling” compares row and column average weights: dS = R_s / (#cols),
##     dT = C_t / (#rows). It removes the side with the smaller average (sparser node).
##   - For each step, it records (#rows, #cols, which-side-removed, removed row/col id)
##     and the lambda-density score. The best step (per lambda) determines the first subgraph.
##   - Then it repeatedly calls subg22_2 (a refinement/extraction step) to pull out
##     additional subgraphs from the remaining index set until convergence.
##
subg22_1 <- function(W1, lambda, cut1, cut2, index0 = list(), gamma0, verbose = FALSE) {
  
  # If no candidate indices provided, default to full matrix
  if (length(index0) == 0) {
    index0 <- list(x = 1:nrow(W1), y = 1:ncol(W1))
  }
  
  # 1) Hard-threshold W1 at cut1 (epsilon1) to form the working matrix W
  #    W00 is a temp copy; entries < cut1 are zeroed (hard thresholding)
  W00 <- W1
  W00[W00 < cut1] <- 0
  W <- W00
  
  # Bookkeeping containers
  deleteS_list <- list() # (not used directly here; kept for API symmetry/extension)
  deleteT_list <- list() # (not used directly here; kept for API symmetry/extension)
  
  nr <- nrow(W) # total rows
  nc <- ncol(W) # total cols
  
  # density_list stores, per peeling iteration i:
  #   [,1]  #rows after step i
  #   [,2]  #cols after step i
  #   [,3]  which side peeled at step i (1 = row, 2 = col)
  #   [,4:5] (unused placeholders)
  #   [,6]  removed row index at step i (if any; else 0)
  #   [,7]  removed col index at step i (if any; else 0)
  #   [,8:(7+length(lambda))]  lambda-density score at step i for each lambda
  density_list <- matrix(0, nrow = (nr + nc - 1), ncol = 7 + length(lambda))
  
  # Current counts (start with full)
  row <- nr
  col <- nc
  
  # Row/col sum vectors over the *current* working matrix
  R <- rowSums(W)  # row weights
  C <- colSums(W)  # col weights
  
  # Sets of indices that have been peeled (removed) so far
  r_which_0 <- integer(0)
  c_which_0 <- integer(0)
  
  # Record initial state (i = 1)
  density_list[1, 1] <- row
  density_list[1, 2] <- col
  # Initial lambda-density: total weight divided by (rows*cols)^lambda
  density_list[1, 8:(7 + length(lambda))] <- sum(R) / ( (density_list[1,1] * density_list[1,2])^lambda )
  
  # 2) Greedy peeling loop: remove one row or one column at a time
  #    We run up to nr + nc - 2 peel steps (at least one row and one col must remain)
  for (i in 2:(nr + nc - 2)) {
    
    # Early stop if almost everything is removed on either side
    if (length(r_which_0) > nr - 2 || length(c_which_0) > nc - 2) break
    
    # Current active sizes
    row <- nr - length(r_which_0)
    col <- nc - length(c_which_0)
    
    # Pick the *active* row with minimal row-sum R
    if (length(r_which_0) > 0) {
      active_rows <- (1:nr)[-r_which_0]
      IndS <- active_rows[ which.min(R[active_rows]) ]
    } else {
      IndS <- which.min(R)
    }
    resultR <- R[IndS]
    dS <- resultR / col  # row’s average weight across current columns
    
    # Pick the *active* column with minimal col-sum C
    if (length(c_which_0) > 0) {
      active_cols <- (1:nc)[-c_which_0]
      IndT <- active_cols[ which.min(C[active_cols]) ]
    } else {
      IndT <- which.min(C)
    }
    resultC <- C[IndT]
    dT <- resultC / row  # col’s average weight across current rows
    
    # Peel the side with the *smaller* average ⇒ sparser node
    if (dS <= dT) {
      # Remove row IndS:
      #   - subtract its contribution from column sums
      #   - zero out the row in W
      #   - mark the row as removed in R
      if (length(IndS) > 1) {
        # (defensive branch; IndS should be scalar)
        C <- C - colSums(W[IndS, ])
      } else {
        C <- C - W[IndS, ]
      }
      r_which_0 <- c(r_which_0, IndS)
      row <- nr - length(r_which_0)
      W[IndS, ] <- 0
      R[IndS] <- 0
      
      # Log this step
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 1       # peeled a row
      density_list[i, 6] <- IndS    # which row was removed
      
    } else {
      # Remove column IndT (symmetrically to the row case)
      if (length(IndT) > 1) {
        R <- R - rowSums(W[, IndT])
      } else {
        R <- R - W[, IndT]
      }
      c_which_0 <- c(c_which_0, IndT)
      col <- nc - length(c_which_0)
      W[, IndT] <- 0
      C[IndT] <- 0
      
      # Log this step
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 2       # peeled a column
      density_list[i, 7] <- IndT    # which column was removed
    }
    
    # Update lambda-density scores after this peel
    density_list[i, 8:(7 + length(lambda))] <- sum(R) / ( (density_list[i,1] * density_list[i,2])^lambda )
  }
  
  if (verbose) cat("k=1 done", "\n")
  
  # Dimensions of (thresholded) W at end of peeling
  xdim <- nrow(W)
  ydim <- ncol(W)
  
  lambda_dim <- length(lambda)
  result <- list()
  
  # 3) For each lambda, locate the step with maximum lambda-density ⇒ define first subgraph
  for (i in 1:lambda_dim) {
    max_ind <- which.max(density_list[, 7 + i])  # best step for this lambda
    
    # Reconstruct sets removed up to max_ind
    x_removed <- unique(density_list[1:max_ind, 6])
    y_removed <- unique(density_list[1:max_ind, 7])
    x_removed <- x_removed[x_removed > 0]
    y_removed <- y_removed[y_removed > 0]
    
    # Indices kept for the first subgraph are the complement of removed ones
    keep_x <- (1:xdim)[-x_removed]
    keep_y <- (1:ydim)[-y_removed]
    
    result[[i]] <- list()
    result[[i]][[1]] <- list(x = keep_x, y = keep_y)
    if (verbose) cat("i:", i, "~", length(keep_x), length(keep_y), "\n")
  }
  
  # 4) Prepare disjoint extraction: convert “kept in W” → “removed in W”
  #    Then map those to original indices via complements.
  #    Here, 'index' holds the *remaining* pool to search for further subgraphs.
  index <- lapply(
    1:lambda_dim,
    function(z) list(
      # Note: result[[z]][[1]]$x are 'kept' indices in 1:xdim;
      #       (1:nr)[- ... ] maps to original row index set that were removed in W,
      #       which becomes the pool to mine next (to avoid overlap)
      x = (1:nr)[-result[[z]][[1]]$x],
      y = (1:nc)[-result[[z]][[1]]$y]
    )
  )
  
  if (verbose) cat("first graph is extracted", "\n")
  
  # 5) For each lambda, iteratively extract additional subgraphs via subg22_2
  is_done <- rep(FALSE, length(index))
  for (i in 1:length(index)) {
    
    # If nothing left on either side, stop for this lambda
    if (length(index[[i]]$x) == 0 || length(index[[i]]$y) == 0) {
      is_done[i] <- TRUE
    }
    
    k <- 2
    while (!is_done[i]) {
      # subg22_2 should return a (possibly refined) subgraph from current pool index[[i]]
      result11 <- subg22_2(W1, index[[i]], lambda[i], cut1, cut2, gamma0)
      
      # If no change (same size), we are done for this lambda
      if (length(result11$x) == length(index[[i]]$x)) {
        is_done[i] <- TRUE
      } else {
        # Record the new subgraph
        result[[i]][[k]] <- list(x = result11$x, y = result11$y)
        
        # Remove the extracted nodes from the pool so subsequent subgraphs are disjoint
        index[[i]]$x <- index[[i]]$x[ !(index[[i]]$x %in% result11$x) ]
        index[[i]]$y <- index[[i]]$y[ !(index[[i]]$y %in% result11$y) ]
        
        if (verbose) cat(i, "th value of lambda,", "for", k, "th graph is extracted", "\n")
        k <- k + 1
      }
    }
  }
  
  return(result)
}

## GIDS Phase 2 (subfunction): Refine/extract one subgraph via greedy peeling + KL test
## ------------------------------------------------------------------------------
## Inputs
##   W       : numeric matrix (|X| x |Y|). Full weight matrix in the *original* coordinates.
##   index0  : list with $x, $y (integer vectors). Candidate row/col indices (in W’s coordinates)
##             that define the current search window W[index0$x, index0$y].
##   lambda0 : numeric scalar lambda used in the density score: sum(R) / (rows*cols)^lambda.
##   cut1    : numeric. First threshold (epsilon1) for hard-thresholding during peeling.
##   cut2    : numeric. Second (higher) threshold (epsilon2) to define edges in the KL test.
##   gamma0  : numeric in (0,1). Baseline/nominal Bernoulli edge probability for the KL test.
##
## Output
##   result : list with fields
##              $x : selected row indices (in W’s coordinate system)
##              $y : selected col indices (in W’s coordinate system)
##            Either the refined subgraph or a fallback to the input index0.
##
subg22_2 <- function(W, index0, lambda0, cut1, cut2, gamma0) {
  
  # 1) Work on the candidate submatrix; hard-threshold at cut1 for peeling.
  W0 <- abs(W[index0$x, index0$y])
  W0 <- as.matrix(W0)
  W0[W0 < cut1] <- 0
  
  # Bookkeeping
  deleteS_list <- list()  # (reserved; not used here)
  deleteT_list <- list()  # (reserved; not used here)
  
  nr <- nrow(W0)
  nc <- ncol(W0)
  
  # density_list layout (nrow ≈ nr+nc-1, ncol = 9):
  #   [,1]  rows after step i
  #   [,2]  cols after step i
  #   [,3]  which side peeled (1 = row, 2 = col)
  #   [,4:5] unused placeholders
  #   [,6]  removed row id at step i (else 0)
  #   [,7]  removed col id at step i (else 0)
  #   [,8]  lambda-density at step i
  density_list <- matrix(0, nrow = (nr + nc - 1), ncol = 8)
  
  # Current sizes and sums
  row <- nr
  col <- nc
  R <- rowSums(W0)
  C <- colSums(W0)
  
  # Sets of removed rows/cols (indices in 1:nr or 1:nc)
  r_which_0 <- integer(0)
  c_which_0 <- integer(0)
  
  # Initial state (i = 1)
  density_list[1, 1] <- row
  density_list[1, 2] <- col
  density_list[1, 8] <- (sum(R)) / ( (density_list[1,1] * density_list[1,2])^lambda0 )
  
  # 2) Greedy peeling loop: remove one row or one column at a time by smaller average weight
  for (i in 2:(nr + nc - 2)) {
    
    # Stop if only one row or column is left on either side
    if (length(r_which_0) > nr - 2 || length(c_which_0) > nc - 2) break
    
    row <- nr - length(r_which_0)
    col <- nc - length(c_which_0)
    
    # Select active row with minimal row-sum
    if (length(r_which_0) > 0) {
      active_rows <- (1:nr)[-r_which_0]
      IndS <- active_rows[ which.min(R[active_rows]) ]
    } else {
      IndS <- which.min(R)
    }
    resultR <- R[IndS]
    dS <- resultR / col  # row’s average weight
    
    # Select active column with minimal col-sum
    if (length(c_which_0) > 0) {
      active_cols <- (1:nc)[-c_which_0]
      IndT <- active_cols[ which.min(C[active_cols]) ]
    } else {
      IndT <- which.min(C)
    }
    resultC <- C[IndT]
    dT <- resultC / row  # col’s average weight
    
    # Peel the sparser side
    if (dS <= dT) {
      # Remove row IndS
      if (length(IndS) > 1) {
        C <- C - colSums(W0[IndS, ])
      } else {
        C <- C - W0[IndS, ]
      }
      r_which_0 <- c(r_which_0, IndS)
      row <- nr - length(r_which_0)
      W0[IndS, ] <- 0
      R[IndS] <- 0
      
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 1       # peeled a row
      density_list[i, 6] <- IndS    # removed row id
      
    } else {
      # Remove column IndT
      if (length(IndT) > 1) {
        R <- R - rowSums(W0[, IndT])
      } else {
        R <- R - W0[, IndT]
      }
      c_which_0 <- c(c_which_0, IndT)
      col <- nc - length(c_which_0)
      W0[, IndT] <- 0
      C[IndT] <- 0
      
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 2       # peeled a column
      density_list[i, 7] <- IndT    # removed col id
    }
    
    # Update lambda-density after this peel
    density_list[i, 8] <- (sum(R)) / ( (density_list[i,1] * density_list[i,2])^lambda0 )
  }
  
  # 3) Identify the best peel step by lambda-density and reconstruct kept indices (local coords)
  dST  <- max(density_list[, 8])
  indST <- which.max(density_list[, 8])
  
  x_removed <- unique(density_list[1:indST, 6])
  y_removed <- unique(density_list[1:indST, 7])
  x_removed <- x_removed[x_removed > 0]
  y_removed <- y_removed[y_removed > 0]
  
  # Kept indices in the *local* (1:nr, 1:nc) frame
  keep_x_local <- (1:length(index0$x))[-x_removed]
  keep_y_local <- (1:length(index0$y))[-y_removed]
  
  result_temp <- list(x = keep_x_local, y = keep_y_local)
  result <- list()
  
  # 4) Guards: tiny subgraphs or no effective reduction → fall back to index0
  if (length(result_temp$x) < 3 || length(result_temp$y) < 3) {
    result$x <- index0$x; result$y <- index0$y
  } else if (length(result_temp$x) == length(index0$x) ||
             length(result_temp$y) == length(index0$y)) {
    result$x <- index0$x; result$y <- index0$y
  } else {
    # 5) KL-style significance test at the *higher* threshold cut2
    #    Define Bernoulli “edge” as (|W| > cut2) within the proposed submatrix.
    #    gamma  = observed edge density in the candidate;
    #    gamma0 = baseline edge probability (null).
    #    The tail bound 'delta' penalizes also by sizes relative to the global dims.
    gamma <- mean(abs(W[index0$x, index0$y])[result_temp$x, result_temp$y] > cut2)
    if (gamma == 1) {
      gamma <- 1 - 0.1^10  # avoid log(0) in KL
    }
    
    # KL divergence for Bernoulli: D(gamma || gamma0) =
    #   gamma * log(gamma/gamma0) + (1-gamma) * log((1-gamma)/(1-gamma0))
    # Composite exponent also includes structural penalties that depend on
    # global dims (dim(W)[1], dim(W)[2]) vs. selected sizes.
    delta <- exp(
      - length(result_temp$x) * length(result_temp$y) * (
        gamma * log(gamma / gamma0) +
          (1 - gamma) * log((1 - gamma) / (1 - gamma0)) -
          log(exp(1) * dim(W)[1] / length(result_temp$x)) / length(result_temp$y) -
          log(exp(1) * dim(W)[2] / length(result_temp$y)) / length(result_temp$x)
      )
    )
    
    
    # Decision rule:
    #  - If delta is NA (numerical issue) or not small enough (> 0.001), reject the refinement (fallback).
    #  - Else accept: map local kept indices to original coordinates.
    if (is.na(delta)) {
      result$x <- index0$x; result$y <- index0$y
    } else if (delta > 0.001) {
      result$x <- index0$x; result$y <- index0$y
    } else {
      result$x <- index0$x[-x_removed]
      result$y <- index0$y[-y_removed]
    }
  }
  
  return(result)
}
