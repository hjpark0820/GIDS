library(Matrix)
library(RSpectra)
library(sRDA)
library(PMA)
library(Ball)
library(truncnorm)


# from the same folder
setwd("folder path") ### your folder path that function2.R is saved
source("function2.R")

## ---------------------------
## Synthetic bipartite setup
## ---------------------------

nx <- 4000   # |X| = total number of X variables (columns in X)
ny <- 1500   # |Y| = total number of Y variables (columns in Y)

# Two planted subgraphs with within- and cross- block correlations cor0[i]
cor0 <- c(0.3, 0.3)
add_cor = -0.3
x_num_sub <- c(20,30)  # sizes of subgraphs along X
y_num_sub <- c(30,40)  # sizes of subgraphs along Y

nx_sub <- sum(x_num_sub)          # total 'signal' X features
ny_sub <- sum(y_num_sub)          # total 'signal' Y features
x_uncor <- nx - nx_sub            # 'noise' X
y_uncor <- ny - ny_sub            # 'noise' Y

# Global correlation structures (block-diagonal within X, within Y, and block-wise cross X-Y)
n <- nx + ny
cormat <- matrix(0, n, n)         # full correlation matrix placeholder (not directly used later)
Sigma_x  <- matrix(0, nx, nx)     # within-X correlation (to fill with block structure)
Sigma_y  <- matrix(0, ny, ny)     # within-Y correlation
Sigma_xy <- matrix(0, nx, ny)     # cross-correlation X-Y


## Fill X blocks with cor0[i]
initial <- 1
for (i in 1:length(x_num_sub)) {
  Sigma_x[initial:cumsum(x_num_sub)[i], initial:cumsum(x_num_sub)[i]] <- cor0[i]
  initial <- cumsum(x_num_sub)[i] + 1
}

## Fill Y blocks with cor0[i]
initial <- 1
for (i in 1:length(y_num_sub)) {
  Sigma_y[initial:cumsum(y_num_sub)[i], initial:cumsum(y_num_sub)[i]] <- cor0[i]
  initial <- cumsum(y_num_sub)[i] + 1
}

## Fill cross X-Y blocks with cor0[i]
initialx <- 1
initialy <- 1
for (i in 1:length(y_num_sub)) {
  Sigma_xy[initialx:cumsum(x_num_sub)[i], initialy:cumsum(y_num_sub)[i]] <- cor0[i]
  initialx <- cumsum(x_num_sub)[i] + 1
  initialy <- cumsum(y_num_sub)[i] + 1
}

## Randomly “punch holes” (zero-out) inside each planted block to increase heterogeneity

xi = c(0,cumsum(x_num_sub))
yi = c(0,cumsum(y_num_sub))

for(i in 1:length(x_num_sub)){
  xindex = (xi[i]+1):xi[i+1]
  yindex = (yi[i]+1):yi[i+1]
  x1 = sample(xindex,size=floor(x_num_sub[i]/5))
  y1 = sample(yindex,size=floor(y_num_sub[i]/5))
  Sigma_x[x1,x1] = Sigma_x[x1,x1]-add_cor
  Sigma_y[y1,y1] = Sigma_y[y1,y1]-add_cor
  Sigma_xy[x1,y1] = Sigma_xy[x1,y1] +add_cor
  
}

diag(Sigma_x) <- 1
diag(Sigma_y) <- 1

## Compose full (signal-only) covariance for first nx_sub + ny_sub variables,
## and also a larger Sigma1 that includes noise positions (not used directly below)
Sigma_x00 <- Sigma_x[1:nx_sub, 1:nx_sub]
Sigma_y00 <- Sigma_y[1:ny_sub, 1:ny_sub]
Sigma_xy00 <- Sigma_xy[1:nx_sub, 1:ny_sub]
Sigma1 <- rbind(cbind(Sigma_x, Sigma_xy), cbind(t(Sigma_xy), Sigma_y))
Sigma  <- rbind(cbind(Sigma_x00, Sigma_xy00), cbind(t(Sigma_xy00), Sigma_y00))

#Sigma_y00[1:15,1:15]
#Sigma_y00[16:35,16:35]

## Storage across replicates
results  <- list()  # stores final Phase-2 picks & metadata
results2 <- list()  # Phase-1 screened set with size 2d
results4 <- list()  # Phase-1 screened set with size 4d

## ---------------------------
## simulation loop
## ---------------------------

#library(matrixcalc)
#is.positive.semi.definite(Sigma_y00)

num_episods = 10

for (k in 1:num_episods) {
  
  ## sample size of X and Y
  sample_n <- 200
  
  ## Draw Gaussian samples for the signal block only; append Gaussian noise features
  ## (Signal block correlates as Sigma; noise features are independent N(0,1))
  results1 <- mvtnorm::rmvnorm(n = sample_n, sigma = Sigma)
  X = cbind(results1[,1:nx_sub],matrix(rnorm(sample_n*x_uncor),sample_n, x_uncor))
  Y = cbind(results1[,(nx_sub+1):(nx_sub+ny_sub)],matrix(rnorm(sample_n*y_uncor),sample_n, y_uncor))
  
  ## Empirical absolute cross-correlation
  W00 <- abs(cor(X, Y))
  
  ## (epsilon1, epsilon2) selection via EM_GMM_counts using histogram of |cor|
  breaks  <- seq(0, 1, by = 0.02)
  bins    <- (breaks[-1] + breaks[-length(breaks)]) / 2
  hist2   <- table(cut(W00, breaks = breaks, right = FALSE))
  
  cut00 <- EM_GMM_counts(as.numeric(hist2), bins, max_iter = 10000)
  ## IMPORTANT: EM_GMM_counts returns list(epsilon1=..., epsilon2=...)
  cut1 <- cut00$epsilon1   # structure threshold epsilon1  (used for peeling & screening)
  cut2 <- cut00$epsilon2   # edge threshold epsilon2      (used for binarization/KL tests)
  
  ## Threshold W00 at epsilon1 to compute row/column sums for Phase-1
  W00[W00 < cut1] <- 0
  rsum <- rowSums(W00)
  csum <- colSums(W00)
  
  ## Screening size parameters: d = n/log n ; Phase-1 will keep about 2d or 4d per side
  d <- ceiling(sample_n / log(sample_n))
  
  ## -------- Phase 1 (subg22_s1): greedy screening using epsilon1 and (rsum, csum) --------
  ## screening size  = 4d
  results4[[k]] <- gids_p1(X, Y, rsum, csum, cut = cut1, ps1 = 4 * d, qs1 = 4 * d, k1 = 13, k2 = 17, verbose=TRUE) 
  ## screening size  = 2d
  results2[[k]] <- gids_p1(X, Y, rsum, csum, cut = cut1, ps1 = 2 * d, qs1 = 2 * d, k1 = 10, k2 = 10)
  
  ## -------- Phase 2 (subg22_1): greedy peeling (lambda-density) + subg22_2 refinement --------
  lambda <- seq(0.5, 0.9, by = 0.1)          # grid of lambda exponents for density score
  ## Work within the Phase-1=2d window; pass γ0 = global edge rate above epsilon2 (within that window)
  W_phase2 <- abs(cor(X[, results2[[k]]$x], Y[, results2[[k]]$y]))
  gamma0   <- sum(W_phase2 > cut2) / length(W_phase2)
  
  result2 <- gids_p2(
    W1     = W_phase2,
    lambda = lambda,
    cut1   = cut1,
    cut2   = cut2,
    index0 = results2[[k]],      # (global indices) used internally for mapping
    gamma0 = gamma0,
    verbose = FALSE
  )
  
  ## -------- Model selection over epsilon⋆ using normalized KL (kld_cross00) --------
  ## We sweep a binarization threshold epsilon⋆ (≥ epsilon2 typically), pick (lambda, epsilon⋆) maximizing nKL.
  threshold0 <- seq(0.05, 0.5, by = 0.05)
  kld_values <- matrix(0, length(lambda), length(threshold0))
  
  for (i in 2:length(threshold0)) {
    W0 <- matrix(0, nrow(W_phase2), ncol(W_phase2))
    W0[W_phase2 > threshold0[i]] <- 1
    # For each lambda’s extracted list of subgraphs, compute normalized KL
    kld_values[, i] <- vapply(result2, function(z) kld_cross00(W0, z), numeric(1))
  }
  
  ## Choose best (lambda, epsilon⋆)
  max_index   <- arrayInd(which.max(kld_values), dim(kld_values))
  best_lambda <- lambda[max_index[1]]
  epsilon_star <- threshold0[max_index[2]]
  best_subgraphs <- result2[[max_index[1]]]
  
  ## Collect all (possibly multiple) disjoint subgraphs selected for the best lambda
  results[[k]] <- list()
  results[[k]]$var <- best_subgraphs
  results[[k]]$x <- best_subgraphs[[1]]$x
  results[[k]]$y <- best_subgraphs[[1]]$y
  
  if (length(best_subgraphs) > 1) {
    for (j in 2:length(best_subgraphs)) {
      results[[k]]$x <- c(results[[k]]$x, best_subgraphs[[j]]$x)
      results[[k]]$y <- c(results[[k]]$y, best_subgraphs[[j]]$y)
    }
  }
  
  results[[k]]$x <- sort(results[[k]]$x)
  results[[k]]$y <- sort(results[[k]]$y)
  results[[k]]$v_num <- c(length(results[[k]]$x), length(results[[k]]$y))
  results[[k]]$lambda <- best_lambda
  
  cat("k: ", k, "  best lambda=", best_lambda, "  epsilon*=", epsilon_star, "\n")
}

## ---------------------------
## Metrics: Sensitivity / Precision / F1
## ---------------------------

# Sensitivity: fraction of true signal features recovered among the first nx_sub / ny_sub indices
sensitivity_x  <- numeric(num_episods)
sensitivity_y  <- numeric(num_episods)
sensitivity_x2 <- numeric(num_episods)  # after Phase-1 (2d)
sensitivity_y2 <- numeric(num_episods)
sensitivity_x4 <- numeric(num_episods)  # after Phase-1 (4d)
sensitivity_y4 <- numeric(num_episods)

for (k in 1:num_episods) {
  
  sensitivity_x[k]  <- sum(results[[k]]$x %in% 1:nx_sub) / nx_sub
  sensitivity_y[k]  <- sum(results[[k]]$y %in% 1:ny_sub) / ny_sub
  sensitivity_x2[k] <- sum(results2[[k]]$x %in% 1:nx_sub) / nx_sub
  sensitivity_y2[k] <- sum(results2[[k]]$y %in% 1:ny_sub) / ny_sub
  sensitivity_x4[k] <- sum(results4[[k]]$x %in% 1:nx_sub) / nx_sub
  sensitivity_y4[k] <- sum(results4[[k]]$y %in% 1:ny_sub) / ny_sub
}

# Precision: fraction of selected features that are truly in the first nx_sub / ny_sub
precision_x  <- numeric(num_episods)
precision_y  <- numeric(num_episods)
precision_x2 <- numeric(num_episods)
precision_y2 <- numeric(num_episods)
precision_x4 <- numeric(num_episods)
precision_y4 <- numeric(num_episods)

for (k in 1:num_episods) {
  precision_x[k]  <- sum((1:nx_sub) %in% results[[k]]$x) / length(results[[k]]$x)
  precision_y[k]  <- sum((1:ny_sub) %in% results[[k]]$y) / length(results[[k]]$y)
  precision_x2[k] <- sum((1:nx_sub) %in% results2[[k]]$x) / length(results2[[k]]$x)
  precision_y2[k] <- sum((1:ny_sub) %in% results2[[k]]$y) / length(results2[[k]]$y)
  precision_x4[k] <- sum((1:nx_sub) %in% results4[[k]]$x) / length(results4[[k]]$x)
  precision_y4[k] <- sum((1:ny_sub) %in% results4[[k]]$y) / length(results4[[k]]$y)
}

# NOTE: The following blocks compare to baselines (dcor/bcor/pcor) that are NOT defined above.
# Keep them only if you compute result_dcor/result_bcor/result_pcor elsewhere.
# precision_dcor2 <- precision_bcor2 <- precision_pcor2 <- numeric(num_episods)
# precision_dcor4 <- precision_bcor4 <- precision_pcor4 <- numeric(num_episods)
# for (k in 1:num_episods) {
#   x_index <- sort(order(result_dcor[[k]])[(ncol(X) - 2*d + 1):ncol(X)])
#   precision_dcor2[k] <- sum((1:nx_sub) %in% x_index) / length(x_index)
#   ...
# }

# F1 scores
f1_x  <- 2 * sensitivity_x  * precision_x  / (sensitivity_x  + precision_x)
f1_y  <- 2 * sensitivity_y  * precision_y  / (sensitivity_y  + precision_y)
f1_x2 <- 2 * sensitivity_x2 * precision_x2 / (sensitivity_x2 + precision_x2)
f1_y2 <- 2 * sensitivity_y2 * precision_y2 / (sensitivity_y2 + precision_y2)
f1_x4 <- 2 * sensitivity_x4 * precision_x4 / (sensitivity_x4 + precision_x4)
f1_y4 <- 2 * sensitivity_y4 * precision_y4 / (sensitivity_y4 + precision_y4)

## Report means across replicates
mean(sensitivity_x);  mean(sensitivity_y)
mean(sensitivity_x2); mean(sensitivity_y2)
mean(sensitivity_x4); mean(sensitivity_y4)

mean(precision_x);  mean(precision_y)
mean(precision_x2); mean(precision_y2)
mean(precision_x4); mean(precision_y4)

mean(f1_x);  mean(f1_y)
mean(f1_x2); mean(f1_y2)
mean(f1_x4); mean(f1_y4)