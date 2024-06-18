null_or_val <- function(x) {

  if (is.null(x)) {

    return(0)

  } else {

    return(x)

  }

}

#' @export
#'
#' @title scGBM: Model-based dimensionality reduction for single-cell
#' RNA-seq with generalized bilinear models
#'
#' @description Fit a Poisson bilinear model to the single-cell
#' count matrix.
#'
#' @param Y A matrix of UMI counts with genes on rows and cells on columns.
#' @param M The number of latent factors to estimate
#' @param init_U Initialization for U
#' @param init_V Initialization for V
#' @param init_d Initialization for d
#' @param max.iter The maximum number of iterations
#' @param max_seconds_run Maximum number of seconds allowed to fit model.
#' @param tol The tolerance for convergence (relative difference in objective function)
#' @param subset If NULL, use the entire data to fit the model. If an integer, take a
#' random subsample of cells to fit the loadings. Then project the remaining
#' cells onto those loadings. If a non-random sample is desired, subset can
#' be an integer vector corresponding to the columns that should be used to
#' fit the loadings.
#' @param ncores If subset is not NULL, then ncores specifies the number of cores to use
#' during the projection step.
#' @param infer.beta If FALSE, beta is set to the log of the number of counts per cell.
#' @param return.W If TRUE, the matrix of weights (which is equal to the estimated mean) is returned. This should
#' be FALSE for large datasets since it requires a lot of memory.
#' @param batch An optional factor containing the assignment of cells to known batches.
#'
#' @return A list with components
#' \itemize{
#' \item \code{V} - A matrix containing the factor scores.
#' \item \code{U} - A matrix contianing the factor loadings
#' \item \code{D} - A vector containing the singular values (scaling factors)
#' \item \code{alpha} - A matrix containing gene-specific intercepts. The number of columns
#' is set to be equal to the number of batches (by default there are no batches so this is 1).
#' \item \code{beta} - A vector of cell-specific intercepts.
#' \item \code{I} - The number of genes
#' \item \code{J} - The number of cells.
#' \item \code{W} - The estimated mean (and by properties of the Poisson distribution, also the variance)
#' for each entry of the count matrix.
#' \item \code{obj} The value of the objective function for each iteration.
#' }
#'
#' @author Phillip B. Nicol <philnicol740@gmail.com>
gbm.sc <- function(Y,
                   M,
                   cluster = FALSE,
                   celltype = NULL,
                   current_iter = NULL,
                   LL = NULL,
                   init_U = NULL,
                   init_V = NULL,
                   init_d = NULL,
                   max.iter=100,
                   max_seconds_run = 10 * 3600,
                   tol=10^{-4},
                   subset=NULL,
                   ncores=1,
                   infer.beta=FALSE,
                   return.W = TRUE,
                   batch=as.factor(rep(1,ncol(Y))),
                   time.by.iter = FALSE,
                   lr=1,
                   min.iter=30) {
  if(!is.null(subset)) {
    out <- gbm.proj.parallel(Y,M,subsample=subset,ncores=ncores,tol=tol,
                             max.iter=max.iter)
    return(out)
  }


  I <- nrow(Y); J <- ncol(Y)
  if (is.null(LL)) {

    LL <- rep(0,max.iter)

  }

  loglik <- c()
  if(time.by.iter) {
    time <- c()
  }

  if(!is.factor(batch)) {
    stop("Batch must be encoded as factor.")
  }
  batch.factor <- batch
  batch <- as.integer(batch)
  nbatch <- max(batch)

  #Precompute relevant quantities
  max.Y <- max(Y)
  nz <- which(Y != 0)

  total_model_time <- 0

  init_start_time <- Sys.time()
  #Starting estimate for alpha and W
  betas <- log(colSums(Y))
  #betas <- betas - mean(betas) Enforce the betas sum to 0
  W <- matrix(0, nrow=I, ncol=J)
  log.rsy <- matrix(0,nrow=nbatch,ncol=I)
  for(j in 1:nbatch) {
    log.rsy[j,] <- log(rowSums(Y[,batch==j]))
  }
  alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
    log.rsy[j,]-log(sum(exp(betas[batch==j])))
  })
  betas <- betas + mean(alphas)
  alphas <- alphas - mean(alphas)

  if (is.null(init_U) || is.null(init_V) || is.null(init_d)) {

    W <- exp(sweep(alphas[,batch], 2, betas, "+"))

    #Starting estimate of X
    Z <- (Y-W)/sqrt(W)
    c <- sqrt(2*log(I*J/0.025))
    #Z[Z > c] <- c
    #Z[Z < -c] <- -c
    LRA <-  irlba::irlba(Z,nv=M,nu=M)
    X <- LRA$u %*%(LRA$d*t(LRA$v))

  } else {

    X <- init_U %*%(init_d*t(init_V))

  }

  X <- sqrt(1/W)*X

  X[X > 8] <- 8
  X[X < -8] <- -8

  #For acceleration, save previous X
  Xt <- matrix(0,nrow=I,ncol=J)

  init_end_time <- Sys.time()
  total_model_time <- total_model_time + as.numeric(
    difftime(
      init_end_time, init_start_time, units = "secs"
    )
  )

  nmi_vec <- numeric(max.iter)
  ari_vec <- numeric(max.iter)

  for(i in 1:max.iter) {

    if (cluster) {

      print("clustering...")
      d <- distances::distances(LRA$v %*% diag(LRA$d))
      dm <- distances::distance_matrix(d)
      clust_tree <- fastcluster::hclust(dm, method="ward.D2")
      clusts <- cutree(clust_tree, k = 10)
      nmi <- aricode::NMI(celltype, clusts)
      ari <- aricode::ARI(celltype, clusts)
      nmi_vec[i] <- nmi
      ari_vec[i] <- ari
      print(nmi)
      print(ari)

    }

    iter_start_time <- Sys.time()

    #Reweight
    alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
      #sweep(X[,batch==j],2,betas[batch==j],"+")
      log.rsy[j,]-log(rowSums(exp(sweep(X[,batch==j],2,betas[batch==j],"+"))))
    })
    betas <- betas + mean(alphas)
    alphas <- alphas - mean(alphas)
    if(infer.beta) {
      betas <- log(colSums(Y))-log(colSums(exp(alphas[,batch]+X)))
    }
    W <- exp(sweep(alphas[,batch]+X, 2, betas, "+"))

    #Prevent W from being too large (stability)
    W[W > max.Y] <- max.Y

    #Compute working variables
    #Z <- X+(Y-W)/W

    ## Compute log likelihood (no normalizing constant)
    LL_i_1 <- tail(LL, 1)
    LLi <- sum(Y[nz]*log(W[nz]))-sum(W)
    LL <- c(LL, LLi)
    if(is.na(LL[i]) | is.infinite(LL[i])) {
      X <- Xt
      lr <- lr/2
      i <- i - 1
      LL <- LL[1:(length(LL) - 1)]
      #next
    }
    if(i >= 3 || (!is.null(current_iter) & null_or_val(current_iter) >= 3)) {
      #tau <- abs((LL[i]-LL[i-2])/LL[i])
      # if(tau < tol & lr <= 1.06 & i >= min.iter) {
      #  break
      # }

      if(LLi <= (LL_i_1+0.1)) {
        lr <- max(lr/2, 1)
        X <- Xt
        #next
      } else {
        lr <- lr*(1.05)
        #lr <- lr
      }
    }

    loglik <- c(loglik,LLi)
    cat("Iteration: ", i, ". Objective=", LLi, "\n")

    ## Gradient Step
    V <- X+((i-1)/(i+2))*(X-Xt)
    Xt <- X
    w.max <- max(W)

    LRA <- irlba::irlba(V+(lr/w.max)*(Y-W),nv=M)
    X <- LRA$u %*%(LRA$d*t(LRA$v))

    if(i == max.iter) {
      warning("Maximum number of iterations reached (increase max.iter).
              Possible non-convergence.")
    }

    end_iter_time <- Sys.time()
    iter_time <- as.numeric(
      difftime(
        end_iter_time, iter_start_time, units = "secs"
        )
      )

    if(time.by.iter) {
      time <- c(time,iter_time)
    }

    total_model_time <- total_model_time + iter_time

    if (total_model_time >= max_seconds_run) {

      warning(
        sprintf("Algorithm reached maximum time without convergence.")
      )
      break

    }


  }

  out <- list()
  if(return.W) {
    out$W <- t(t(exp(alphas[,batch]+X))*exp(betas))
  }
  out$V <- LRA$v; rownames(out$V) <- colnames(Y); colnames(out$V) <- 1:M
  out$D <- LRA$d; names(out$D) <- 1:M
  out$U <- LRA$u; rownames(out$U) <- rownames(Y); colnames(out$U) <- 1:M
  out$alpha <- alphas; rownames(out$alpha) <- rownames(Y)
  out$beta <- betas; names(out$beta) <- colnames(Y)
  out$M <- M
  out$I <- nrow(out$W); out$J <- ncol(out$W)
  out$loglik <- loglik
  out$alpha <- drop(out$alpha)
  out$beta <- drop(out$beta)
  out$V <- out$V %*% diag(1/out$D)
  if(time.by.iter)
    out$time <- cumsum(time)
  out <- process.results(out)
  out$LL <- LL
  out$lr <- lr
  out$ari <- ari_vec
  out$nmi <- nmi_vec
  return(out)
}

#' Initialize scGBM fit with PCA on Rank-2 Approximation
#'
#' @param Y data
#' @param M rank of approximation
#' @param batch batch for each variable
#'
#' @return fit
#' @export
#'
gbm.init <- function(Y, M, batch=as.factor(rep(1,ncol(Y)))) {

  I <- nrow(Y); J <- ncol(Y)

  if(!is.factor(batch)) {
    stop("Batch must be encoded as factor.")
  }
  batch.factor <- batch
  batch <- as.integer(batch)
  nbatch <- max(batch)

  #Precompute relevant quantities
  max.Y <- max(Y)
  nz <- which(Y != 0)

  full_model_start_time <- Sys.time()

  #Starting estimate for alpha and W
  betas <- log(colSums(Y))
  #betas <- betas - mean(betas) Enforce the betas sum to 0
  W <- matrix(0, nrow=I, ncol=J)
  log.rsy <- matrix(0,nrow=nbatch,ncol=I)
  for(j in 1:nbatch) {
    log.rsy[j,] <- log(rowSums(Y[,batch==j]))
  }
  alphas <- vapply(1:nbatch, FUN.VALUE=numeric(I), function(j) {
    log.rsy[j,]-log(sum(exp(betas[batch==j])))
  })
  betas <- betas + mean(alphas)
  alphas <- alphas - mean(alphas)
  W <- exp(sweep(alphas[,batch], 2, betas, "+"))

  #Starting estimate of X
  Z <- (Y-W)/sqrt(W)
  c <- sqrt(2*log(I*J/0.025))
  #Z[Z > c] <- c
  #Z[Z < -c] <- -c
  LRA <-  irlba::irlba(Z,nv=M,nu=M)

  out <- list()

  out$V <- LRA$v; rownames(out$V) <- colnames(Y); colnames(out$V) <- 1:M
  out$D <- LRA$d; names(out$D) <- 1:M
  out$U <- LRA$u; rownames(out$U) <- rownames(Y); colnames(out$U) <- 1:M

  return(out)

}

gbm.proj.parallel <- function(Y,M,subsample=2000,min.counts=5,
                              ncores,tol=10^{-4},max.iter=max.iter) {

  J <- ncol(Y); I <- nrow(Y)
  alphas.full <- log(rowSums(Y))
  if(length(subsample)==1) {
    jxs <- sample(1:J,size=subsample,replace=FALSE)
    Y.sub <- Y[,jxs]
  } else {
    Y.sub <- Y[,subsample]
  }
  Y.sub <- as.matrix(Y.sub)
  ixs <- which(rowSums(Y.sub) > 5)
  Y.sub <- Y.sub[ixs,]
  out <- gbm.sc(Y.sub,M=M,tol=tol,max.iter=max.iter)

  U <- out$U
  #U <- as.data.frame(U)
  U <- cbind(rep(1,length(ixs)),U)
  colnames(U) <- c("intercept",paste0("U",1:M))
  #alpha <- alphas.full
  alpha <- out$alpha[,1]
  #alpha <- alphas.full[ixs]

  Y <- Y[ixs,]
  V <- matrix(0,nrow=J,ncol=2*M+1)
  split.len <- 200000
  max.iter <- ceiling(J/split.len)
  for(i in 1:max.iter) {
    start <- (i-1)*split.len + 1
    if(i == max.iter) {
      stop <- J
    } else {
      stop <- i*split.len
    }
    V.sub <- parallel::mclapply(start:stop, function(j) {
      cell <- as.vector(Y[,j])
      o <- alpha#+log(sum(cell))
      val <- matrix(rep(0,2*M+1),nrow=1)
      try({
        fit <- fastglm::fastglm(x=U,y=cell,offset=o,
                       family=poisson(),
                       method=3)
        #fit <- glm(cell~0+offset(o)+.,
        #data=U,
        #family=poisson(link="log"))
        val <- matrix(c(fit$coefficients,fit$se[-1]),
                      nrow=1)
      })
      val
    },
    mc.cores = ncores)
    V.sub <- matrix(unlist(V.sub),nrow=stop-start+1,ncol=2*M+1,byrow=TRUE)
    V[start:stop,] <- V.sub
  }
  out$se_V <- V[,(M+2):(2*M+1)]
  out$V <- V[,2:(M+1)]
  out$beta <- V[,1]

  alpha <- rep(0, I)
  alpha[ixs] <- out$alpha[,1]
  alpha[-ixs] <- alphas.full[-ixs] - log(sum(exp(out$beta)))
  out$beta <- out$beta + mean(alpha)
  alpha <- alpha - mean(alpha)
  out$alpha <- alpha
  U <- matrix(0, nrow=I,ncol=M)
  U[ixs,] <- out$U
  out$U <- U
  return(out)
}

process.results <- function(gbm) {
  #Enforce identifiability in U
  M <- gbm$M
  for(m in 1:M) {
    if(gbm$U[1,m] < 0) {
      gbm$U[,m] <- -1*gbm$U[,m]
      gbm$V[,m] <- -1*gbm$V[,m]
    }
  }

  gbm$V <- t(diag(gbm$D) %*% t(gbm$V))

  return(gbm)
}
