group_matrices <- function(y,X,grouping=2, M=25e3, sim_type="mc" ) {
  N <- length(y)
  groups <- (1 + (1:N %%  (N/grouping)))
  groups <- groups[order(groups)]
  groups <- as.integer(groups)
  groups <- list( y=groups )
  groups $ Sigma <- Matrix::bdiag(replicate(ceiling(N/grouping),matrix(1,grouping,grouping), simplify=F))[1:length(groups$y),1:length(groups$y)]
  groups $ Sigma <- as.matrix(diag(groups$y) %*%  groups$Sigma)
  groups $ X <- groups$X <- t(apply(matrix(groups$y),1, function(x) rep(x,ncol(X))))
  table <- table(groups$y)-1
  if (sim_type=="mc") {
    rdm <- runif(M*sum(table))
    rdm <- matrix(rdm, M)
  } else {
    rdm <- qrng::ghalton(M,sum(table),method =  "generalized")
    }
  splitmat <- rep(as.integer(names(table)),table) #questa va replicata M volte
  splitmat <- as.matrix(splitmat) #questa va replicata M volte (ogni riga)
  rdm <- split(rdm, splitmat)
  rdm <- lapply(rdm, function(x) matrix(x, ncol=M))
  groups$rdm <- t(rdm)
  if (length(groups$rdm)!= max(groups$y))
    length(groups$rdm) <- max(groups$y)
  groups
}
################################f_sar
f_sar <- function (rho_tilde,W,X,eig,qu=Inf,method_inv='solve') {
  rho  <- to_natural(rho_tilde, eig$W )
  n <- dim(W)[1]

  if (rho==0) {
    Inv_rho <- Sigma_rho <- diag(n)
    X_rho <- X
    dotSigma_rho <- W+t(W)
    dotX_rho <- W%*%X
    return(list(Inv_rho=Inv_rho,
                X_rho=X_rho,
                Sigma_rho=Sigma_rho,
                dotSigma_rho=dotSigma_rho,dotX_rho=dotX_rho))
  }

  # if (method_inv=="solve") {
  #   A <- as(diag(n)-rho*W,"CsparseMatrix")
  #   Inv_rho <-  solve(A, tol=tol.solve)
  #   #Inv_rho <-  solve(as(diag(n)-rho*W,"Matrix"), tol=tol.solve)
  #   # Sigma_rho <- tcrossprod(Inv_rho)
  # }
  # else {
  #   if (method_inv=="chol") {
  #     A <- as(diag(n)-rho*W,"CsparseMatrix")
  #     Inv_rho <- Matrix::chol2inv( chol(A))
  #     # Sigma_rho <- tcrossprod(Inv_rho)
  #     # Inv_Sigma <- crossprod(diag(n)-rho*W)
  #     # Sigma_rho <-  Matrix::chol2inv(chol(as(Inv_Sigma,"CsparseMatrix")))
  #     # Inv_rho <- tcrossprod(Sigma_rho,diag(n)-rho*W)
  #   }
  #   else {
  #     Inv_rho <- addendum <-  diag(n)
  #     for (k in 1:qu) #thanks fagiant
  #     {
  #       addendum <- rho*W%*%addendum
  #       Inv_rho <- Inv_rho + addendum
  #     }
  #
  #   }
  # }
  if ((method_inv %in% c('solve','chol'))) {

    Inv_rho <- spatialreg::invIrW(x=W, rho=rho, method=method_inv ,feasible=T)
  } else {
    Inv_rho <- addendum <-  diag(n)
    for (k in 1:qu) #thanks fagiant
    {
      addendum <- rho*W%*%addendum
      Inv_rho <- Inv_rho + addendum
    }
  }

        Sigma_rho <- tcrossprod(Inv_rho)
        Sigma_rho <- .5*(Sigma_rho + t(Sigma_rho))
#  assign("my_result", list(Sigma_rho=Sigma_rho, rho=rho), envir = globalenv())
  jacob <- ((eig$W[2]^(-1)-eig$W[1]^(-1))*exp(-rho_tilde))/((1+exp(-rho_tilde))^2) ##derivative of rho wrt rho tilde
  dotSigma_rho <- jacob * Inv_rho %*% W %*% Sigma_rho #now this is expressed as a function of rho tilde (chain rule)
  dotSigma_rho <- dotSigma_rho + t(dotSigma_rho)
    X_rho <- Inv_rho%*%X
  dotX_rho <- jacob *Inv_rho%*%W%*%X_rho

  list(Inv_rho=Inv_rho, X_rho=X_rho, Sigma_rho=Sigma_rho, dotSigma_rho=dotSigma_rho,dotX_rho=dotX_rho)

}

logLIK_SAR <- function(theta,y,W,X,eig,qu=Inf,method_inv="solve", groups,
                       mvtnorm_control=list(M=25e3, tol = .Machine$double.eps, fast = FALSE), bobyqa=F) {

  beta <- theta[ - length(theta)  ]
  rho_tilde <- tail(theta,1)

  f_rho <-  f_sar(rho_tilde,W,X,eig,qu,method_inv)
  Sigma <- f_rho$Sigma_rho
  xb <- f_rho$X_rho%*%beta
  y_list <- split(y,groups$y)
  dotXb_rho <- f_rho$dotX_rho%*%beta
  dotXb_rho_list <- split(dotXb_rho,groups$y)

  X_rho_list <- lapply(split( f_rho$X_rho, groups$X), function(x) matrix(x, ncol=ncol(X)) )
  xb_list <- split(xb, groups$y)
  assign("theta", theta, envir = globalenv())
  Sigma_list <- lapply(split(f_rho$Sigma_rho, groups$Sigma)[-1],
                       function(x) matrix(x, ncol=sqrt(length(x))))
  dotSigma_list <- lapply(split(f_rho$dotSigma_rho, groups$Sigma)[-1],
                          function(x) matrix(x, ncol=sqrt(length(x))))
  #ll <-  logLIK_g(xb_list[[1]], Sigma_list[[1]], y_list[[1]],X_rho_list[[1]],dotSigma_list[[1]])
  ll <- t(mapply(logLIK_g_SAR, xb_list, Sigma_list, y_list,X_rho_list,dotXb_rho_list,dotSigma_list,groups$rdm,
               MoreArgs = list(mvtnorm_control=mvtnorm_control),  SIMPLIFY =T))
  score <- ll[,-1]
  ll <- ll[,1]
  attr(ll,"f")  <- f_rho
    if (!bobyqa) {
    attr(ll,"gradient") <-  -score
    return(ll)
  } else
    return(-sum(ll))

}

logLIK_g_SAR <- function(xb,Sigma,y,  X_rho ,dotXb_rho, dotSigma, simw,
                          mvtnorm_control=list(M=25e3, tol = .Machine$double.eps, fast = FALSE)  ) {
  n <- length(y)
  upper <- rep(Inf, n  )
  lower <- -upper
  upper[y==0] <- -xb[y==0]
  lower[y==1] <- -xb[y==1]
  Sigma <- .5*(Sigma  + t(Sigma) )
  assign("Sigma", Sigma, envir = globalenv())
  # if (!matrixcalc::is.positive.definite(Sigma, tol=1e-8)) {
  #   #Sigma <- Matrix::nearPD(Sigma, ensureSymmetry = T)$mat
  #   assign("my_result", Sigma, envir = globalenv())
  #   stopifnot(matrixcalc::is.positive.definite(Sigma, tol=1e-8))
  #   }
  C <- t(chol(Sigma))
  chol <- mvtnorm::ltMatrices(C[lower.tri(C, diag=T)], diag = T)
  p <- mvtnorm::slpmvnorm(lower=lower, upper=upper, chol=chol, M=mvtnorm_control$M, w=simw,
                        tol=mvtnorm_control$tol, fast=mvtnorm_control$fast, logLik = TRUE)
  score_beta <- crossprod(-p$mean, X_rho)
  score_rho1 <-  crossprod(-p$mean, dotXb_rho)
  if (n>1) {
    #dSigma_dC <-  elimination.matrix(n)%*% ( (C%x%diag(n)) + (diag(n)%x%C)%*%commutation.matrix(n))%*%duplication.matrix(n)
    L <- matrixcalc::elimination.matrix(n)
    dSigma_dC <- L %*% (diag(n^2) + matrixcalc::commutation.matrix(n)) %*% (C %x% diag(n)) %*% t(L)
    dC_dSigma <-  solve(dSigma_dC)
  } else dC_dSigma <- as.matrix(0.5/C)

  score_rho <- t(mvtnorm::Lower_tri(p$chol, diag = T))%*%dC_dSigma%*%dotSigma[lower.tri(dotSigma,diag=T)]
  c(p$logLik, score_beta,score_rho1-score_rho)

}


