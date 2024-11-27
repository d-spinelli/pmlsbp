f_sarar <- function (rho_tilde,lambda_tilde,W,M,X,eig,qu=Inf,method_inv="solve") {
  rho  <- eig$W[1]^(-1) + ((eig$W[2]^(-1)-eig$W[1]^(-1))/(1+exp(-rho_tilde)))
  lambda  <- eig$M[1]^(-1) + ((eig$M[2]^(-1)-eig$M[1]^(-1))/(1+exp(-lambda_tilde)))
  n<-dim(W)[1]
  if (method_inv=="solve") {
    Inv_rho <- solve(diag(n)-rho*W)
    Inv_lambda <- solve(diag(n)-lambda*M)
    Sigma_rho <- tcrossprod(Inv_rho%*%Inv_lambda)
  }
  else {
    if (method_inv=="chol") {
      Inv_Sigma <- crossprod(diag(n)-rho*W-lambda*M+rho*lambda*M%*%W)
      Sigma_rho <- chol2inv(chol(as(Inv_Sigma,"CsparseMatrix")))
      S_rho <- crossprod(diag(n)-rho*W)
      Inv_rho <- chol2inv(chol(as(S_rho,"CsparseMatrix")))%*%(diag(n)-rho*W)
      S_lambda <- crossprod(diag(n)-lambda*M)
      Inv_lambda <- chol2inv(chol(as(S_lambda,"CsparseMatrix")))%*%(diag(n)-lambda*M)

    }
    else {
      Inv_rho <- Inv_lambda <- addendum_rho <- addendum_lambda  <-  diag(n)
      for (k in 1:qu) #thanks fagiant
      {
        addendum_rho <- rho*W%*%addendum_rho
        Inv_rho <- Inv_rho + addendum_rho
        addendum_lambda <- lambda*M%*%addendum_lambda
        Inv_lambda <- Inv_lambda + addendum_lambda
      }
      Sigma_rho <- tcrossprod(Inv_rho%*%Inv_lambda)
    }
  }
  dotSigma_rho <- Inv_rho %*% W %*% Sigma_rho
  dotSigma_rho <- dotSigma_rho + t(dotSigma_rho)
  dotSigma_lambda=matrix(0,n,n)
  if (lambda==0) {
    dotSigma_lambda=Inv_rho%*%(M+t(M))%*%t(Inv_rho)
  } else {
    dotSigma_lambda <- Inv_rho%*%Inv_lambda%*%M%*%Inv_lambda%*%t(Inv_lambda)%*%t(Inv_rho)
    dotSigma_lambda <- dotSigma_lambda + t(dotSigma_lambda)
  }
  X_rho <- Inv_rho%*%X

  jacob_rho <- ((eig$W[2]^(-1)-eig$W[1]^(-1))*exp(-rho_tilde))/((1+exp(-rho_tilde))^2)
  dotSigma_rho <- jacob_rho*dotSigma_rho
  dotX_rho <- jacob_rho*Inv_rho%*%W%*%X_rho
  jacob_lambda <- ((eig$M[2]^(-1)-eig$M[1]^(-1))*exp(-lambda_tilde))/((1+exp(-lambda_tilde))^2)
  dotSigma_lambda <- jacob_lambda*dotSigma_lambda

  list(Inv_rho=Inv_rho, X_rho=X_rho, Sigma_rho=Sigma_rho,
       dotSigma_rho=dotSigma_rho,dotSigma_lambda=dotSigma_lambda,
       Inv_lambda=Inv_lambda,dotX_rho=dotX_rho)

}

#########

logLIK_SARAR <- function(theta,y,W,M,X,eig,qu=Inf,method_inv="solve", groups,
                         mvtnorm_control=list(M=25e3, w=NULL, tol = .Machine$double.eps, fast = FALSE), bobyqa=F) {
  beta <- theta[ 1:(length(theta)-2)  ]
  lambda_tilde <- tail(theta,1)
  rho_tilde <-  theta[length(theta)-1]
  f_rho<- f_sarar(rho_tilde,lambda_tilde,W,M,X,eig,qu,method_inv)
  Sigma <- f_rho$Sigma_rho
  xb <- f_rho$X_rho%*%beta
  y_list<-split(y,groups$y)
  dotXb_rho<-f_rho$dotX_rho%*%beta
  dotXb_rho_list<-split(dotXb_rho,groups$y)

  X_rho_list <- lapply(split( f_rho$X_rho, groups$X), function(x) matrix(x, ncol=ncol(X)) )
  xb_list <- split(xb, groups$y)

  Sigma_list <- lapply(split(f_rho$Sigma_rho, groups$Sigma)[-1],
                       function(x) matrix(x, ncol=sqrt(length(x))))
  dotSigma_rho_list <- lapply(split(f_rho$dotSigma_rho, groups$Sigma)[-1],
                          function(x) matrix(x, ncol=sqrt(length(x))))
  dotSigma_lambda_list <- lapply(split(f_rho$dotSigma_lambda, groups$Sigma)[-1],
                          function(x) matrix(x, ncol=sqrt(length(x))))
  # ll<- logLIK_g_SARAR(xb_list[[1]], Sigma_list[[1]], y_list[[1]],X_rho_list[[1]],dotSigma_rho_list[[1]],dotSigma_lambda_list[[1]])
  ll<-t(mapply(logLIK_g_SARAR, xb_list, Sigma_list, y_list,X_rho_list, dotXb_rho_list ,dotSigma_rho_list,
                dotSigma_lambda_list, groups$rdm,
               MoreArgs = list(mvtnorm_control=mvtnorm_control),
                SIMPLIFY =T))
  score<-ll[,-1]
  ll<-ll[,1]
  attr(ll,"f") <-f_rho
  attr(ll,"gradient")<- -score
  if (!bobyqa) {
    attr(ll,"gradient")<- -score
    return(ll)
  } else
    return(-sum(ll))
}


logLIK_g_SARAR <- function(xb, Sigma, y,  X_rho, dotXb_rho, dotSigma_rho, dotSigma_lambda, simw,
                           mvtnorm_control=list(M=25e3,tol = .Machine$double.eps, fast = FALSE)  ) {
  upper <- rep(Inf, length(y)  )
  lower <- -upper
  n <- length(y)
  upper[y==0] <- -xb[y==0]
  lower[y==1] <- -xb[y==1]
  C<-t(chol(Sigma))
  chol<-mvtnorm::ltMatrices(C[lower.tri(C, diag=T)], diag = T)
  p<-mvtnorm::slpmvnorm(lower=lower, upper=upper, chol=chol, M=mvtnorm_control$M, w=simw, tol=mvtnorm_control$tol, fast=mvtnorm_control$fast)
  score_beta <- crossprod(-p$mean, X_rho)
  score_rho1<- crossprod(-p$mean, dotXb_rho)
  if (n>1) {
    L <- matrixcalc::elimination.matrix(n)
    dSigma_dC <- L %*% (diag(n^2) + matrixcalc::commutation.matrix(n)) %*% (C %x% diag(n)) %*% t(L)
    dC_dSigma<- solve(dSigma_dC) } else dC_dSigma <- as.matrix(0.5/C)
  score_rho <- t(mvtnorm::Lower_tri(p$chol, diag = T))%*%dC_dSigma%*%dotSigma_rho[lower.tri(dotSigma_rho,diag=T)]
  score_lambda <- t(mvtnorm::Lower_tri(p$chol, diag = T))%*%dC_dSigma%*%dotSigma_lambda[lower.tri(dotSigma_lambda,diag=T)]
  c(p$logLik, score_beta,score_rho1-score_rho,-score_lambda)
}

#################### bobyqa

#
# logLIK_SARAR_bobyqa <- function(theta,y,W,M,X,eig,qu=Inf,method_inv="solve", groups, mvtnorm_control=list(M=25e3, w=NULL, tol = .Machine$double.eps, fast = FALSE)) {
#   beta <- theta[ 1:(length(theta)-2)  ]
#   lambda_tilde <- tail(theta,1)
#   rho_tilde <-  theta[length(theta)-1]
#   f_rho<- f_sarar(rho_tilde,lambda_tilde,W,M,X,eig,qu,method_inv)
#   Sigma <- f_rho$Sigma_rho
#   xb <- f_rho$X_rho%*%beta
#   y_list<-split(y,groups$y)
#   dotXb_rho<-f_rho$dotX_rho%*%beta
#   dotXb_rho_list<-split(dotXb_rho,groups$y)
#
#   X_rho_list <- lapply(split( f_rho$X_rho, groups$X), function(x) matrix(x, ncol=ncol(X)) )
#   xb_list <- split(xb, groups$y)
#
#   Sigma_list <- lapply(split(f_rho$Sigma_rho, groups$Sigma)[-1],
#                        function(x) matrix(x, ncol=sqrt(length(x))))
#   dotSigma_rho_list <- lapply(split(f_rho$dotSigma_rho, groups$Sigma)[-1],
#                               function(x) matrix(x, ncol=sqrt(length(x))))
#   dotSigma_lambda_list <- lapply(split(f_rho$dotSigma_lambda, groups$Sigma)[-1],
#                                  function(x) matrix(x, ncol=sqrt(length(x))))
#   ll<-t(mapply(logLIK_g_SARAR_bobyqa, xb_list, Sigma_list, y_list,X_rho_list,  groups$rdm,
#                MoreArgs = list(mvtnorm_control=mvtnorm_control),
#                SIMPLIFY =T))
#   sum(ll)
#
# }
#
#
# logLIK_g_SARAR_bobyqa <- function(xb, Sigma, y,  X_rho,  simw,
#                            mvtnorm_control=list(M=25e3,tol = .Machine$double.eps, fast = FALSE)  ) {
#   upper <- rep(Inf, length(y)  )
#   lower <- -upper
#   n <- length(y)
#   upper[y==0] <- -xb[y==0]
#   lower[y==1] <- -xb[y==1]
#   C<-t(chol(Sigma))
#   chol<-mvtnorm::ltMatrices(C[lower.tri(C, diag=T)], diag = T)
#   p<-mvtnorm::slpmvnorm(lower=lower, upper=upper, chol=chol, M=mvtnorm_control$M, w=simw, tol=mvtnorm_control$tol, fast=mvtnorm_control$fast)
#   p$logLik
# }
