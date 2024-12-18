#' Predicted values based on pmlsprobit object
#'
#' @param object a \code{"pmlsprobit"} object
#' @param type \code{"prob"} for the probability \eqn{P(Y_i=1)} or \code{"me"}
#' for the direct, indirect and total average marginal effect for each observation
#' @param newdata 	An optional data frame in which to look for variables with which to predict.
#' If omitted, the fitted values are used.
#' @param variables a string including the variable names of the covariates for
#' which the marginal effects should be calculated using numeric derivatives.  It is valid only if \code{type="me"}.
#' This option should be used only if the model contains non-linear terms in the specification \eqn{X \beta}.
#' @param delta.method a boolean. If \code{TRUE} the standard errors of the predictions are calculated.
#' It is valid only if \code{type="me"}
#'
#' @return
#' @export
#'
#' @examples
predict.pmlsprobit<- function(object,
                                 type='prob',
                                 newdata=NULL, variables=NULL,delta.method=F,...) {
  if(!inherits(object, "pmlsprobit"))
    stop("'predict.pmlsprobit' called on a non-'pmlsprobit' object")
stopifnot(type %in% c('prob','response','margins','xb','me','mfx','p'))
if (type %in% c('prob','response','p'))
  ret <-  predict.prob.pmlsprobit(object, newdata)[,1]
if (type=='xb')
  ret <- predict.xb.pmlsprobit(object, newdata )[,1]
if (type %in% c('me','margins','mfx')) {
ret <- predict.me.pmlsprobit(object, variables = variables, delta.method=delta.method)
class(ret) <- "pe.pmlsprobit"
}
  ret
  }

predict.me.pmlsprobit<- function(mod, variables=NULL, delta.method=T) {
  variables.factor <- names(which(attr(mod$terms,"dataClasses")[-1]=="factor"))
  variables.numeric <- names(which(attr(mod$terms,"dataClasses")[-1]=="numeric"))
  passthru.factors <- passthru.numeric <- FALSE
  if (!is.null(variables)) {
    if ( any(variables.factor %in% variables) )
      variables.factor <-  variables.factor[variables.factor %in% variables]
    else
      passthru.factors <- TRUE

    if ( any(variables.numeric %in% variables) )
      variables.numeric <-  variables.numeric[variables.numeric %in% variables]
    else
      passthru.numeric <- TRUE
  } else {
    passthru.numeric <- length(variables.numeric) == 0
    passthru.factors <- length(variables.factor) == 0
    variables.numeric <- variables.factor <- NULL
  }
  if (passthru.factors*passthru.numeric)
    stop('the specified variables were not used to estimate the model')
  me.numeric <- me.factor <- NULL
  if (!passthru.numeric) {
    me.numeric <- predict.me.numeric.pmlsprobit.v2(mod, variables.numeric, delta.method)
    if (passthru.factors)
      return(me.numeric)
  }
  if (!passthru.factors) {
    me.factor <- predict.me.factor.pmlsprobit.v2(mod, variables.factor, delta.method )
    if (passthru.numeric)
      return(me.factor)
  }
  f <- abind::abind(me.numeric$Indiv,me.factor$Indiv, along=3)
  me.total <- cbind(me.numeric$Total,me.factor$Total)
  me.direct <- cbind(me.numeric$Direct,me.factor$Direct)
  me.indirect <- cbind(me.numeric$Indirect,me.factor$Indirect)
 if (delta.method) {
    attr(f, 'se') <- abind::abind(attr(me.numeric$Indiv,'se'),
                                  attr(me.factor$Indiv,'se'), along=3 )
    attr(me.total, 'se') <- cbind(attr(me.numeric$Total,'se'),
                                  attr(me.factor$Total,'se'))
    attr(me.direct, 'se') <- cbind(attr(me.numeric$Direct,'se'),
                                   attr(me.factor$Direct,'se'))
    attr(me.indirect, 'se') <- cbind(attr(me.numeric$Indirect,'se'),
                                     attr(me.factor$Indirect,'se'))


  }
  list(
    #Indiv=f,
     Direct=me.direct, Indirect=me.indirect, Total=me.total)
}

############# predict xb
predict.xb.pmlsprobit<- function(mod, newdata=NULL, inv_rho=T) {
  #sigma <- sqrt(diag(mod$f$Sigma_rho))
  if (is.null(newdata)) {
    X <- model.matrix(mod$terms, mod$model$mf)
    if (mod$slx) {
      Xlag <-   model.matrix(mod$terms.xlag, mod$model$mf.xlag)
      colnames(Xlag) <-  paste0("W2*",colnames(Xlag))
      Xlag <- mod$W2%*%Xlag[,colnames(Xlag) %in% names(mod$beta), drop=F]
      X <- cbind(X,Xlag)
    }
  } else {
   X <-  model.matrix(mod$terms, model.frame(mod$terms, newdata))
   if (mod$slx) {
     Xlag <-  model.matrix(mod$terms.xlag, model.frame(mod$terms.xlag, newdata))
     colnames(Xlag) <-  paste0("W2*",colnames(Xlag))
     Xlag <- mod$W2%*%Xlag[,colnames(Xlag) %in% names(mod$beta), drop=F]
     X <- cbind(X,Xlag)
   }
     }
if (inv_rho==T)
  xb <- (mod$f$Inv_rho%*%X)%*%mod$beta
else
  xb <- X%*%mod$beta
dimnames(xb)[1] <- dimnames(X)[1]
xb
}

predict.prob.pmlsprobit<- function(mod, newdata=NULL) {
  xb <- predict.xb.pmlsprobit(mod, newdata)
  sigma <- sqrt(diag(mod$f$Sigma_rho))
  p <- pnorm(xb/ sigma)
  p
}


####################### FACTOR VARIABLES ########################
predict.me.factor.pmlsprobit.v2 <- function(mod, variables=NULL, delta.method=F) {
   factor.vars <- which(attr(mod$terms,"dataClasses")[-1]=="factor")
  if (!is.null(variables)) {
    stopifnot( all(is.character(variables)) )
    if (!all(variables %in% names(factor.vars)  ))
      stop('some of the elements in variables are not in the covariates of the model')
    factor.vars <- factor.vars[which(names(factor.vars) %in% variables)]
  }
  attrX <- attributes(model.matrix(mod$terms, mod$model$mf))
  nobs <- nrow(mod$model$mf)
  sigma <- sqrt(diag(mod$f$Sigma_rho))
  invA <- mod$f$Inv_rho
  if (mod$model.type!= "SARAR")
    invB <- diag(nobs) else {  invB <- mod$f$Inv_lambda }
  nmargin <- sum(sapply(mod$model$mf[names(factor.vars)], function(x){length(levels(x))}, simplify=T))-length(factor.vars)
  f <- array(dim=c(nobs,nobs,nmargin),
             dimnames = list(dimnames(mod$model[[1]])[[1]],dimnames(mod$model[[1]])[[1]],attrX$dimnames[[length(attrX$dimnames)]][which(attrX$assign %in% factor.vars)]))
  if (delta.method) {
    grad <- list()
    gradd <- array(dim=c(nobs,length(mod$beta),nobs) )
    grad.rho <- grad.lambda <- matrix(nrow=nobs, ncol=nobs)
    dSigma.drho <- invA%*%mod$W%*%invA%*%invB%*%t(invB%*%invA)
    dSigma.drho <- dSigma.drho + t(dSigma.drho)
    if (mod$model.type == "SARAR") {
      dSigma.dlambda <-  invA%*%invB%*%mod$M%*%invB%*%t(invB%*%invA)
      dSigma.dlambda <-  dSigma.dlambda + t(dSigma.dlambda)
    }
  }
  j <- 1

  for (facvar in factor.vars) {
    fv <- names(factor.vars)[factor.vars==facvar]
    # me.total <- me.indirect <- me.direct
    for (l in levels(mod$model$mf[,fv])[-1]  ) {
      for (i in 1:nobs) {
        df1 <- df0 <- mod$model$mf
        df1[i,fv] <- l
        df0[i,fv] <- levels(mod$model$mf[,fv])[1] ##sto mettendo il livello 1 ma in realtà dovrebbe essere il livello scelto dall'utente
        X1 <- model.matrix(mod$terms, df1)
        X0 <- model.matrix(mod$terms, df0)
        #browser()
        if (mod$slx) {
          df1lag <- df0lag <- mod$model$mf.xlag
          df1lag[i,fv] <- l
          df0lag[i,fv] <- levels(mod$model$mf[,fv])[1]
          Xlag1 <- model.matrix(mod$terms.xlag, df1lag)
          Xlag0 <- model.matrix(mod$terms.xlag, df0lag)
          Xlag1 <- mod$W2 %*% Xlag1[,-which(colnames(Xlag1)=="(Intercept)"), drop=F]
          Xlag0 <- mod$W2 %*% Xlag0[,-which(colnames(Xlag0)=="(Intercept)"), drop=F]
          colnames(Xlag1) <- colnames(Xlag0) <- paste0("W2*",colnames(Xlag1))
          X1 <- cbind(X1, Xlag1)
          X0 <- cbind(X0, Xlag0)
        }
        f[i,,j] <- predict.prob.pmlsprobit(mod, df1)-predict.prob.pmlsprobit(mod, df0)
        if (delta.method) {

          z1 <- predict.xb.pmlsprobit(mod,df1)
          z0 <- predict.xb.pmlsprobit(mod,df0)

          gradd[i,,] <-  (sweep((invA%*%X1),1,dnorm(z1)/sigma,"*") -  sweep((invA%*%X0),1,dnorm(z0)/sigma,"*"))

          xb1 <- predict.xb.pmlsprobit(mod, df1,inv_rho=F)
          xb0 <- predict.xb.pmlsprobit(mod, df0,inv_rho=F)

          grad.rho[i,] <- dnorm(z1)/sigma*(invA%*%mod$W%*%invA%*%xb1 - 0.5*sigma^(-2)*diag(dSigma.drho)*invA%*%xb1) -
            dnorm(z0)/sigma*(invA%*%mod$W%*%invA%*%xb0 - 0.5*sigma^(-2)*diag(dSigma.drho)*invA%*%xb0)

          if (mod$model.type == "SARAR") {
            dnorm(z1)*( - 0.5*sigma^(-3)*diag(dSigma.dlambda)*invA%*%xb1) -
              dnorm(z0)*( - 0.5*sigma^(-3)*diag(dSigma.dlambda)*invA%*%xb0)
          }
        } #end if delta.method
      }
      if (delta.method) {

        grad[[j]] <- gradd
        grad[[j]] <- abind::abind(grad[[j]], grad.rho, along=2)
        if (mod$model.type == "SARAR")
          grad[[j]] <- abind::abind(grad[[j]], grad.lambda, along=2)
        dimnames(grad[[j]]) <- list(NULL, names(mod$estimate),NULL)
        }

      j <- j+1
    }

  }
 dimnames(f)[1] <- dimnames(f)[1] <- dimnames(mod$model[[1]])[[1]]
  me.direct <- apply(f,3,diag)
  me.total <- apply(f,3,rowSums)
  me.indirect <- (me.total - me.direct)/(nobs-1)
  me.total <- me.total/nobs
  if (delta.method)  {

    grad.dir <- lapply(grad, function(x) { apply(x, 2,diag)})
    grad.total <- lapply(grad, function(x) { apply(x, 2,rowSums)})
    grad.indirect <-  lapply(grad, function(x) { apply(x, 2,rowSums) -  apply(x, 2,diag)})

    se.dir <- sqrt(sapply(grad.dir, function(x) {apply(x ,1, function(g) {t(g)%*%mod$vcov%*%(g)}  )} ))
    se.indirect <- sqrt(sapply(grad.total, function(x) {apply(x ,1, function(g) {t(g)%*%mod$vcov%*%(g)/((nobs-1)^2)}  )} ))
    se.total <- sqrt(sapply(grad.indirect, function(x) {apply(x ,1, function(g) {t(g)%*%mod$vcov%*%(g)/(nobs^2)}  )} ))


    se.f <-  simplify2array(lapply(grad, function(X) apply(X,  c(1,3), function(x) sqrt(t(x)%*%mod$vcov%*%x ))))

    dimnames(se.f)[[3]]  <- colnames(se.total) <- colnames(se.dir) <- colnames(se.indirect) <-   colnames(me.direct)
    attr(me.direct,"gradient") <- grad.dir
    attr(me.total,"gradient") <- lapply(grad.total ,"/",nobs)
    attr(me.indirect,"gradient") <- lapply(grad.indirect ,"/",nobs-1)
    attr(f, "se") <- se.f
    attr(me.direct,"se") <- se.dir
    attr(me.total,"se") <- se.indirect
    attr(me.indirect,"se") <- se.total
  }
  #
  list(
    #Indiv=f,
    Direct=me.direct, Indirect=me.indirect, Total=me.total)

}



##############################
predict.me.numeric.pmlsprobit.v2 <- function(mod, variables=NULL, delta.method=T) {
  sigma <- sqrt(diag(mod$f$Sigma_rho))
  xb <- predict.xb.pmlsprobit(mod,inv_rho=F)
  z <- predict.xb.pmlsprobit(mod) /sigma
  invA  <- mod$f$Inv_rho
  me <- dnorm(z)/sigma
  nobs <- length(z)

  if (is.null(variables)) {
    X <- model.matrix(mod$terms, mod$model$mf)
    numeric.vars <- names(which(attr(mod$terms,"dataClasses")[-1]=="numeric"))
    all.beta <- mod$beta
    all.names <- colnames(X)[colnames(X) %in% numeric.vars]
    if (mod$slx) {
      numeric.vars.slx <- names(which(attr(mod$terms.xlag,"dataClasses")[-1]=="numeric"))
      Xlag <- model.matrix(mod$terms.xlag, mod$model$mf.xlag)
      names.Xlag <- colnames(Xlag)[colnames(Xlag) %in% numeric.vars.slx]
      all.names <- unique(c(all.names,names.Xlag))
      all.beta <- rep(0, 2*length(all.names))
      names(all.beta) <- c(all.names, paste0("W2*",all.names))
      for (i in names(mod$beta)[names(mod$beta)!="(Intercept)"] )
        all.beta[i] <- mod$beta[i]
      colnames(Xlag) <-  paste0("W2*",colnames(Xlag))
      Xlag <- mod$W2%*%Xlag[,colnames(Xlag) %in% names(mod$beta), drop=F]
      X <- cbind(X,Xlag)
    }
    if (delta.method) {
      dq.dbeta  <- X

    }
  } else {
    if (!any(variables %in% colnames(mod$model$mf)))
      stop('variables must contain at least a function one of the variables used to estimate the model')
    all.names <- variables
    if (delta.method) {
      dq.dbeta <- numDeriv::jacobian(
        function(beta, model=mod){
          model$beta <- beta
          predict.xb.pmlsprobit(model, inv_rho=F)}
        , mod$beta,model=mod)
      #d2q.dxdbeta <- array(0, dim=c(nobs, length(mod$beta),nobs),dimnames = list(NULL, names(mod$beta),NULL))
    }

  }
  f <- array(dim=c(nobs,nobs, length(all.names) ),dimnames = list(dimnames(mod$model[[1]])[[1]],dimnames(mod$model[[1]])[[1]],all.names) )
  mfx <- array(dim=c(nobs,3, length(all.names) ),dimnames = list(dimnames(mod$model[[1]])[[1]],c("Direct","Indirect","Total"),all.names) )
  grad <- list()
  for (i in all.names){  ##looping over all the variables for which mfx are calculated
    if (is.null(variables)) {
      dq.dx <- all.beta[i]*diag(nobs)
      if (mod$slx)
        dq.dx <-  dq.dx + all.beta[paste0('W2*',i)]*mod$W2
    } else {
      dq.dx <- numDeriv::jacobian(
        function(x, index,model=mod){
          newdata <- model$model$mf
          newdata[,index] <- x
          predict.xb.pmlsprobit(model, newdata=newdata, inv_rho=F)},
        mod$model$mf[,i], index=which( colnames(mod$model$mf)==i))

    }
    f[,,i] <- sweep(invA  ,1,me,'*')%*%dq.dx
    mfx[,"Direct",i] <- diag(f[,,i])
    mfx[,"Total",i] <- rowSums(f[,,i])
    mfx[,"Indirect",i] <-  mfx[,"Total",i] -  mfx[,"Direct",i]
    if (delta.method) {
      grad[[i]] <-  sweep(invA%*%dq.dbeta , 1,-z*dnorm(z)/(sigma^2),'*')
      grad[[i]] <- apply(invA%*%dq.dx,2, function(x) {sweep(grad[[i]],1,x,'*') })
      grad[[i]] <- array(grad[[i]], c(nobs,length(mod$beta),nobs))
      if (is.null(variables)) {

        d2q.dxdbeta <- array(0, dim=c(nobs, length(mod$beta),nobs),dimnames = list(NULL, names(mod$beta),NULL))
        for (r in 1:nobs) {
          d2q.dxdbeta[r,i,r] <- 1
          if (mod$slx && (paste0('W2*',i) %in% names(mod$beta)) )
            d2q.dxdbeta[r,paste0('W2*',i),] <- t(mod$W2[r,])
        }
        d2q.dxdbeta <- lapply(seq(nobs), function(x) t(d2q.dxdbeta[x,,])  )
      }
      else {

        xbeta.bidi <- function(param, model=mod, which.var=1, which.obs=1) {
          x <- head(param, nrow(model$model$mf))
          beta <- tail(param,length(model$beta))
          model$beta <- beta
          model$model$mf[,which.var] <- x
          predict.xb.pmlsprobit(model, newdata=model$model$mf, inv_rho=F)[which.obs]
        }
        # if (!parallel){
        d2q.dxdbeta <- lapply(seq(nobs), function(x) {
          numDeriv::hessian(xbeta.bidi, c(mod$model$mf[,i], mod$beta), which.var = which( colnames(mod$model$mf)==i) , which.obs=x)[1:nobs,-(1:nobs)]
        } )
        # } else {
        # clust <-   makeCluster(detectCores()-1)
        # clusterExport(clust, ls(environment()), envir=environment()  )
        # d2q.dxdbeta <- parallel::parLapply(clust,seq(nobs), function(x) {
        #   numDeriv::hessian(function(param, model=mod, which.var=1, which.obs=1) {
        #     x <- head(param, nrow(model$model$mf))
        #     beta <- tail(param,length(model$beta))
        #     model$beta <- beta
        #     model$model$mf[,which.var] <- x
        #     predict.xb.pmlsprobit(model, newdata=model$model$mf, inv_rho=F)[which.obs]
        #   }
        #   , c(mod$model$mf[,i], mod$beta), which.var = which( colnames(mod$model$mf)==i) , which.obs=x)[1:nobs,-(1:nobs)]
        # } )
        # }
      }#
      Ad2q  <- lapply(d2q.dxdbeta, function(x) t(invA%*%x) )
      Ad2q <- aperm(simplify2array(Ad2q),perm=c(3,1,2))
      Ad2q <-  sweep(Ad2q , 1,dnorm(z)/(sigma^2),'*')
      grad[[i]] <- grad[[i]] + Ad2q

      if (mod$model.type!= "SARAR")
        invB <- diag(nobs) else {  invB <- mod$f$Inv_lambda }

      dSigma.drho <- invA%*%mod$W%*%invA%*%invB%*%t(invB%*%invA)
      dSigma.drho <- dSigma.drho + t(dSigma.drho)

      dz.drho <- -0.5*sigma^(-3)*diag(dSigma.drho)*(invA%*%xb) +
        1/sigma * (invA%*%mod$W%*%invA%*%xb)

      df.drho <- -0.5*sigma^(-3)*diag(dSigma.drho)*(invA%*%dq.dx) +
        1/sigma*(-mod$W%*%dq.dx)


      grad.rho <- sweep(f[,,i] ,1 ,-z*dnorm(z)*dz.drho,"*") + sweep(df.drho, 1, dnorm(z),"*")
      grad[[i]] <- abind::abind(grad[[i]], grad.rho, along=2)
      if (mod$model.type== "SARAR") {
        dSigma.dlambda <- invA%*%invB%*%mod$M%*%invB%*%t(invB%*%invA)
        dSigma.dlambda <- dSigma.dlambda + t(dSigma.dlambda)

        dz.dlambda <- -0.5*sigma^(-3)*diag(dSigma.dlambda)*(invA%*%xb)
        df.dlambda  <- -0.5*sigma^(-3)*diag(dSigma.drho)*(invA%*%dq.dx)

        grad.lambda <- sweep(f[,,i] ,1 ,-z*dnorm(z)*dz.dlambda,"*") + sweep(df.dlambda, 1, dnorm(z),"*")
        grad[[i]] <- abind::abind(grad[[i]], grad.lambda, along=2)
      }
    }
  }

  me.direct <-  matrix(mfx[,"Direct",,drop=F], ncol=length(all.names), dimnames = list(dimnames(mod$model[[1]])[[1]], all.names) )
  me.total <- matrix(mfx[,"Total",,drop=F]/nobs, ncol=length(all.names), dimnames = list(dimnames(mod$model[[1]])[[1]], all.names) )
  me.indirect <- matrix(mfx[,"Indirect",,drop=F]/(nobs-1), ncol=length(all.names), dimnames = list(dimnames(mod$model[[1]])[[1]], all.names) )

    if (delta.method)  {
    grad.dir <- lapply(grad, function(x) { (apply(x, 2,diag))})
    grad.total <- lapply(grad, function(x) { (apply(x, 2,rowSums))})
    grad.indirect <-  lapply(grad, function(x) { (apply(x, 2,rowSums) -  apply(x, 2,diag))})
    se.direct <- se.indirect <- se.total <- matrix(NA, nrow=nobs, ncol= length(all.names), dimnames = list(NULL,all.names) )
    for (v in all.names) {

      se.direct[,v] <- sqrt(apply(grad.dir[[v]], 1, function(x) { t(x)%*%mod$vcov%*%x }   ))
      se.indirect[,v] <- sqrt(apply(grad.indirect[[v]], 1, function(x) { t(x)%*%mod$vcov%*%x }   ))
      se.total[,v] <- sqrt(apply(grad.total[[v]], 1, function(x) { t(x)%*%mod$vcov%*%x }   ))
    }

    se.f <-  simplify2array( lapply(grad, function(X) apply(X,  c(1,3), function(x) sqrt(t(x)%*%mod$vcov%*%x ))))
    attr(me.direct,"gradient") <- grad.dir
    attr(me.total,"gradient") <- lapply(grad.total ,"/",nobs)
    attr(me.indirect,"gradient") <- lapply(grad.indirect ,"/",nobs-1)
    attr(f, "se") <- se.f
    attr(me.direct,"se") <- se.direct
    attr(me.total,"se") <- se.indirect
    attr(me.indirect,"se") <- se.total
  }
  ret <- list(
    #Indiv=f ,
    Direct=me.direct, Indirect=me.indirect, Total=me.total)
  ret
}





#' Average partial (i.e., marginal) total, direct and indirect effects from a pmlsprobit object
#'
#' @param mod 	Object of class inheriting from "pmlsbp"
#' @param variables a string including the variable names of the covariates for
#' which the marginal effects should be calculated using numeric derivatives.  It is valid only if \code{type="me"}.
#' This option should be used only if the model contains non-linear terms in the specification \eqn{X \beta}
#'
#' @return
#' @export
#'
#' @examples
ape <- function(mod,variables=NULL){
  pe <- predict(mod, 'me',variables, delta.method = T)
  pe$Indiv <- NULL
  ape <- lapply(pe, colMeans)

  grad <- attr(pe$Direct,'gradient')
  simplify2array(lapply(grad, colMeans))
  g <- lapply(pe,
              function(x) {
                simplify2array(lapply(attr(x,'gradient'), colMeans))
              } )
  V <- lapply(g, function(x){ t(x)%*%mod$vcov%*%x } )
  results <- list()
  for (j in names(g))
    results[[j]] <- cbind("Estimate"=ape[[j]], "Std. error"=sqrt(diag(V[[j]])) )

  results <- lapply(results, function(x) {
    z <- x[,1]/x[,2]
    p <- 2*pnorm( -abs( z))
    cbind(x,"z value" = z,"Pr(> z)" = p)
  } )
  class(results) <- "ape.pmlsprobit"
  results
}


#' @export
print.ape.pmlsprobit<-function(ape, digits = max(3, getOption("digits") - 3)) {
  for (i in 1:length(ape)) {
    cat("Average", paste( names(ape)[[i]], "Partial Effects\n"))
    stats::printCoefmat(ape[[i]], digits=digits, signif.legend = (i==length(ape)))
    cat("\n\n")
    }
}


#' @describeIn Predicted values based on pmlsprobit object
#' @export
print.pe.pmlsprobit<-function(pe, digits = max(3, getOption("digits") - 3)) {
lapply(pe, print.default, digits=digits)
}
