#' pmlsbp: Fitting Spatial Probit Models
#' @description Estimation of Partial Maximum Likelihood Spatial Autoregressive (SAR) Probit Models and Spatial Autoregressive Probit Models with Autoregressive Disturbances (SARAR).
#' The model structure for the SARAR model is: \cr
#' \tabular{c}{
#' \eqn{y^* = \rho W y + X \beta + W_2 Z \gamma} + u  \cr
#' \eqn{u = \lambda M u + \epsilon} \cr
#' }
#' where \eqn{\epsilon \sim N(0,\sigma^2)} and \eqn{y=1(y^*>0)}. The SAR model assumes \eqn{\lambda=0}
#' @usage pmlsbp(formula,data, model="SAR", grouping=2, W=NULL,
#'  zero.policy =spatialreg::get.ZeroPolicyOption(), M=NULL, formula_xlag=NULL,
#'  W2=NULL, y=NULL, X=NULL, method_inv="chol", start=NULL, subset=NULL,
#'  na.action=na.fail,qu=Inf,
#'  mvtnorm_control=list(M=25e3, sim_type="mc" , tol = .Machine$double.eps, fast = FALSE),
#'  finalHessian=TRUE, method="bfgsr",print.level=2, vce.type=c("asy") ,
#'  Conley=list(coords=NULL,LM=2),nBoot=1e3 , ...)
#' @param formula an object of class "formula": a symbolic description of the model to be fit. The details of model specification are given for lm()
#' @param data an optional data frame containing the variables in the model. By default the variables are taken from the environment which the function is called.
#' @param model a character (\code{"SAR"} or \code{"SARAR"}) defining the model to be estimated. The default is \code{"SAR"}.
#' @param grouping an integer defining the number of observations to include in the tuples
#' to be used for estimation, the default is 2.
#' @param W a listw object created for example by \code{nb2listw}
#' @param zero.policy default NULL, use global option value; if TRUE assign zero to the lagged
#' value of zones without neighbours, if FALSE (default) assign NA - causing
#' \code{pmlsbp} to terminate with an error
#' (see also \code{\link[spatialreg]{set.ZeroPolicyOption}}) and the same option
#' in \code{\link[spatialreg]{GMerrorsar}})
#' @param M a listw object created for example by \code{nb2listw} to be used
#' for the model with spatial autoregressive disturbances. This is relevant only if \code{model="SARAR"}.
#' @param formula_xlag an object of class "formula": a symbolic description of the variables to be spatially lagged
#' @param W2 a listw object created for example by \code{nb2listw} to be used
#' to create the spatial weight matrix to multiply the covariates. This is relevant only if \code{formula_xlag} is not \code{NULL}.
#' @param method_inv a character for setting the inversion of \eqn{I-\rho*W}, either \code{"solve"}, \code{"chol"} or \code{"fast"}.
#' \code{"solve"} uses \code{solve}, "chol" performs the inversion of the matrix from its Choleski decomposition (see \code{\link[Matrix]{chol2inv}}), \code{"fast"} uses an approximation that requires option \code{qu}
#' @param start numeric vector, initial value of parameters.
#' @param subset an optional vector specifying a subset of observations to be used for fitting
#' @param na.action a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude}
#' with consequences for residuals and fitted values - in these cases the weights
#' list will be subsetted to remove NAs in the data. It may be necessary to set
#' \code{zero.policy} to TRUE because this subsetting may create no-neighbour observations.
#' Note that only weights lists created without using the \code{glist} argument to nb2listw may be subsetted.
#' @param qu an integer to be used only if \code{method_inv="fast"}
#' @param mvtnorm_control an optional object of class \code{"list"} that includes the controls for
#' the calculation of the multivariate normal probabilities.  The list members should include
#' \code{M}, \code{rmtol} and \code{fast} (see \code{\link[mvtnorm]{lmpvnorm}}) and
#' \code{simtype}, either \code{"mc"} for Monte-Carlo procedure or \code{"qmc"}, quasi-Monte-Carlo procedure)
#' The default is \code{list(M=25e3, sim_type="mc" , tol = .Machine$double.eps, fast = FALSE)}
#' @param finalHessian how (and if) to calculate the final Hessian. Either \code{FALSE} (not calculate), \code{TRUE}
#' (use analytic/numeric Hessian) or \code{"bhhh"/"BHHH"} for information equality approach
#'  (see \code{\link[maxLik]{maxLik}}).
#' @param method maximisation method to be used by \code{maxLik}. The default is \code{method="bfgsr"}
#' Broyden-Fletcher-Goldfarb-Shanno as implemented in R. See also the same option in (see \code{\link[maxLik]{maxLik}}).
#' @param print.level Integer, debugging information. Larger number prints more details (see \code{\link[maxLik]{maxLik}}).
#' @param vce.type  an optional character specifying the type of variance covariance
#' estimation  (VCE) to be performed. It can be either \code{"asy"} for the asymptotic VCE,
#' \code{"bootstrap"} for bootstrap  VCE or \code{"mConley"} for the modified Conley approach.
#' The default is "asy".
#' @param Conley a list with two members, which is relevant only for \code{vce.type="mConley"}.
#' The first member (\code{coords}) is a matrix with the coordinates of each observation;
#' the second member (\code{LM}) is the number of neighbors to be used.
#' The default is \code{list(coords=NULL,LM=2)}
#' @param nBoot an integer, the number of bootstrap samples to be used for vce.type="bootstrap".
#' It is relevant only if \code{vce.type="bootstrap"}.
#' @param ...
#'
#' @return \code{pmlsbp} returns an object of class \code{"pmlsprobit"}.
#' The function \code{summary} is used to obtain and print a summary of the results.
#' A object of class \code{"pmlsprobit"} is list containing the following components:
#' \tabular{ll}{
#' \code{beta} \tab a named vector including the \eqn{\beta} and \eqn{\gamma}, the coefficents attached to the variables resulting from \code{formula} and  \code{formula_xlag} \cr
#' \code{call} \tab the call \cr
#' \code{code} \tab \cr
#' \code{estimate} \tab a named vector of model estimates \cr
#' \code{f} \tab a list including useful matrices to be used in \code{predict} \cr
#' \code{maximum} \tab the value of the partial log-likelihood given the data and the model estimates \cr
#' \code{message} \tab convergence message \cr
#' \code{model} \tab \cr
#' \code{model.type} \tab either \code{"SAR"} or  \code{"SARAR"}  \cr
#' \code{rho} \tab a named scalar for the \eqn{\rho} parameter \cr
#' \code{slx} \tab \code{TRUE} if the model includes spatially lagged covariates \cr
#' \code{start} \tab a named vector, the starting values used in estimation \cr
#' \code{terms} \tab 	the terms object used. \cr
#' \code{vcov} \tab the variance-covariance matrix of the estimates \cr
#' \code{W} \tab matrix \eqn{W}  \cr
#' \code{lambda} \tab a named scalar for parameter \eqn{lambda} (only for SARAR models)  \cr
#' \code{M} \tab matrix \eqn{M} (only for SARAR models)  \cr
#' \code{W2} \tab matrix \eqn{W_2} (only if \code{slx=TRUE})  \cr
#' }

#' @export
#' @references Bille', A. G., & Leorato, S. (2020). Partial ML estimation for spatial autoregressive nonlinear probit models with autoregressive disturbances. Econometric Reviews, 39(5), 437-475.
#' @examples
pmlsbp <- function(formula,data, model="SAR", grouping=2, W=NULL,zero.policy =spatialreg::get.ZeroPolicyOption(),
                         M=NULL, formula_xlag=NULL, W2=NULL,  method_inv="solve",
                         start=NULL, subset=NULL, na.action=na.fail,qu=Inf, iterlim=1000,
                         mvtnorm_control=list(M=25e3, sim_type="mc" , tol = .Machine$double.eps, fast = FALSE),
                         finalHessian=TRUE, method="bhhh",print.level=2, vce.type='asy' , Conley=list(coords=NULL,LM=2),nBoot=1e3 , spectral=F, tol.solve= .Machine$double.eps, version=0, ...) {
  match.arg(model, c("SAR","SARAR"), several.ok = FALSE)
  match.arg(method_inv, c("chol","solve","fast"), several.ok = FALSE)
  stopifnot(is.numeric(qu))
  stopifnot(inherits(formula,"formula"))
  match.arg(vce.type, c("asy", "mConley", "bootstrap"), several.ok = F)



  if (is.null(mvtnorm_control$M))
    mvtnorm_control$M <- 25e3
  if (is.null(mvtnorm_control$tol))
    mvtnorm_control$tol <- .Machine$double.eps
  if (is.null(mvtnorm_control$fast))
    mvtnorm_control$fast <- FALSE
  if (is.null(mvtnorm_control$sim_type))
    mvtnorm_control$sim_type <- "mc"
  stopifnot(mvtnorm_control$sim_type %in% c("mc","qmc"))

  #serve per chiamare model.frame con relative opzioni, senza generare errori
  call <- match.call()
  mmf <- match.call(expand.dots = FALSE)
  mmf$na.action <- na.action
  mmf$subset <- subset
  m <- match(c("formula", "data","na.action","subset"), names(mmf), 0L)
  mf <- mmf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)

  gav <- mf
  gav.formula <- all.vars(formula)



  if (!is.null(formula_xlag)) {
    mf.xlag <- mf
    mf.xlag$formula <- formula_xlag
    gav.formula <- c(gav.formula,all.vars(formula_xlag))
    mf.xlag <- eval(mf.xlag, parent.frame())
    mt.xlag <- attr(mf.xlag, "terms")

  }

  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  gav$formula <- reformulate(gav.formula)
  gav <- eval(gav,parent.frame())

  if (is.empty.model(mt))
    stop("'formula' is an empty model")
  y <- stats::model.response(mf, "any")
  if (is.character(y) | is.factor(y))
    stop ('y must be a numeric variable')
  stopifnot(length(unique(y))==2)
  y <- as.integer(as.logical(y))
  if (!all(y %in% c(0,1))) {
    stop('spatial_probit failed to convert the response variable into (0,1) values,
         please check the response variable')
  }
  X <- stats::model.matrix(mt, mf)
  if (nrow(X)!=length(y))
    stop("independent variables and response variable have different dimensions")
  if (is.null(subset))
    subset <- rep(T, nrow(X))
  if (is.null(W))
    stop("W is mandatory")
  if (!inherits(W,"matrix") & !inherits(W,"listw"))
    stop("W must be a matrix or a listw")
  if (inherits(W, "listw")) {
    if (!any(W$style==c("W","B","C")))
      stop("the style of W must be W or B")
    styleW<-W$style
    W <- subset(W, !(1:length(W$neighbours) %in% c(attr(mf, "na.action"), which(!subset))),zero.policy = zero.policy)
    W <- spdep::listw2mat (W)
  } else {
    W <- subset(W, !(1:dim(W)[2] %in% c(attr(mf, "na.action"), which(!subset))))
    }
  eig <-  list(W=range(eigen(W, symmetric = F)$values))
  if (styleW=='B' & spectral) {
    W <- W/max(abs(eig$W))
    eig <-  list(W=range(eigen(W, symmetric = F)$values))
  }
  if (nrow(W)!=length(y))
    stop("spatial matrix W and response variable have different dimensions")

  if (model=="SARAR" )   {
    if ( is.null(M) )
      stop("M is mandatory for SARAR model")
    if (!inherits(M,"matrix") & !inherits(M,"listw"))
      stop("M must be a matrix or a listw")
    if (inherits(M, "listw")) {
      if (!any(M$style==c("W","B","C")))
        stop("the style of M must be W or B")
      M <- subset(M, !(1:length(M$neighbours) %in% c(attr(mf, "na.action"), which(!subset))),zero.policy = zero.policy)
      M <- spdep::listw2mat(M)
    } else {
      M <- subset(M, !(1:dim(M)[2] %in% c(attr(mf, "na.action"), which(!subset))))
      }
    eig$M <- range(eigen(M, symmetric = F)$values)

    if (nrow(M)!=length(y))
      stop("spatial matrix M and response variable have different dimensions")
  }

  if (!is.null(formula_xlag)) {
    if ( is.null(W2) )
      W2 <- W
    else {
      if (!inherits(W2,"matrix") & !inherits(W2,"listw"))
      stop("W2 must be a matrix or a listw")
    if (inherits(W2, "listw")) {
       if (!any(W2$style==c("W","B","C")))
        stop("the style of W2 must be W or B")
      styleW2 <-  W2$style
      W2 <- subset(W2, !(1:length(W2$neighbours) %in% c(attr(mf, "na.action"), which(!subset))))
      W2  <- spdep::listw2mat(W2)
    } else   {
      W2 <- subset(W2, !(1:dim(W2)[2] %in% c(attr(mf, "na.action"), which(!subset))))
      if (all(rowSums(W2)==1))
        styleW2 <-  "W"
        }
    }


    Xlag <- stats::model.matrix(mt.xlag, mf.xlag)
    attr.Xlag <- attributes(Xlag)
    attr.names <- names(attributes(Xlag))
    attr.names <- attr.names[attr.names != 'names']


    Xlag <- W2%*%stats::model.matrix(mt.xlag, mf.xlag)
    #cn <-  colnames(Xlag)
    if (styleW2=="W" && any(colnames(Xlag)=="(Intercept)")) {
       #cn <- cn[-which(colnames(Xlag)=="(Intercept)")]
       Xlag <- as.matrix(Xlag[,-which(colnames(Xlag)=="(Intercept)")])
       attr.Xlag$dimnames[[2]] <- (attr.Xlag$dimnames[[2]])[-which(attr.Xlag$dimnames[[2]]=="(Intercept)")]
       attr.Xlag$dim[2] <- attr.Xlag$dim[2]-1
       attr.Xlag$assign <- attr.Xlag$assign[-which(attr.Xlag$assign==0) ]
       #attr(attributes(mod$model$mf.xlag)[[2]], 'intercept')
    }

      attributes(Xlag)[attr.names] <- attr.Xlag[attr.names]
      attr(Xlag, 'dimnames')[[2]] <-   paste0("W2*",attr(Xlag, 'dimnames')[[2]])
    if (nrow(Xlag)!=length(y))
      stop("spatial lags and response variable have different dimensions")
    Xall <- cbind(X,Xlag)
    if (nrow(W2)!=length(y))
      stop("spatial matrix W2 and response variable have different dimensions")
  }
  else Xall <- X
  if (is.null(start)) {
    start <- suppressWarnings(glm(y~Xall-1, family = binomial(link = "probit"))$coefficients)
    names(start) <- colnames(Xall)

    start <- c(start,to_tilde(0, eig$W)[1])
    names(start)[length(start)] <- "rho"
    if (model=="SARAR") {
      start <- c(start,to_tilde(0, eig$M)[1])
      names(start)[length(start)] <- "lambda"
    }
  }
  if ("mConley" %in% vce.type) {
    if (is.null(Conley$coords))
      stop("Conley score requires coordinates data")
    stopifnot(is.data.frame(Conley$coords) | is.matrix(Conley$coords) )
    stopifnot(Conley$LM>1)
    Conley$coords <- subset(Conley$coords, !(1:nrow(Conley$coords) %in% c(attr(mf, "na.action"), which(!subset))))
  }


 groups <- group_matrices(y,Xall,grouping, M=mvtnorm_control$M, sim_type= mvtnorm_control$sim_type)

   if (method!="bobyqa") {
     m.maxlik <- match(c("method", "iterlim","tol","reltol","gradtol", "steptol","lambdatol","qac",
                         "qrtol", "marquardt_lambda0", "marquardt_lambdaStep" , "marquardt_maxLambda"  , "print.level" ,
                         "nm_alpha", "nm_beta", "nm_gamma", "sann_cand", "sann_temp","sann_tmax","sann_randomSeed","finalHessian","constraints" ), names(call), 0L)
     optimizer.call <- call[c(1L,m.maxlik)]
     optimizer.call[[1L]] <- quote(maxLik::maxLik)
     if (version==0 & grouping==2 & nrow(W)%%2==0 )
       optimizer.call$logLik <- ifelse(model=="SAR" ,quote(logLIK_SAR3), quote(logLIK_SARAR))
     else
       optimizer.call$logLik <- ifelse(model=="SAR" ,quote(logLIK_SAR), quote(logLIK_SARAR))
     optimizer.call$start <- quote(start)
     #optimizer.call$tol.solve <-quote(tol.solve)
   } else {
     optimizer.call <- call[c(1L)]
     optimizer.call[[1L]] <- quote(minqa::bobyqa)
     optimizer.call$par <- quote(start)
     #optimizer.call$feval <-quote(iterlim)
     optimizer.call$fn <- ifelse(model=="SAR" ,quote(logLIK_SAR), quote(logLIK_SARAR))
     optimizer.call$bobyqa <- quote(TRUE)
     optimizer.call$control <- quote(list(iprint=print.level , maxfun=iterlim) )
   }

   optimizer.call$X <-  quote(Xall)
   optimizer.call$y <-  quote(y)
   optimizer.call$eig <-  quote(eig)
   optimizer.call$W <-  quote(W)
   optimizer.call$qu <-  quote(qu)
   optimizer.call$method_inv <- quote(method_inv)
   optimizer.call$groups <- quote(groups)
   optimizer.call$mvtnorm_control <-  quote(mvtnorm_control)
   if (model!="SAR")
     optimizer.call$M <- quote(M)
   return <- eval(optimizer.call)
   # browser()

   return$rho <- to_natural(tail(return$estimate,1), eig$W )
  ## return(return)
   if (method=="bobyqa") {
     return$code<- return$ierr
     return$estimate <- return$par
   }
   # if (return$code!=0) {
   #   cat("----------------------------------------------\n WARNING \n")
   #   cat(paste0("converge not achieved: ",return$message,"\n"))
   #   cat("----------------------------------------------\n")
   #   return(NULL)
   # }

  return$groups <- groups
  return$start <- start

  obFun.call <- optimizer.call
  m.obFun <- match(c("X","y","eig","W","qu","method_inv","groups","mvtnorm_control","M"  ), names(optimizer.call), 0L)
  obFun.call <- optimizer.call[c(1L,m.obFun)]
  obFun.call[[1L]] <- ifelse(model=="SAR" ,quote(logLIK_SAR), quote(logLIK_SARAR))
  obFun.call$theta <- quote(return$estimate)


  obFun <- eval(obFun.call)

  return$beta <- return$estimate[1:ncol(Xall)]
  if (model=="SAR")  {
    return$rho <- to_natural(tail(return$estimate,1), eig$W )
    return$rho_tilde <- return$estimate[length(return$estimate)]
    return$estimate[length(return$estimate)] <- return$rho
    G <- diag(c(rep(1,length(return$estimate)-1), attr(return$rho,"jacobian")))
  }
  else {
    return$rho <- to_natural(return$estimate[length(return$estimate)-1], eig$W )
    return$lambda <- to_natural(tail(return$estimate,1), eig$M )
    return$rho_tilde <- return$estimate[length(return$estimate)-1]
    return$lambda_tilde <- return$estimate[length(return$estimate)]
    return$estimate[length(return$estimate)-1] <- return$rho
    return$estimate[length(return$estimate)] <- return$lambda
    G <- diag(c(rep(1,length(return$estimate)-2), attr(return$rho,"jacobian"), attr(return$lambda,"jacobian")))
  }

  ##maxlik uses rho_tilde for the maximization problem, to obtain the standard error the gradient must be adjusted to work with rho using the cain rule
if (method=='bobyqa') {
  return$gradientObs<-attr(obFun, 'gradient')
  return$hessian <- quote(-numericHessian(fn, t0=return$estimate, X=Xall, y = y, eig = eig, W = W, qu = qu, method_inv = method_inv,
                                    groups = groups, mvtnorm_control = mvtnorm_control, bobyqa=T))
  return$hessian[[2]][[2]] <- obFun.call[[1]]
  return$hessian<- eval(return$hessian)
}
  invH <- solve(as(return$hessian,"Matrix")) ##### devo tirare fuori hessiano
  V <-  list()
  if ("asy" %in% vce.type) {
    J <- crossprod(return$gradientObs)
    V <- G%*%invH%*%J%*%invH%*%t(G)
  }

  if ("mConley" %in% vce.type) {
      mean_coord <- aggregate(Conley$coords, list(return$groups$y), FUN=mean)
      knn <- spdep::knearneigh(mean_coord, k=Conley$LM-1)$nn     #matrix of knn for each site (pair of points)
      knn <-  cbind(matrix(1:nrow(knn), ncol=1), knn)
      KM <- 2*(1-abs(0:(Conley$LM-1)/Conley$LM))
      J <- matrix(0,ncol(return$gradientObs),ncol(return$gradientObs))
      for (j in 2:ncol(knn)) {
        J <- J+KM[j]*t(return$gradientObs)%*%return$gradientObs[knn[,j],]
      }
      J <-  (J+crossprod(return$gradientObs))
      V <- G%*%invH%*%J%*%invH%*%t(G)
    }     ## M, LM e tau?
  if ("bootstrap" %in% vce.type) { #bootstrap qui
    eta <- t(mvtnorm::rmvnorm(nBoot,sigma=attr(obFun, "f")$Sigma_rho))
    ystar <- t(tcrossprod(rep(1, ncol(eta)), xb) ) + eta
    ystar <- ifelse(ystar>0,1,0)
    bootstrap.call <- maxlik.call
    bootstrap.call$start <- quote(return$estimate)
    bootstrap.call$print.level <- quote(0)
    boot_estimate <- matrix(NA, nrow=nBoot, ncol=length(return$estimate))
    cat("bootstrap iterations\n")
    pb <- txtProgressBar(min = 0, max = nBoot, initial = 0,style = 3)
           for (i in 1:nBoot) {
            setTxtProgressBar(pb,i)
            bootstrap.call$y <- quote(ystar[,i])
            bi <- try(eval(bootstrap.call), silent = T)
            if (!( "try-error" %in% class(bi)))
             boot_estimate[i,] <- bi$estimate
           }
    boot_estimate <- na.omit(boot_estimate)
    V <- var(boot_estimate)
    #colnames(V$bootstrap) <- rownames(V$bootstrap) <- names(return$estimate)
    close(pb)
  }
  colnames(V) <- rownames(V) <- names(return$estimate)
  return$vcov <- V
  attr(return$vcov,"type") <- vce.type
  return$eig<-eig
  return$f <- attr(obFun, "f")
  return$model.type <- model
  return$slx <- !is.null(formula_xlag)
  return$W <- W
  return$terms <- mt
  return$model <- list(mf=mf)
  if (model=="SARAR")
    return$M <- M
  if (!is.null(formula_xlag)) {
    return$terms.xlag <- mt.xlag
    return$W2 <- W2
    return$model$mf.xlag <- mf.xlag
  }
  return$call <- call

  keep.results <- c('beta','estimate','call','rho','code','f','maximum','message','model','model.type','slx','start','terms','vcov','W', 'eig')
  if (!is.null(formula_xlag))  keep.results <- c(keep.results, 'W2','terms.xlag')
  if (model=="SARAR") keep.results <- c(keep.results,'M','lambda')
  ret<-return[keep.results]
  class(ret) <- "pmlsprobit"
  ret
  }

to_natural <- function(tilde, eig) {
  ret <-  eig[1]^(-1)+ ((eig[2]^(-1)-eig[1]^(-1))/(1+exp(-tilde)))
  attr(ret,"jacobian") <- (1/(ret-eig[1]^-1)-1/(ret-eig[2]^-1))^-1 ##drho_tilde/drho
  ret
}
to_tilde <- function(x,eig) {
  ret <-  -log((eig[2]^(-1)-eig[1]^(-1))/(x-eig[1]^(-1))-1)
  attr(ret,"jacobian") <-  1/(x-eig[1]^-1)-1/(x-eig[2]^-1) ##drho/drho-tilde (dret/dx)
  ret

}

#' @describeIn pmlsbp computes and returns a list of summary statistics of the fitted spatial probit model
#' @export
summary.pmlsprobit <- function(object, eigentol=1e-12,... ) {
  ## object      object of class "pmlsprobit"
  ##
  ## RESULTS:
  ## list of class "summary.maxLik" with following components:
  ## maximum    : function value at optimum
  ## estimate   : estimated parameter values at optimum
  ## gradient   :           gradient at optimum
  ## code       : code of convergence
  ## message    : message, description of the code
  ## iterations : number of iterations
  ## type       : type of optimisation
  ##
  if(!inherits(object, "pmlsprobit"))
    stop("'summary.pmlsprobit' called on a non-'pmlsprobit' object")
  vcov <- object$vcov

  se <- sqrt(diag(vcov))
  z <- object$estimate/se

  results <- cbind("Estimate"=object$estimate, "Std. error"=se, "z value"=z, "Pr(> z)" = 2*pnorm( -abs( z)) )
  summary <- list(  estimate=results,
                   rho=object$rho,
                   lambda=object$lambda,
                   model.type=object$model.type,
                   loglik=object$maximum,
                   iteration=object$iteration,
                   returnCode=object$code,
                   returnMessage=object$message ,
                   vcov=vcov)
  class(summary) <- "summary.pmlsprobit"
  summary
}

#' @export
print.summary.pmlsprobit <- function( x,
                                          digits = max( 3L, getOption("digits") - 3L ), ... ) {
   cat("Partial Maximum Likelihood estimation\n")
  cat(x$model," model\n")
  cat("Partial Log-Likelihood=", x$loglik, "\n")
  cat(x$returnMessage, "\n")
  cat("---\n")
  printCoefmat( x$estimate, digits = digits )
}

#' @describeIn pmlsbp returns the variance covariance matrix of a \code{"pmlsprobit"} object
#' @export
vcov.pmlsprobit <- function(x) {
  V<-x$vcov
  colnames(V) <- rownames(V) <- names(x$estimate)
  V
}

#' @describeIn pmlsbp returns the coefficients estimate from a \code{"pmlsprobit"} object
#' @export
coef.pmlsprobit <- function(x) {
  x$estimate
}

#' @describeIn pmlsbp prints the call and the coefficient estimates of a \code{"pmlsprobit"} object
#' @export
print.pmlsprobit<-function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n",
                         collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    cat(paste0("Coefficients (",x$model.type," model) :\n"))
    print.default(format(coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}


