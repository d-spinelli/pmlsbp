
logLIK_SAR3 <- function(theta,y,W,X,eig,qu=Inf,method_inv="solve", groups,
                        mvtnorm_control=list(M=25e3, tol = .Machine$double.eps, fast = FALSE), bobyqa=F) {
  method_inv="InvW"
  Method_inv=method_inv

  theta_tilde=theta
  quas=qu
  mvtnorm_control=NULL
  bobyqa=T
  groups=2
  beta      <- as.matrix(theta_tilde[1:(length(theta_tilde)-1)])
  rho_tilde <- theta_tilde[length(theta_tilde)]
  n=length(y)
  n2<-1:n
  n_dispari<-n2[n2%%2!=0]

  func_rho <- f_rho_sar3(rho_tilde,W,X,eig,quas,Method_inv)

  logl <- array(0,dim=c(1,n))
  for (i in n_dispari) {

    j = i+1

    X_rho_i_beta <- as.numeric(t(as.matrix(func_rho$X_rho[i,]))%*%beta)
    X_rho_j_beta <- as.numeric(t(as.matrix(func_rho$X_rho[j,]))%*%beta)
    sigma_ij_ii  <- as.numeric(func_rho$Sigma_rho[i,j]/func_rho$Sigma_rho[i,i])
    sigma_sqrt2  <- as.numeric(func_rho$Sigma_rho[j,j]-func_rho$Sigma_rho[i,j]^2/func_rho$Sigma_rho[i,i])
    sigma_ii     <- as.numeric(sqrt(func_rho$Sigma_rho[i,i]))
    sigma_rho_ij <- as.matrix(func_rho$Sigma_rho[i:j,i:j])

    psi2 <- function (u){
      (X_rho_j_beta + sigma_ij_ii*u)/sqrt(sigma_sqrt2)
    }

    k0=2*y[i]+y[j]+1

    p<-prob.v0(k0,X_rho_i_beta=X_rho_i_beta,X_rho_j_beta=X_rho_j_beta,sigma_ii=sigma_ii,sigma_rho_ij=sigma_rho_ij)


    logl[i] <- log(p)
  }

  # logl<- logl
  logl<-logl[logl!=0]
  logl<-t(logl)
  attr(logl,"gradient")<-negscore_sar(theta_tilde,eig=eig,func_rho=func_rho,y=y,quas=quas,Method_inv=Method_inv)
  return(logl)
}



negscore_sar <- function (theta_tilde,eig,func_rho,y,quas,Method_inv) {

  n=length(y)
  beta<-as.matrix(theta_tilde[1:(length(theta_tilde)-1)])
  rho_tilde<-theta_tilde[length(theta_tilde)]

  jacob<-((eig$W[2]^(-1)-eig$W[1]^(-1))*exp(-rho_tilde))/((1+exp(-rho_tilde))^2)
  n2<-1:n
  n_dispari<-n2[n2%%2!=0]

  score <- array(0,dim=c(length(n_dispari),length(theta_tilde)))
  for (i in n_dispari) {

    j=i+1

    X_rho_i_beta   <- as.numeric(t(as.matrix(func_rho$X_rho[i,]))%*%beta)
    X_rho_j_beta   <- as.numeric(t(as.matrix(func_rho$X_rho[j,]))%*%beta)
    sigma_ii       <- as.numeric(sqrt(func_rho$Sigma_rho[i,i]))
    sigma_jj       <- as.numeric(sqrt(func_rho$Sigma_rho[j,j]))
    sigma_rho_ij   <- as.matrix(func_rho$Sigma_rho[i:j,i:j])
    dotX_rho_i     <- as.numeric(t(as.matrix(func_rho$dotX_rho[i,]))%*%beta/sigma_ii - func_rho$dotSigma_rho[i,i]/(2*func_rho$Sigma_rho[i,i])*(X_rho_i_beta/sigma_ii))
    dotX_rho_j     <- as.numeric(t(as.matrix(func_rho$dotX_rho[j,]))%*%beta/sigma_jj - func_rho$dotSigma_rho[j,j]/(2*func_rho$Sigma_rho[j,j])*(X_rho_j_beta/sigma_jj))
    sigma_ij_ii    <- as.numeric(func_rho$Sigma_rho[i,j]/func_rho$Sigma_rho[i,i])
    sigma_ij_jj    <- as.numeric(func_rho$Sigma_rho[i,j]/func_rho$Sigma_rho[j,j])
    dotsigma_ij_ii <- as.numeric(2*func_rho$dotSigma_rho[i,j]-func_rho$dotSigma_rho[i,i]*func_rho$Sigma_rho[i,j]/func_rho$Sigma_rho[i,i]+
                                   -func_rho$dotSigma_rho[j,j]*func_rho$Sigma_rho[i,j]/func_rho$Sigma_rho[j,j])
    sigma_sqrt1    <- as.numeric(func_rho$Sigma_rho[i,i]-func_rho$Sigma_rho[i,j]^2/func_rho$Sigma_rho[j,j])
    sigma_sqrt2    <- as.numeric(func_rho$Sigma_rho[j,j]-func_rho$Sigma_rho[i,j]^2/func_rho$Sigma_rho[i,i])

    psi1 <- function (u){
      (X_rho_i_beta + sigma_ij_jj*u)/sqrt(sigma_sqrt1)
    }

    psi2 <- function (u){
      (X_rho_j_beta + sigma_ij_ii*u)/sqrt(sigma_sqrt2)
    }

    biv_fun <- function (X_rho_i_beta,X_rho_j_beta) {
      if (det(func_rho$Sigma_rho[i:j,i:j]) > exp(-20))
      {res=1/(det(func_rho$Sigma_rho[i:j,i:j])^(1/2)*2*pi)*exp(-(1/2)*t(c(X_rho_i_beta,X_rho_j_beta))%*%solve(func_rho$Sigma_rho[i:j,i:j])%*%as.matrix(c(X_rho_i_beta,X_rho_j_beta)))      }
      else
      {Sigmaij=func_rho$Sigma_rho[i:j,i:j]
      Sigmaij[1,1]=(1.1)*Sigmaij[1,1]+0.0000000001
      Sigmaij[2,2]=(1.1)*Sigmaij[2,2]+0.0000000001
      res=1/(det(Sigmaij)^(1/2)*2*pi)*exp(-(1/2)*t(c(X_rho_i_beta,X_rho_j_beta))%*%solve(Sigmaij)%*%as.matrix(c(X_rho_i_beta,X_rho_j_beta)))
      }
      return(res)}

    sy=c(2*(y[i]-.5),2*(y[j]-.5))
    k0=2*y[i]+y[j]+1
    p<-prob.v0(k0,X_rho_i_beta=X_rho_i_beta,X_rho_j_beta=X_rho_j_beta,sigma_ii=sigma_ii,sigma_rho_ij=sigma_rho_ij)


    deriv_beta <- (sy[1]/sigma_ii)*dnorm(X_rho_i_beta/sigma_ii)*pnorm(sy[2]*psi2(-X_rho_i_beta))*as.matrix(func_rho$X_rho[i,]) +
      +(sy[2]/sigma_jj)*dnorm(X_rho_j_beta/sigma_jj)*pnorm(sy[1]*psi1(-X_rho_j_beta))*as.matrix(func_rho$X_rho[j,])

    deriv_rho  <- sy[1]*sy[2]*(biv_fun(X_rho_i_beta=X_rho_i_beta,X_rho_j_beta=X_rho_j_beta)/2)*(dotsigma_ij_ii) +
      + sy[1]*dnorm(X_rho_i_beta/sigma_ii)*pnorm(sy[2]*psi2(-X_rho_i_beta))*(dotX_rho_i) +
      + sy[2]*dnorm(X_rho_j_beta/sigma_jj)*pnorm(sy[1]*psi1(-X_rho_j_beta))*(dotX_rho_j)


    score[(i+1)/2,] <- c(matrix(deriv_beta/p,1,length(beta)),as.numeric(deriv_rho/p))

  }
  return(score)
  score<- -colSums(score)                               #score(theta)
  score_rho<-jacob*score[length(theta_tilde)]          #score_rho(theta_tilde)
  score<-c(score[1:(length(theta_tilde)-1)],score_rho) #score(theta_tilde)
  return(score)
}


f_rho_sar3 <- function (rho_tilde,W,X,eig,qu,Method_inv) {

  rho  <- eig$W[1]^(-1) + ((eig$W[2]^(-1)-eig$W[1]^(-1))/(1+exp(-rho_tilde))) #EQ C1 pag 33

  n=length(W[1,])
  n2<-1:n
  Q=is.finite(qu)
  tempo=proc.time()
  if (Q==FALSE)
  {if (Method_inv=='fast') {tempo=proc.time()
  if (rho==0) {
    Inv_rho=diag(n)
    X_rho=X
    dotX_rho=W%*%X
    Sigma_rho=diag(n)
    dotSigma_rho=W+t(W)}
  else {
    Inv_rho <- invIrW(x=W,rho=rho,method='solve',feasible=T)
    X_rho <- Inv_rho%*%X
    dotX_rho <- W%*%X_rho
    dotX_rho <- Inv_rho%*%X_rho #equazione 16 in proof (fondo pag 15)
    n_dispari<-n2[n2%%2!=0]
    Sigma_rho=matrix(0,n,n)
    dotSigma_rho=matrix(0,n,n)
    Inv_rho_s=Inv_rho+t(Inv_rho) #Invrho+ la sua trasposta
    for (i in n_dispari) {
      Inv_rho_g <- t(Inv_rho[i:(i+1),]) #prendo due RIGHE consecutive di Inv_rho ma tutte le colonne (perch??)
      Sigma_rho_g <- tcrossprod(t(Inv_rho_g)) #calcolo matrice varianza covarianza per il gruppo g (matrice 2x2)
      Sigma_rho[i,i]=Sigma_rho_g[1,1] #crea una copia , perch??
      Sigma_rho[i+1,i+1]=Sigma_rho_g[2,2]
      Sigma_rho[i,i+1]=Sigma_rho_g[1,2]
      Sigma_rho[i+1,i]=Sigma_rho_g[2,1]
      dot_Sigma_rho_g <- tcrossprod(Inv_rho_s,t(Inv_rho_g))
      dot_Sigma_rho_g <- crossprod(Inv_rho_g,dot_Sigma_rho_g)
      dotSigma_rho[i:(i+1),i:(i+1)]=as.matrix(dot_Sigma_rho_g)/rho-Sigma_rho[i:(i+1),i:(i+1)]/rho ##calcola a pezzi dotsigmarho?
    }
  }
  tempo=proc.time()-tempo
  }
    if (Method_inv=="chol") {tempo=proc.time()
    Inv_Sigma <- crossprod(diag(n)-rho*W)
    Sigma_rho <- chol2inv(chol(as(Inv_Sigma,"CsparseMatrix")))
    Inv_rho <- tcrossprod(Sigma_rho,diag(n)-rho*W)
    X_rho <- Inv_rho%*%X
    dotX_rho <- W%*%X_rho
    dotX_rho     <- Inv_rho%*%X_rho
    dotSigma_rho=matrix(0,n,n)
    if (rho==0) {dotSigma_rho=W+t(W)} else { n2<-1:n
    n_dispari<-n2[n2%%2!=0]
    for (i in n_dispari)
    {dotSigma_rho[i,i]=2*(Inv_rho[i,]%*%Sigma_rho[,i]-Sigma_rho[i,i])/rho
    dotSigma_rho[i+1,i+1]=2*(Inv_rho[i+1,]%*%Sigma_rho[,i+1]-Sigma_rho[i+1,i+1])/rho
    dotSigma_rho[i,i+1]=(Inv_rho[i,]%*%Sigma_rho[,i+1]+Sigma_rho[i,]%*%Inv_rho[i+1,]-2*Sigma_rho[i,i+1])/rho
    dotSigma_rho[i+1,i]=dotSigma_rho[i,i+1]}
    }
    tempo=proc.time()-tempo}
    if (Method_inv=="InvW") {tempo=proc.time()
    if (rho==0) {
      Inv_rho=diag(n)
      X_rho=X
      dotX_rho=W%*%X
      Sigma_rho=diag(n)
      dotSigma_rho=W+t(W)}
    else {
      Inv_rho      <- invIrW(x=W,rho=rho,method='solve',feasible=T)
      X_rho        <- Inv_rho%*%X
      Sigma_rho    <- tcrossprod(Inv_rho)                                   #directly approximate the inverse
      dotX_rho     <- W%*%X_rho
      dotX_rho     <- Inv_rho%*%dotX_rho
      #derivative X_rho with respect to rho   #equation C5
      dotSigma_rho=matrix(0,n,n)
      n2<-1:n
      n_dispari<-n2[n2%%2!=0]
      for (i in n_dispari)
      {dotSigma_rho[i,i]=2*(Inv_rho[i,]%*%Sigma_rho[,i]-Sigma_rho[i,i])/rho
      dotSigma_rho[i+1,i+1]=2*(Inv_rho[i+1,]%*%Sigma_rho[,i+1]-Sigma_rho[i+1,i+1])/rho
      dotSigma_rho[i,i+1]=(Inv_rho[i,]%*%Sigma_rho[,i+1]+Sigma_rho[i,]%*%Inv_rho[i+1,]-2*Sigma_rho[i,i+1])/rho
      dotSigma_rho[i+1,i]=dotSigma_rho[i,i+1]}
    }
    tempo=proc.time()-tempo}
  } else   {
    rhoW=diag(n)
    for (l in 1:qu) {rhoW=rho*rhoW%*%W}
    Inv_rho=rhoW
    X_rho        <- Inv_rho%*%X
    Sigma_rho    <- tcrossprod(Inv_rho)                                   #directly approximate the inverse
    dotX_rho     <- W%*%X_rho
    dotX_rho     <- Inv_rho%*%dotX_rho
    dotSigma_rho=matrix(0,n,n)
    if (rho==0) {dotSigma_rho=W+t(W)} else { n2<-1:n
    n_dispari<-n2[n2%%2!=0]
    for (i in n_dispari)
    {dotSigma_rho[i,i]=2*(Inv_rho[i,]%*%Sigma_rho[,i]-Sigma_rho[i,i])/rho
    dotSigma_rho[i+1,i+1]=2*(Inv_rho[i+1,]%*%Sigma_rho[,i+1]-Sigma_rho[i+1,i+1])/rho
    dotSigma_rho[i,i+1]=(Inv_rho[i,]%*%Sigma_rho[,i+1]+Sigma_rho[i,]%*%Inv_rho[i+1,]-2*Sigma_rho[i,i+1])/rho
    dotSigma_rho[i+1,i]=dotSigma_rho[i,i+1]}
    }
  }


  m_rho<-list()
  m_rho$Inv_rho      <- Inv_rho
  m_rho$X_rho        <- X_rho
  m_rho$Sigma_rho    <- Sigma_rho
  m_rho$dotX_rho     <- dotX_rho
  m_rho$dotSigma_rho <- dotSigma_rho
  m_rho$time<- c(tempo)
  return(m_rho)
}

prob.v0 <- function (k0,X_rho_i_beta,X_rho_j_beta,sigma_ii,sigma_rho_ij) {
  #ARRIVANO DA PROOF OF THEOREM 3.1
  #forse posso evitare di calcolare pmvnorm due volte
  #
  p_01 <-max(1e-135,pmvnorm(lower=c(-Inf,-X_rho_j_beta),upper=c(-X_rho_i_beta,Inf),sigma=sigma_rho_ij))
  #  if (p_01 <= 1e-135) {p_01 = 1e-135}

  p_00 <- max(1e-135,(1-pnorm(X_rho_i_beta/sigma_ii)) - p_01)
  #  if (p_00 <= 1e-135) {p_00 = 1e-135}

  p_11 <- max(1e-135,pmvnorm(lower=c(-X_rho_i_beta,-X_rho_j_beta),upper=c(Inf,Inf),sigma=sigma_rho_ij))
  #  if (p_11 <= 1e-135) {p_11 = 1e-135}

  p_10 <- max(1e-135,pnorm(X_rho_i_beta/sigma_ii) - p_11)
  #  if (p_10 <= 1e-135) {p_10 = 1e-135}

  p <- c(p_00,p_01,p_10,p_11)
  #p = p/sum(p)
  #print(p)
  return(p[k0])
}
