  mstn.varp 	<- function(Data, p, max.iter) {
  library(mvtnorm)  
  library(truncnorm)
  library(pracma)
  
  N 		<- nrow(Data)
  K 		<- ncol(Data)
  n		<- nrow(Data)
  create_matrix_A <- function(B, p) {
    N 		<- nrow(B)
    K 		<- ncol(B)
    
    # Calculate dimensions of matrix A
    rows_A 	<- N - p
    cols_A 	<- K * p
    
    # Create an empty matrix A
    A 		<- matrix(0, nrow = rows_A, ncol = cols_A)
    
    # Fill matrix A with values from matrix B
    for (i in 1:p) {
      start_row <- p - i + 1
      end_row 	<- N - i
      
      start_col <- (i - 1) * K + 1
      end_col 	<- i * K
      
      A[, start_col:end_col] <- B[start_row:end_row, ]
    }
    
    return(A)
  }
  
  X 		<- create_matrix_A(B=Data, p=p)
  X 		<- cbind(1, X)
  Y 		<- Data
    
  # initial values
  B.t 		<- matrix(c(rep(0,K), rep(diag(K),p)), nrow = K, ncol = (K*p + 1))
  C.t 		<- diag(K)
  BETA.t 	<- matrix(runif(K,0,1), nrow = K, ncol=1)
  S.t 		<- sqrtm(C.t)$B%*%BETA.t
  vt 		<- 4

  begin <- proc.time()[1]
  iter <- 0
  
  niter=1
  
  # Define log-likelihood function
  ll<-function(Y,X,vk,B.t,C.t,S.t) 
  {
    n 		<- nrow(Y)
    lik		<- array(0,(n-p))
    for(i in 1:(n-p)) {
      lik[i]	<- dmvt((as.numeric(as.matrix(Y[(i+p),])) - X[i,]%*%t(B.t)),sigma = C.t, df = vt)*ptruncnorm(t(S.t)%*%solve(sqrtm(C.t)$B)%*%t((as.numeric (as.matrix(Y[(i+p),])) - X[i,]%*%t(B.t))))
    }
    
    sum(lik) 
  }
  # Compute the log-likelihood for the initial values
  lltm1		<-ll(Y,X,vt,B.t,C.t,S.t) 
  
  repeat {  
    # E-Step1
    iter 	<- iter + 1
    Tau 	<- array(0,n)
    Kay 	<- array(0,n)
    Gamma 	<- array(0,n)
    Dnorm 	<- array(0,n)
    Pnorm 	<- array(0,n)
    for(i in 1:(n-p)) {
      Tau[i] 	<- (vt+K)/(vt+t(solve(sqrtm(C.t)$B)%*%(Y[(i+p),] - B.t%*%X[i,]))%*%(solve(sqrtm(C.t)$B)%*%(Y[(i+p),] - B.t%*%X[i,])))
      Kay[i] 	<- digamma(0.5*(vt+K)) - log(0.5*(vt+t(solve(sqrtm(C.t)$B)%*%(Y[(i+p),] - B.t%*%X[i,]))%*%(solve(sqrtm(C.t)$B)%*%(Y[(i+p),] - B.t%*%X[i,]))))
      Dnorm[i] 	<- dnorm(t(S.t)%*%(solve(sqrtm(C.t)$B)%*%(Y[(i+p),] - B.t%*%X[i,])))
      Pnorm[i] 	<- pnorm(t(S.t)%*%(solve(sqrtm(C.t)$B)%*%(Y[(i+p),] - B.t%*%X[i,])))
    }
    Dnorm[Dnorm == (0)] = 1.097221e-314
    Pnorm[Pnorm == (0)] = 5.725571e-300
    
    for(i in 1:(n-p)) {
      Gamma[i] 	<- t(S.t)%*%(solve(sqrtm(C.t)$B)%*%(Y[(i+p),] - B.t%*%X[i,])) + Dnorm[i]/Pnorm[i]
    }
    W 		<- diag(Tau[-1])
    
    # M-Step1

    R_0 		<- array(NA,c(K,K,p))
    W0 		<- matrix(0, nrow=K, ncol=K*p+1)
    r0 		<- t(Y[(p+1):n,])%*%Tau[-c(1:p)]
    for(j in 1:p) {
    R_0[,,j] 	<- t(Y[(p+1):n,])%*%diag(Tau[-c(1:p)])%*%Y[(p+1-j):(n-j),]
    }
    W0[, 1] 	<- r0
    W0[, 2:(K * p + 1)] <- as.matrix(R_0)

    r_i 	<- array(NA,c(K,1,p))
    for(j in 1:p) {
    r_i[,,j]	<- t(Y[(p+1-j):(n-j),])%*%Tau[-c(1:p)]
    }
    
    R_ij	<- array(NA,c(K,K,p,p))
    for(i in 1:p) {
    for(j in 1:p) {
    R_ij[,,i,j]	<- t(Y[(p+1-i):(n-i),])%*%diag(Tau[-c(1:p)])%*%Y[(p+1-j):(n-j),] 
    }
    }  
    
    r 				<- sum(Tau[(p+1):n])
    
    W1 				<- matrix(0, nrow = K*p+1, ncol = K*p+1)
    W1[1,1] 			<- r 
    W1[2:(K*p+1),1] 		<- as.matrix(r_i) 
    W1[1,2:(K*p+1)] 		<- t(as.matrix(r_i))
    
    block_size 			<- K  # Size of each block
    for (i in 1:p) {
    for (j in 1:p) {
    row_start 			<- (i - 1) * block_size + 2
    row_end 			<- i * block_size + 1
    col_start 			<- (j - 1) * block_size + 2
    col_end 			<- j * block_size + 1
    
    W1[row_start:row_end, col_start:col_end] <- R_ij[, , i, j]
    }
    }

    G1 <- t(Y[(p+1):n,])%*%Gamma[-c(1:p)]
    G2 <- t(X[1:(n-p),])%*%Gamma[-c(1:p)]
    R00 <- t(Y[(p+1):n,])%*%diag(Tau[-c(1:p)])%*%Y[(p+1):n,]
    
    B.t1 <- 4*solve(solve(C.t,tol=1e-40)+S.t%*%t(S.t),tol=1e-40)%*%S.t%*%t(G2)%*%solve(t(W1)+W1,tol=1e-44)+2*W0%*%solve(t(W1)+W1,tol=1e-46)
    BETA.t1 <- 4*(solve(2*R00 -2*W0%*%t(B.t1) - 2*B.t1%*%t(W0) + B.t1%*%W1%*%t(B.t1) + B.t1%*%t(W1)%*%t(B.t1),tol=1e-40))%*%(B.t1%*%G2 - G1)
    C.t1 <- (R00 - W0%*%t(B.t1) - B.t1%*%t(W0) + B.t1%*%W1%*%t(B.t1))/(n-p)
    Yhat <- X%*%t(B.t)
    #C.t1 <- (t(Y - Yhat)%*%(Y - Yhat))/(n - p)
    S.t1 <- sqrtm(C.t1)$B%*%BETA.t1
    
    u.t1 <- (as.numeric(as.matrix(Y[(p+1),])) - X%*%t(B.t1))
    df.ll <- function(vt,u.t1) 
    {
      df.lik <- c()
      for(i in 1:(n-p)) {
        df.lik[i]<- dmvt(u.t1[i,], sigma = C.t1, df = vt)*ptruncnorm(t(S.t1)%*%solve(sqrtm(C.t1)$B)%*%(u.t1[i,]))
      }
      sum(df.lik) 
    }
    v.t1 = optim(vt, u.t1 = u.t1, df.ll, lower = 0.5121, upper = 5.5, control=list(fnscale=-1), method = "Brent")$par
    #u.t1 = as.matrix(u.t1)
    # Update parameter values
    B.t<-B.t1
    C.t<-C.t1
    S.t<-S.t1
    vt <-v.t1
    
    # Compute log-likelihood using current estimates
    
    llt<-ll(Y, X, vt, B.t, C.t, S.t)
    
    # Stop if converged
    #cat(llt, "\n")
    
    if(abs(lltm1-llt) < 1e-04  | iter == max.iter) break
    
    lltm1<-llt
    niter=niter+1
    
  }
  Yhat <- X%*%t(B.t)
  T <- n-p
  f1 <- 1 + 76/T
  f2 <- 1 - 76/T
  FPE <- det(C.t)*(f1/f2)
  AIC <- log(det(C.t)) + 2*(p*K^2 + K + 1)/T 
  HQ <- log(det(C.t)) + 2*log(log(T))*p*K^2/T
  SC <- log(det(C.t)) + log(T)*p*K^2/T
  
  results <-  list("Phi"   	    	= B.t,
                   "Sigma"  	    	= C.t, 
                   "Alpha"  	    	= S.t, 
                   "DF"     	    	= vt,
                   "Loglik"         = lltm1,
                   "FPE"        	= FPE,
                   "AIC"        	= AIC,
                   "HQ"         	= HQ,
                   "SC"         	= SC,
                   "iteration"  	= niter)
			 #"Residuals"	= u.t1)
  return(results)
}
