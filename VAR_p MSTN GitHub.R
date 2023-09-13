rm(list=ls(all=TRUE))
library(mvtnorm)
#library(expm)
library(truncnorm)
library(pracma)

# Data generation
VAR.MSTN_simulator <- function(K,p,n,v) {
  MSTN.sim<-function(n,v,K) {
    Sigma = matrix(0,K,K)
    for(i in 1:K) {
      for(j in 1:K) {
        Sigma[i,j] = 0.64^abs(i-j)
      }
    }
    Alpha<-as.matrix(c(s1=0.2, s2=0.2, s3=0.2, s4=0.2, s5=0.2))
    mu<-c(0,0,0,0,0)
    U1 <- rnorm(1,0,1)
    U2 <- rmvnorm(1,mu,diag(K))
    tau <- rgamma(1,0.5*v,0.5*v)
    error <- sqrtm(Sigma)$B%*%((Alpha*abs(U1))/as.numeric(sqrt(tau*(tau + t(Alpha)%*%Alpha)))) + solve(sqrtm(tau*diag(K) + Alpha%*%t(Alpha))$B)%*%t(U2)
    return(error)
  }
  # Generate coefficient matrices
  A.1 <- matrix(c(.25, .010, -.009,  .091, -.009, .036, .25, .015, -.079, .080, -.051, -.092, .25, -.034, .091, .078, .039, .028, 0.25, .099, .031, .042, .009, .019, 0.25), K) # Coefficient matrix of lag 1
  A.2 <- matrix(c(.15, -.042, -.071, .093, .080, .038, 0.15, .059, -.095, -.004, .052, -.057, 0.15, -.036, -.054, -.071, -.017, -.017, 0.15, -.026, -.070, -.072, -.053, -.007, 0.15), K) # Coefficient matrix of lag 2
  mu <- rbind(-.056,  .086, -.109, -.014,  .072)
  A <- cbind(mu, A.1, A.2) # Companion form of the coefficient matrices
  
  
  # Generate series
  series <- matrix(0, K, n + 2*p) # Raw series with zeros
  for (i in (p + 1):(n + 2*p)){ # Generate series with e ~ N(0,0.5)
    series[, i] <- mu + A.1%*%series[, i-1] + A.2%*%series[, i-2] + MSTN.sim(1,v,K)
  }
  series <- ts(t(series[, -(1:p)])) # Convert to time series format
  return(series)
}

set.seed(1212)
Data <- array(0,c(200,5,100))
for(i in 1:100) {
  Data[, , i] <- VAR.MSTN_simulator(K=5, p=2, n=198, v=2.5)
}


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
      Kay[i] 	<- digamma(0.5*(vt+K)) - log(0.5*(vt+t(solve(sqrtm(C.t)$B)%*%(Y[(i+p),] - B.t%*%X[i,]))%*%(solve(sqrtm(C.t)$B)%*%(Y[(i+p),] 
                                                                                                                                - B.t%*%X[i,]))))
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
    
    R_0 	<- array(NA,c(K,K,p))
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
    
    r 		<- sum(Tau[(p+1):n])
    
    W1 				<- matrix(0, nrow = K*p+1, ncol = K*p+1)
    W1[1,1] 			<- r 
    W1[2:(K*p+1),1] 		<- as.matrix(r_i) 
    W1[1,2:(K*p+1)] 		<- t(as.matrix(r_i))
    
    block_size 			<- K  # Size of each block
    for (i in 1:p) {
      for (j in 1:p) {
        row_start 			<- (i - 1) * block_size + 2
        row_end 			  <- i * block_size + 1
        col_start 			<- (j - 1) * block_size + 2
        col_end 			  <- j * block_size + 1
        
        W1[row_start:row_end, col_start:col_end] <- R_ij[, , i, j]
      }
    }
    
    G1 <- t(Y[(p+1):n,])%*%Gamma[-c(1:p)]
    G2 <- t(X[1:(n-p),])%*%Gamma[-c(1:p)]
    R00 <- t(Y[(p+1):n,])%*%diag(Tau[-c(1:p)])%*%Y[(p+1):n,]
    
    B.t1 <- 4*solve(solve(C.t,tol=1e-20)+S.t%*%t(S.t),tol=1e-20)%*%S.t%*%t(G2)%*%solve(t(W1)+W1,tol=1e-24)+2*W0%*%solve(t(W1)+W1,
                                                                                                                        tol=1e-26)
    BETA.t1 <- 4*(solve(2*R00 -2*W0%*%t(B.t1) - 2*B.t1%*%t(W0) + B.t1%*%W1%*%t(B.t1) + B.t1%*%t(W1)%*%t(B.t1),tol=1e-20))%*%(B.t1%*%G2 - G1)
    C.t1 <- (R00 - W0%*%t(B.t1) - B.t1%*%t(W0) + B.t1%*%W1%*%t(B.t1))/(n-p)
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
    v.t1 = optim(vt, u.t1 = u.t1, df.ll, lower = 0.09, upper = 5.5, control=list(fnscale=-1), method = "Brent")$par
    
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
  
  A.1 <- matrix(c(.25, .010, -.009,  .091, -.009, .036, .25, .015, -.079, .080, -.051, -.092, .25, -.034, .091, .078, .039, .028, 0.25, .099, .031, .042, .009, .019, 0.25), K) # Coefficient matrix of lag 1
  A.2 <- matrix(c(.15, -.042, -.071, .093, .080, .038, 0.15, .059, -.095, -.004, .052, -.057, 0.15, -.036, -.054, -.071, -.017, -.017, 0.15, -.026, -.070, -.072, -.053, -.007, 0.15), K) # Coefficient matrix of lag 2
  mu <- rbind(-.056,  .086, -.109, -.014,  .072)
  A <- cbind(mu, A.1, A.2) # Companion form of the coefficient matrices
  norm_phi <- norm((A - B.t), type = "F")
  norm_phi <- norm_phi^2
  
  Sigma = matrix(0,K,K)
  for(i in 1:K) {
    for(j in 1:K) {
      Sigma[i,j] = 0.64^abs(i-j)
    }
  }
  norm_Sigma <- norm((Sigma - C.t), type = "F")
  norm_Sigma <- norm_Sigma^2
  Alpha<-as.matrix(c(s1=0.2, s2=0.2, s3=0.2, s4=0.2, s5=0.2))
  norm_Alpha <- norm((Alpha - S.t), type = "F")
  norm_Alpha <- norm_Alpha^2
  df <- 3.5
  norm_df <- (df - vt)^2
  
  Yhat <- X%*%t(B.t)
  T <- n-p
  f1 <- 1 + 76/T
  f2 <- 1 - 76/T
  FPE <- det(C.t)*(f1/f2)
  AIC <- log(det(C.t)) + 2*(p*K^2 + K + 1)/T 
  HQ <- log(det(C.t)) + 2*log(log(T))*p*K^2/T
  SC <- log(det(C.t)) + log(T)*p*K^2/T
  
  results <-    c( "Phinorm"        = norm_phi,
                   "Sigmanorm"      = norm_Sigma, 
                   "Alphanorm"      = norm_Alpha, 
                   "DFnorm"         = norm_df,
                   "Loglik"         = lltm1,
                   "FPE"        	  = FPE,
                   "iteration"  	  = niter)
  return(as.numeric(results))
}

#Output1 <- mstn.varp(Data=Data, p=2, max.iter=2000)
#Output1

library(future)
future(seed=T)
plan(multiprocess)

Res1 %<-% mstn.varp(Data[,,1], p=2, max.iter=2000)
Res2 %<-% mstn.varp(Data[,,2], p=2, max.iter=2000)
Res3 %<-% mstn.varp(Data[,,3], p=2, max.iter=2000)
Res4 %<-% mstn.varp(Data[,,4], p=2, max.iter=2000)
Res5 %<-% mstn.varp(Data[,,5], p=2, max.iter=2000)
Res6 %<-% mstn.varp(Data[,,6], p=2, max.iter=2000)
Res7 %<-% mstn.varp(Data[,,7], p=2, max.iter=2000)
Res8 %<-% mstn.varp(Data[,,8], p=2, max.iter=2000)
Res9 %<-% mstn.varp(Data[,,9], p=2, max.iter=2000)
Res10 %<-% mstn.varp(Data[,,10], p=2, max.iter=2000)
Res11 %<-% mstn.varp(Data[,,11], p=2, max.iter=2000)
Res12 %<-% mstn.varp(Data[,,12], p=2, max.iter=2000)
Res13 %<-% mstn.varp(Data[,,13], p=2, max.iter=2000)
Res14 %<-% mstn.varp(Data[,,14], p=2, max.iter=2000)
Res15 %<-% mstn.varp(Data[,,15], p=2, max.iter=2000)
Res16 %<-% mstn.varp(Data[,,16], p=2, max.iter=2000)
Res17 %<-% mstn.varp(Data[,,17], p=2, max.iter=2000)
Res18 %<-% mstn.varp(Data[,,18], p=2, max.iter=2000)
Res19 %<-% mstn.varp(Data[,,19], p=2, max.iter=2000)
Res20 %<-% mstn.varp(Data[,,20], p=2, max.iter=2000)
Res21 %<-% mstn.varp(Data[,,21], p=2, max.iter=2000)
Res22 %<-% mstn.varp(Data[,,22], p=2, max.iter=2000)
Res23 %<-% mstn.varp(Data[,,23], p=2, max.iter=2000)
Res24 %<-% mstn.varp(Data[,,24], p=2, max.iter=2000)
Res25 %<-% mstn.varp(Data[,,25], p=2, max.iter=2000)
Res26 %<-% mstn.varp(Data[,,26], p=2, max.iter=2000)
Res27 %<-% mstn.varp(Data[,,27], p=2, max.iter=2000)
Res28 %<-% mstn.varp(Data[,,28], p=2, max.iter=2000)
Res29 %<-% mstn.varp(Data[,,29], p=2, max.iter=2000)
Res30 %<-% mstn.varp(Data[,,30], p=2, max.iter=2000)
Res31 %<-% mstn.varp(Data[,,31], p=2, max.iter=2000)
Res32 %<-% mstn.varp(Data[,,32], p=2, max.iter=2000)
Res33 %<-% mstn.varp(Data[,,33], p=2, max.iter=2000)
Res34 %<-% mstn.varp(Data[,,34], p=2, max.iter=2000)
Res35 %<-% mstn.varp(Data[,,35], p=2, max.iter=2000)
Res36 %<-% mstn.varp(Data[,,36], p=2, max.iter=2000)
Res37 %<-% mstn.varp(Data[,,37], p=2, max.iter=2000)
Res38 %<-% mstn.varp(Data[,,38], p=2, max.iter=2000)
Res39 %<-% mstn.varp(Data[,,39], p=2, max.iter=2000)
Res40 %<-% mstn.varp(Data[,,40], p=2, max.iter=2000)
Res41 %<-% mstn.varp(Data[,,41], p=2, max.iter=2000)
Res42 %<-% mstn.varp(Data[,,42], p=2, max.iter=2000)
Res43 %<-% mstn.varp(Data[,,43], p=2, max.iter=2000)
Res44 %<-% mstn.varp(Data[,,44], p=2, max.iter=2000)
Res45 %<-% mstn.varp(Data[,,45], p=2, max.iter=2000)
Res46 %<-% mstn.varp(Data[,,46], p=2, max.iter=2000)
Res47 %<-% mstn.varp(Data[,,47], p=2, max.iter=2000)
Res48 %<-% mstn.varp(Data[,,48], p=2, max.iter=2000)
Res49 %<-% mstn.varp(Data[,,49], p=2, max.iter=2000)
Res50 %<-% mstn.varp(Data[,,50], p=2, max.iter=2000)
Res51 %<-% mstn.varp(Data[,,51], p=2, max.iter=2000)
Res52 %<-% mstn.varp(Data[,,52], p=2, max.iter=2000)
Res53 %<-% mstn.varp(Data[,,53], p=2, max.iter=2000)
Res54 %<-% mstn.varp(Data[,,54], p=2, max.iter=2000)
Res55 %<-% mstn.varp(Data[,,55], p=2, max.iter=2000)
Res56 %<-% mstn.varp(Data[,,56], p=2, max.iter=2000)
Res57 %<-% mstn.varp(Data[,,57], p=2, max.iter=2000)
Res58 %<-% mstn.varp(Data[,,58], p=2, max.iter=2000)
Res59 %<-% mstn.varp(Data[,,59], p=2, max.iter=2000)
Res60 %<-% mstn.varp(Data[,,60], p=2, max.iter=2000)
Res61 %<-% mstn.varp(Data[,,61], p=2, max.iter=2000)
Res62 %<-% mstn.varp(Data[,,62], p=2, max.iter=2000)
Res63 %<-% mstn.varp(Data[,,63], p=2, max.iter=2000)
Res64 %<-% mstn.varp(Data[,,64], p=2, max.iter=2000)
Res65 %<-% mstn.varp(Data[,,65], p=2, max.iter=2000)
Res66 %<-% mstn.varp(Data[,,66], p=2, max.iter=2000)
Res67 %<-% mstn.varp(Data[,,67], p=2, max.iter=2000)
Res68 %<-% mstn.varp(Data[,,68], p=2, max.iter=2000)
Res69 %<-% mstn.varp(Data[,,69], p=2, max.iter=2000)
Res70 %<-% mstn.varp(Data[,,70], p=2, max.iter=2000)
Res71 %<-% mstn.varp(Data[,,71], p=2, max.iter=2000)
Res72 %<-% mstn.varp(Data[,,72], p=2, max.iter=2000)
Res73 %<-% mstn.varp(Data[,,73], p=2, max.iter=2000)
Res74 %<-% mstn.varp(Data[,,74], p=2, max.iter=2000)
Res75 %<-% mstn.varp(Data[,,75], p=2, max.iter=2000)
Res76 %<-% mstn.varp(Data[,,76], p=2, max.iter=2000)
Res77 %<-% mstn.varp(Data[,,77], p=2, max.iter=2000)
Res78 %<-% mstn.varp(Data[,,78], p=2, max.iter=2000)
Res79 %<-% mstn.varp(Data[,,79], p=2, max.iter=2000)
Res80 %<-% mstn.varp(Data[,,80], p=2, max.iter=2000)
Res81 %<-% mstn.varp(Data[,,81], p=2, max.iter=2000)
Res82 %<-% mstn.varp(Data[,,82], p=2, max.iter=2000)
Res83 %<-% mstn.varp(Data[,,83], p=2, max.iter=2000)
Res84 %<-% mstn.varp(Data[,,84], p=2, max.iter=2000)
Res85 %<-% mstn.varp(Data[,,85], p=2, max.iter=2000)
Res86 %<-% mstn.varp(Data[,,86], p=2, max.iter=2000)
Res87 %<-% mstn.varp(Data[,,87], p=2, max.iter=2000)
Res88 %<-% mstn.varp(Data[,,88], p=2, max.iter=2000)
Res89 %<-% mstn.varp(Data[,,89], p=2, max.iter=2000)
Res90 %<-% mstn.varp(Data[,,90], p=2, max.iter=2000)
Res91 %<-% mstn.varp(Data[,,91], p=2, max.iter=2000)
Res92 %<-% mstn.varp(Data[,,92], p=2, max.iter=2000)
Res93 %<-% mstn.varp(Data[,,93], p=2, max.iter=2000)
Res94 %<-% mstn.varp(Data[,,94], p=2, max.iter=2000)
Res95 %<-% mstn.varp(Data[,,95], p=2, max.iter=2000)
Res96 %<-% mstn.varp(Data[,,96], p=2, max.iter=2000)
Res97 %<-% mstn.varp(Data[,,97], p=2, max.iter=2000)
Res98 %<-% mstn.varp(Data[,,98], p=2, max.iter=2000)
Res99 %<-% mstn.varp(Data[,,99], p=2, max.iter=2000)
Res100 %<-% mstn.varp(Data[,,100], p=2, max.iter=2000)

B4=rbind(Res1,Res2,Res3,Res4,Res5,Res6,Res7,Res8,Res9,Res10,
         Res11,Res12,Res13,Res14,Res15,Res16,Res17,Res18,Res19,Res20,
         Res21,Res22,Res23,Res24,Res25,Res26,Res27,Res28,Res29,Res30,
         Res31,Res32,Res33,Res34,Res35,Res36,Res37,Res38,Res39,Res40,
         Res41,Res42,Res43,Res44,Res45,Res46,Res47,Res48,Res49,Res50,
         Res51,Res52,Res53,Res54,Res55,Res56,Res57,Res58,Res59,Res60,
         Res61,Res62,Res63,Res64,Res65,Res66,Res67,Res68,Res69,Res70,
         Res71,Res72,Res73,Res74,Res75,Res76,Res77,Res78,Res79,Res80,
         Res81,Res82,Res83,Res84,Res85,Res86,Res87,Res88,Res89,Res90,
         Res91,Res92,Res93,Res94,Res95,Res96,Res97,Res98,Res99,Res100)
Result4 <- colMeans(B4)
names(Result4) <- c("Phinorm", "Sigmanorm", "Alphanorm", "DFnorm", "Loglik", "FPE", "iter")
Result4

