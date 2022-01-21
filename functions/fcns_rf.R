rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

log_p_alpha=function(Alpha,A,X,b_eta,G1,Tsq,O,E,K,m,Gden){
  #This Alpha is a vector without the first 0
  Ind=matrix(0,nrow=length(E),ncol=K)
  
  for (i in 1: length(E)){
    Ind[i,E[i]]=1;
  }
  
  if (length(A)==length(Alpha)){ #intercept only
    log_p=
      sum(
        rowSums(Ind*(X%*%t(c(0,Alpha)) + cbind(0,b_eta[G1,])- log(rowSums(exp(X%*%t(c(0,Alpha))+cbind(0,b_eta[G1,]) )))))
      ) - 0.5*(t(Alpha)%*%O%*%Alpha)/Tsq
  }else{
    Xm=X[,m]; X1=X[,-m];
    
    if(dim(A)[1]==2){ #int and X1 only
      AA=cbind(0,t(A[-m,]));
    }else{ #more covariates; note: I don't think this is needed here		
      AA=cbind(rep(0,(dim(A)[1]-1)),A[-m,]);
    }		
    # alpha_log_p = c()
    # for(i in 1:nrow(AA)){
    #   alpha_log_p = c(alpha_log_p, -0.5*(t(AA[i,2:K])%*%O%*%AA[i,2:K])/Tsq[i])
    # }
    log_p=sum(
      rowSums(Ind*((Xm%*%t(c(0,Alpha))+(X1%*%AA)+cbind(0,b_eta[G1,]))-
                     log(rowSums(exp(Xm%*%t(c(0,Alpha))+(X1%*%AA)+cbind(0,b_eta[G1,]))))))
    ) + sum(dnorm(x = Alpha, mean = 0, sd = 2.5, log = T))
    #+ sum(dcauchy(x = Alpha,scale=2.5, log = T))
    # log_p = sum(log(
    #   rowSums(
    #     Gden*exp(Xm%*%t(c(0,Alpha))+(X1%*%AA)+cbind(0,b_eta[G1,]))/rowSums(exp(Xm%*%t(c(0,Alpha))+(X1%*%AA)+cbind(0,b_eta[G1,]))) 
    #   )
    # ))- 0.5*(t(Alpha)%*%O%*%Alpha)/Tsq[m]+sum(alpha_log_p)
  }
  return(log_p)
}



## Draw Alpha,log-scale MH alg:multiple-try;
Draw_Alpha_MMH=function(Alpha,A,X1,b_eta,G1,Tau_sq,df,Multiple,Omega_star,SSA,SSAP,Eta,K,m,iter,burnIn, cov = NULL, Gden)
{					
		log_p_Alpha_star=rep(0,Multiple)
		log_p_Alpha_2star=rep(0,Multiple)
		p=rep(0,Multiple)
		den=0;num=0;

		# if(is.null(cov)){
		alpha_opt <-tryCatch({optim(Alpha,A=A,X=X1,b_eta=b_eta,G1=G1,Tsq=Tau_sq,O=Omega_star,E=Eta,K=K,m=m,Gden=Gden,
		                      log_p_alpha,method="Nelder-Mead",
		                      hessian=TRUE,control=list(maxit=10,fnscale=-1))},error=function(e) NULL)
		
		if(is.null(alpha_opt)){
		  alpha = Alpha
		  sigma = cov
		  Alpha_star=t(rmvt(n=Multiple,delta = alpha,sigma=sigma,df=df))
		}else{
		  sigma=tryCatch(solve(-alpha_opt$hessian), error=function(e) NULL)
		  
		  if(is.null(sigma)){
		    alpha=alpha_opt$par	
		    sigma = cov
		    Alpha_star=t(rmvt(n=Multiple,delta = alpha,sigma=sigma,df=df))
		  }else{
		    alpha=alpha_opt$par	
		    
		    if(iter > burnIn & sample(c(0,1),1,prob=c(1-SSAP,SSAP))==1){
		      diag(sigma)=diag(sigma)*SSA
		    }
		    if(!isSymmetric(sigma)){
		      sigma = diag(diag(sigma));
		    }
		    Alpha_star=t(rmvt(n=Multiple,delta = alpha,sigma=sigma,df=df))
		  }
		  
		}
		  
		# }else{
		#   alpha = Alpha
		#   sigma = cov
		#   Alpha_star=t(MASS::mvrnorm(n=Multiple,alpha,Sigma=sigma))
		# }
		
			
			for (i in 1:Multiple){
				log_p_Alpha_star[i]=log_p_alpha(Alpha_star[,i],A,X1,b_eta,G1,Tau_sq,Omega_star,Eta,K,m, Gden)
			}
			d_p_Alpha_star = dmvt(t(Alpha_star), delta = alpha,sigma=sigma,df=df, log =T)
      if(any(is.na(log_p_Alpha_star))){
        print("overflow in Alpha density")
      }
			# browser()
			max_p_alpha = median(log_p_Alpha_star)
			  # log_p_Alpha_star[which(abs(log_p_Alpha_star) == max(abs(log_p_Alpha_star)),log_p_Alpha_star)]
			#control overfloat, -max(log_p_Alpha_star);
			p=exp(log_p_Alpha_star-d_p_Alpha_star-max_p_alpha)/sum(exp(log_p_Alpha_star - d_p_Alpha_star-max_p_alpha))
			#in case there is still overfloat;
			p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
			j=sample(c(1:Multiple),1,prob=p);							
			
			Alpha_2star=t(rmvt(n=Multiple-1,delta = Alpha_star[,j],sigma=sigma,df=df))
			Alpha_2star <-cbind(Alpha_2star,Alpha)
			
			for (i in 1:Multiple){
				log_p_Alpha_2star[i]=log_p_alpha(Alpha_2star[,i],A,X1,b_eta,G1,Tau_sq,Omega_star,Eta,K,m,Gden)
			}
			d_p_Alpha_2star = dmvt(t(Alpha_2star), delta = Alpha_star[,j],sigma=sigma,df=df,log =T)
			if(any(is.na(log_p_Alpha_2star))){
			  print("overflow in Alpha density")
			}
			
			#control overfloat
			num=sum(exp(log_p_Alpha_star -d_p_Alpha_star-max_p_alpha))
			den=sum(exp(log_p_Alpha_2star -d_p_Alpha_2star-max_p_alpha))
											
			rho=min(1,num/den)
			
			#in case overfloat again
			if(is.na(rho)) {
			  print("alpha rho is na")
			  rho=0.05};
			
			accp=0;
			u=runif(1)
			if(u<rho){
				Alpha=Alpha_star[,j]
				accp=1;
			}
		
		
		return(list(par=Alpha, accp=accp))
}

## Draw Alpha,log-scale MH alg:multiple-try;
Draw_Alpha_PG=function(Alpha_k,A,X1,b_eta,G1,Eta,k, Sigma_alpha)
{					
  
  # Ind=matrix(0,nrow=length(Eta),ncol=K)
  # 
  # for (i in 1: length(Eta)){
  #   Ind[i,Eta[i]]=1;
  # }
  X1_no_int = X1[,-1]
  #Parameters
  kappa = as.numeric(Eta==k)-1/2
  C_k = log(rowSums(exp(X1 %*% A + cbind(0,b_eta[G1,-(k-1)]))))-Alpha_k[1]-b_eta[G1,(k-1)]
  S_k = X1_no_int %*% Alpha_k[-1] - C_k
  Omega = c()
  for(j in 1: length(S_k)){
    Omega[j] = BayesLogit::rpg(1, h = 1, z = S_k[j])
  }
  Omega = diag(Omega)
  V_k = solve(t(X1_no_int) %*% Omega %*% X1_no_int + solve(Sigma_alpha))
  m_k = V_k %*% (t(X1_no_int)%*%(kappa+Omega%*%C_k))
  Alpha_k = mvrnorm(1, mu = m_k, Sigma = V_k)
  return(list(par=Alpha_k, accp=1))
}

Draw_Alpha_PG_full=function(Alpha_k,A,X1,b_eta,G1,Eta,k, Sigma_alpha)
{					
  
  # Ind=matrix(0,nrow=length(Eta),ncol=K)
  # 
  # for (i in 1: length(Eta)){
  #   Ind[i,Eta[i]]=1;
  # }
  X1_no_int = X1[,-1]
  #Parameters
  kappa = as.numeric(Eta==k)-1/2
  C_k = log(rowSums(exp(X1 %*% A + cbind(0,b_eta[G1,-(k-1)]))))-b_eta[G1,(k-1)]
  S_k = X1 %*% Alpha_k - C_k
  Omega = c()
  for(j in 1: length(S_k)){
    Omega[j] = BayesLogit::rpg(1, h = 1, z = S_k[j])
  }
  Omega = diag(Omega)
  V_k = solve(t(X1) %*% Omega %*% X1 + solve(Sigma_alpha))
  m_k = V_k %*% (t(X1)%*%(kappa+Omega%*%C_k))
  Alpha_k = mvrnorm(1, mu = m_k, Sigma = V_k)
  return(list(par=Alpha_k, accp=1))
}

Draw_Gamma_log_M = function(y, X, G, SG, Phi,sigma_b_delta_sq,same_length_G){
  T_1 = list()
  T_2 = list()
  V = list()
  i = 1
  if(same_length_G){
    indx = G == 1 
    n = sum(indx)
    V_i = diag(n)+rep(1,n)%*%t(rep(1,n))*sigma_b_delta_sq
    V_i_inv = solve(V_i)
    for(g in unique(G)){
      indx = G == g
      X_i = X[indx,]
      y_i = y[indx]
      V[[i]]=V_i
      if(n>1){
        T_1[[i]] = t(X_i)%*%V_i_inv%*%X_i
        T_2[[i]] = t(X_i)%*%V_i_inv%*%y_i
      }else{
        T_1[[i]] = X_i%*%V_i_inv%*%t(X_i)
        T_2[[i]] = X_i%*%V_i_inv%*%y_i
      }
      i = i+1
    }
  }else{
    for(g in unique(G)){
      indx = G == g
      n = sum(indx)
      X_i = X[indx,]
      y_i = y[indx]
      V_i = diag(n)+rep(1,n)%*%t(rep(1,n))*sigma_b_delta_sq
      V[[i]]=V_i
      V_i_inv = solve(V_i)
      if(n>1){
        T_1[[i]] = t(X_i)%*%V_i_inv%*%X_i
        T_2[[i]] = t(X_i)%*%V_i_inv%*%y_i
        
      }else{
        T_1[[i]] = X_i%*%V_i_inv%*%t(X_i)
        T_2[[i]] = X_i%*%V_i_inv%*%y_i
      }
      i = i+1
    }
  }
  
  var = solve(solve(SG)+Reduce("+",T_1))
  mean = var%*%Reduce("+",T_2)
  Gamma = MASS::mvrnorm(n=1,mu = mean,Sigma=var)
  return(list(par=Gamma, accp=1, V = V))
}

# Draw_Gamma_log_M = function(y, X, G, SG, Phi,b_delta){
#   SG_inv = solve(SG)
#   y_prime = y - b_delta[G]
#   var = solve(SG_inv+t(X)%*%X)
#   mean = var%*%t(X)%*%y_prime
#   Gamma = MASS::mvrnorm(n=1,mu = mean,Sigma=var)
#   return(list(par=Gamma, accp=1))
# }
# log_p_gamma=function(Gamma,Z,X,b_delta,G,SG,phi){
# 	
# 	log_p=sum(phi*(X%*%Gamma+b_delta[G])-log(1+exp(X%*%Gamma+b_delta[G])))-0.5*t(Gamma)%*%solve(SG)%*%Gamma
# 
# 	
# 	return(log_p)
# }
# 
# ##Draw Gamma, log scale, multiple-try MH;
# Draw_Gamma_log_M=function(Gamma,Z,X,b_delta,G,Phi,df,Sigma_Gamma,SSG,Multiple,useCov, Var = NULL)
#   
# {
#   
#   log_p_Gamma_star=rep(0,Multiple)
#   log_p_Gamma_2star=rep(0,Multiple)
#   p=rep(0,Multiple)
#   
#   
#   if(is.null(Var)){
#     if (useCov==1){
#       gamma_opt <-tryCatch(optim(Gamma,Z=Z,X=X,b_delta=b_delta,G=G,SG=Sigma_Gamma,phi=Phi,log_p_gamma,method="Nelder-Mead",
#                                  hessian=TRUE,control=list(maxit=10,fnscale=-1)), error=function(e) NULL)
#     }else{
#       gamma_opt <-tryCatch(optim(Gamma,Z=Z,X=X,b_delta=b_delta,G=G,SG=Sigma_Gamma,phi=Phi,log_p_gamma,method="Brent",
#                                  hessian=TRUE,lower=-1000-abs(Gamma),upper=1000+abs(Gamma),control=list(maxit=10,fnscale=-1)), error=function(e) NULL)
#       
#     }
#     if(is.null(gamma_opt)){
#       return(list(par=rnorm(M,0,1), accp=0))
#     }
#     
#     gamma=gamma_opt$par
#     Var=solve(-gamma_opt$hessian)
#     diag(Var)=diag(Var)*SSG
#     Gamma_star=t(rmvt(n=Multiple,gamma,sigma=Var,df=df)) #use the optim algorithm
#   }else{
#     diag(Var)=diag(Var)*SSG
#     Gamma_star=t(MASS::mvrnorm(n=Multiple,Gamma,Sigma=Var)) #use the adaptive algorithm
#   }
#   
#   
#   
#   for (i in 1:Multiple){
#     log_p_Gamma_star[i]=log_p_gamma(Gamma_star[,i],Z,X,b_delta,G,Sigma_Gamma,Phi)
#     
#   }
#   
#   if(any(is.na(log_p_Gamma_star))){
#     print("overflow in Gamma density")
#   }
#   
#   #control overfloat, -max(log_p_Gamma_star);
#   max_log_p = median(log_p_Gamma_star)
#   p=exp(log_p_Gamma_star-max_log_p)/sum(exp(log_p_Gamma_star -max_log_p))
#   
#   #in case there is still overfloat;
#   p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
#   j=sample(c(1:Multiple),1,prob=p);	
#   
#   
#   Gamma_2star=t(rmvt(n=Multiple-1,Gamma_star[,j],sigma=Var,df=df))
#   Gamma_2star <-cbind(Gamma_2star,Gamma)
#   
#   for (i in 1:Multiple){
#     log_p_Gamma_2star[i]=log_p_gamma(Gamma_2star[,i],Z,X,b_delta,G,Sigma_Gamma,Phi)
#   }
#   
#   
#   #control overfloat
#   num=sum(exp(log_p_Gamma_star -max_log_p))
#   den=sum(exp(log_p_Gamma_2star -max_log_p))
#   
#   rho=min(1,num/den)
#   
#   #in case overfloat again
#   if(is.na(rho)) {rho=0.5};
#   
#   
#   accp=0;
#   u=runif(1)
#   if(u<rho){
#     Gamma=Gamma_star[,j]
#     accp=1;
#     
#   }
#   
#   
#   return(list(par=Gamma, accp=accp))
# }

Draw_b_delta_i = function(b_delta_i, i, G, X, gamma, sigma_b_delta_sq, Phi, var = 1, Multiple = 3,y)
{
  indx = G==i
  J = sum(indx)
  Phi_i = Phi[indx]
  y_i = y[indx]
  X_i = X[indx,]
  v_square = sigma_b_delta_sq/(1+J*sigma_b_delta_sq)
  mu = v_square * sum(y_i - X_i%*%gamma)
  b_delta_i_new = rnorm(1, mean = mu, sd = sqrt(v_square))
  accp = 1
  return(list(par=b_delta_i_new, accp=accp))
} 



log_p_b_delta = function(b_delta_i, i, G, X, gamma, sigma_b_delta_sq, Phi,y){
  indx = G==i
  Phi_i = Phi[indx]
  y_i = y[indx]
  if(ncol(X)==1){# intercept only
    X_i = X[indx]
    log_p = -0.5 * b_delta_i^2/sigma_b_delta_sq +
      sum(pnorm(X_i[indx_phi]%*%gamma + b_delta_i,log.p = T))+sum(log((1-pnorm(X_i[!indx_phi]%*%gamma + b_delta_i))))
  }else{
    X_i = X[indx,]
    indx_phi = Phi_i==1
    log_p = -0.5 * b_delta_i^2/sigma_b_delta_sq +
      sum(dnorm(y_i - X_i%*%gamma - b_delta_i, log=T))
      #sum(pnorm(X_i[indx_phi,]%*%gamma + b_delta_i,log.p = T))+sum(log((1-pnorm(X_i[!indx_phi,]%*%gamma + b_delta_i))))
      #
  }
  return(log_p)
}
# 
# 
# 
# Draw_b_delta_i = function(b_delta_i, i, G, X, gamma, sigma_b_delta_sq, Phi, var = 1, Multiple = 3,y)
#   
# {
#   
#   # log_p_b_delta_star=rep(0,Multiple)
#   # log_p_b_delta_2star=rep(0,Multiple)
#   # p=rep(0,Multiple)
#   # b_delta_star=rnorm(n=Multiple,b_delta_i,sd = sqrt(var)) #use the adaptive algorithm
#   # d_b_delta_star = dnorm(b_delta_star, mean = b_delta_i, sd = sqrt(var))
#   # for (i in 1:Multiple){
#   #   log_p_b_delta_star[i]=log_p_b_delta(b_delta_star[i], i, G, X, gamma, sigma_b_delta_sq, Phi)
#   # }
#   # 
#   # if(any(is.na(log_p_b_delta_star))){
#   #   print("overflow in Gamma density")
#   # }
#   # 
#   # p=exp(log_p_b_delta_star-log(d_b_delta_star))/sum(exp(log_p_b_delta_star -log(d_b_delta_star)))
#   # 
#   # #in case there is still overfloat;
#   # p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
#   # j=sample(c(1:Multiple),1,prob=p);
#   # 
#   # 
#   # b_delta_2star=rnorm(n=Multiple-1,b_delta_star[j],sd = sqrt(var))
#   # b_delta_2star <-c(b_delta_2star,b_delta_i)
#   # d_b_delta_2star = dnorm(b_delta_2star, mean = b_delta_star[j], sd = sqrt(var))
#   # for (i in 1:Multiple){
#   #   log_p_b_delta_2star[i]=log_p_b_delta(b_delta_2star[i], i, G, X, gamma, sigma_b_delta_sq, Phi)
#   # }
#   # 
#   # 
#   # #control overfloat
#   # num=sum(exp(log_p_b_delta_star -log(d_b_delta_star)))
#   # den=sum(exp(log_p_b_delta_2star -log(d_b_delta_2star)))
#   # 
#   # rho=min(1,num/den)
#   # 
#   # #in case overfloat again
#   # if(is.na(rho)) {rho=0.5};
#   # 
#   # 
#   # accp=0;
#   # u=runif(1)
#   # if(u<rho){
#   #   b_delta_i=b_delta_star[j]
#   #   accp=1;
#   # }
#   # return(list(par=b_delta_i, accp=accp))
# 
#   
#   b_delta_i_star=rnorm(n=1,b_delta_i,sqrt(var))
#   #browser()
#   log_p_b_delta_i_star = log_p_b_delta(b_delta_i_star, i, G, X, gamma,
#                                        sigma_b_delta_sq, Phi,y=y)
#   log_p_b_delta_i_old = log_p_b_delta(b_delta_i, i, G, X, gamma,
#                                       sigma_b_delta_sq, Phi,y=y)
#   rho = min(0, log_p_b_delta_i_star - log_p_b_delta_i_old)
#   accp = 0
#   u = log(runif(1, min = 0, max = 1))
#   if(is.na(rho)) {
#     rho=log(0.5)
#     print("b_delta overflow")
#   };
#   if(u < rho){
#     b_delta_i = b_delta_i_star
#     accp = 1
#   }
#   return(list(par=b_delta_i, accp=accp))
# } 


log_p_b_eta = function(b_eta_i, i, G, X_i, gamma, Eta_i,
                       Alpha, Sigma_b_eta){
  Ind=matrix(0,nrow=length(Eta_i),ncol=K)
  for (i in 1: length(Eta_i)){
    Ind[i,Eta_i[i]]=1;
  }
  if(ncol(X)==1){# intercept only
    X_i = X[G==i]
    log_p = -0.5 * t(b_eta_i) %*% solve(Sigma_b_eta) %*% b_eta_i + 
      sum(log(
        (exp(X_i * gamma + b_delta_i)/(1 + exp(X_i * gamma + b_delta_i))) *
          rowSums(
            Gden_i*exp(t(t(X_i %*% Alpha) + c(0,b_eta_i)))/colSums(exp(t(X_i %*% Alpha) + c(0,b_eta_i))) 
          )*as.numeric(Z_i>0) + 
          as.numeric(Z_i==0)/(1 + exp(X_i * gamma + b_delta_i))
      ))
  }else{
    log_p = -0.5 * t(b_eta_i) %*% solve(Sigma_b_eta) %*% b_eta_i + 
      sum(
        rowSums(Ind*(t(t(X_i %*% Alpha) + c(0,b_eta_i))- log(rowSums( exp(t(t(X_i%*%Alpha)+c(0,b_eta_i))) )) ))
      )
    #checking the calculation of rowSums(Gden_i*exp(t(t(X_i %*% Alpha) + c(0,b_eta_i)))/colSums(exp(t(X_i %*% Alpha) + c(0,b_eta_i))))
    # temp = rep(0, length(Z_i))
    # for(j in 1:length(Z_i)){
    #   X_ij = X_i[j,]
    #   Gden_ij = Gden_i[j,]
    #   temp[j] = sum(exp(t(X_ij)%*%Alpha + c(0,b_eta_i))/sum(exp(t(X_ij)%*%Alpha + c(0,b_eta_i))) * Gden_ij)
    # }
  }
  return(log_p)
}



Draw_b_eta_i = function(b_eta_i, i, G, X, gamma, Eta = Eta,Z,
                        Alpha, Sigma_b_eta, K = 5, iter,burnIn, individual = F,
                        cov = diag(0.5, nrow = K-1)){
  # X_i = X[G==i & Z > 0,]
  # Eta_i = Eta[G==i & Z>0]
  # if(length(Eta_i)==0){
  #   return(list(par=MASS::mvrnorm(n = 1, mu = rep(0, ncol(Alpha)-1),Sigma = Sigma_b_eta), accp=1))
  # }
  # b_eta_opt <-tryCatch({optim(b_eta_i, i=i, G=G, X=X_i, gamma=gamma, Eta_i = Eta_i,
  #                       Alpha=Alpha, Sigma_b_eta=Sigma_b_eta,
  #                       log_p_b_eta,method="Nelder-Mead",
  #                       hessian=TRUE,control=list(maxit=10,fnscale=-1))},error=function(e) NULL)
  # if(is.null(b_eta_opt) == T){
  #   return(list(par=rnorm(K-1,0,1), accp=0))
  # }
  # b_eta_opt_mean=b_eta_opt$par
  # sigma=solve(-b_eta_opt$hessian)
  # 
  # if(iter > burnIn & sample(c(0,1),1,prob=c(0.5, 0.5))==1){
  #   diag(sigma)=diag(sigma)*1
  # }
  # 
  # b_eta_i_star=tryCatch({MASS::mvrnorm(n = 1, mu = b_eta_opt_mean,Sigma = sigma)}, error=function(e) NULL)# sometimes the sigma is negative definite...
  # if(is.null(b_eta_i_star) == T){
  #   b_eta_i_star = MASS::mvrnorm(n = 1, mu = b_eta_opt_mean,Sigma = diag(1, nrow = K-1, ncol = K-1))
  # }
  # log_p_b_eta_i_star = log_p_b_eta(b_eta_i_star, i, G, X_i, gamma,Eta_i,
  #                                  Alpha, Sigma_b_eta)
  # log_p_b_eta_i_old = log_p_b_eta(b_eta_i, i, G, X_i, gamma,Eta_i,
  #                                 Alpha, Sigma_b_eta)
  # rho = min(0, log_p_b_eta_i_star - log_p_b_eta_i_old)
  # accp = 0
  # u = log(runif(1, min = 0, max = 1))
  # if(is.na(rho)) {rho=log(0.5)};
  # if(u < rho){
  #   b_eta_i = b_eta_i_star
  #   accp = 1
  # }
  # return(list(par=b_eta_i, accp=accp))
  
  #browser()
  X_i = X[G==i & Z > 0,]
  Eta_i = Eta[G==i & Z>0]
  if(length(Eta_i)==0){
    return(list(par=MASS::mvrnorm(n = 1, mu = rep(0, ncol(Alpha)-1),Sigma = Sigma_b_eta), accp=1))
  }
  # multiple = 2
  # log_p_Alpha_star=rep(0,multiple)
  # log_p_Alpha_2star=rep(0,multiple)
  # # b_eta_opt <-tryCatch({optim(b_eta_i, i=i, G=G, X=X_i, gamma=gamma, Eta_i = Eta_i,
  # #                             Alpha=Alpha, Sigma_b_eta=Sigma_b_eta,
  # #                             log_p_b_eta,method="Nelder-Mead",
  # #                             hessian=TRUE,control=list(maxit=10,fnscale=-1))},error=function(e) NULL)
  # # if(is.null(b_eta_opt) == T){
  # #   return(list(par=rnorm(K-1,0,1), accp=0))
  # # }
  # # b_eta_opt_mean=b_eta_opt$par
  # # sigma=solve(-b_eta_opt$hessian)
  # # cov = sigma
  # # b_eta_i = b_eta_opt_mean
  # b_eta_i_star = rmvt(n=multiple,delta = b_eta_i,sigma=cov,df=5)
  # #b_eta_i_star = mvtnorm::rmvnorm(n = multiple, mean = b_eta_i,sigma = cov)
  # 
  # for (j in 1:multiple){
  #   log_p_Alpha_star[j]=log_p_b_eta(b_eta_i_star[j,], i, G, X_i, gamma,Eta_i, Alpha, Sigma_b_eta)
  # }
  # d_p_Alpha_star = mvtnorm::dmvnorm(b_eta_i_star, mean = b_eta_i, sigma = cov, log = T)
  # if(any(is.na(log_p_Alpha_star))){
  #   print("overflow in b_eta density")
  # }
  # # browser()
  # max_p_alpha = median(log_p_Alpha_star)
  # # log_p_Alpha_star[which(abs(log_p_Alpha_star) == max(abs(log_p_Alpha_star)),log_p_Alpha_star)]
  # #control overfloat, -max(log_p_Alpha_star);
  # p=exp(log_p_Alpha_star-d_p_Alpha_star-max_p_alpha)/sum(exp(log_p_Alpha_star - d_p_Alpha_star-max_p_alpha))
  # #in case there is still overfloat;
  # p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
  # v=sample(c(1:multiple),1,prob=p);
  # 
  # #b_eta_i_2star= mvtnorm::rmvnorm(n = multiple-1, mean = b_eta_i_star[v,],sigma = cov)
  # b_eta_i_2star = rmvt(n=multiple-1,delta = b_eta_i_star[v,],sigma=cov,df=5)
  # b_eta_i_2star <-rbind(b_eta_i_2star,b_eta_i)
  # 
  # for (j in 1:multiple){
  #   log_p_Alpha_2star[j]=log_p_b_eta(b_eta_i_2star[j,], i, G, X_i, gamma,Eta_i, Alpha, Sigma_b_eta)
  # }
  # d_p_Alpha_2star = mvtnorm::dmvnorm(b_eta_i_2star, mean = b_eta_i_star[v,], sigma = cov, log = T)
  # if(any(is.na(log_p_Alpha_2star))){
  #   print("overflow in b_eta density")
  # }
  # 
  # #control overfloat
  # num=sum(exp(log_p_Alpha_star -d_p_Alpha_star-max_p_alpha))
  # den=sum(exp(log_p_Alpha_2star -d_p_Alpha_2star-max_p_alpha))
  # 
  # rho=min(1,num/den)
  # 
  # #in case overfloat again
  # if(is.na(rho)) {
  #   print("b_eta rho is na")
  #   rho=0.05};
  # 
  # accp=0;
  # u=runif(1)
  # if(u<rho){
  #   b_eta_i=b_eta_i_star[v,]
  #   accp=1;
  # }
  # return(list(par=b_eta_i, accp=accp))
  
  
  # X_i = X[G==i & Z > 0,]
  # Eta_i = Eta[G==i & Z>0]
  # if(length(Eta_i)==0){
  #   return(list(par=MASS::mvrnorm(n = 1, mu = rep(0, ncol(Alpha)-1),Sigma = Sigma_b_eta), accp=1))
  # }
  b_eta_i_star = MASS::mvrnorm(n = 1, mu = b_eta_i,Sigma = cov)


  log_p_b_eta_i_star = log_p_b_eta(b_eta_i_star, i, G, X_i, gamma,Eta_i,
                                   Alpha, Sigma_b_eta)
  log_p_b_eta_i_old = log_p_b_eta(b_eta_i, i, G, X_i, gamma,Eta_i,
                                  Alpha, Sigma_b_eta)
  rho = min(0, log_p_b_eta_i_star - log_p_b_eta_i_old)
  accp = 0
  u = log(runif(1, min = 0, max = 1))
  if(is.na(rho)) {rho=log(0.5)};
  if(u < rho){
    b_eta_i = b_eta_i_star
    accp = 1
  }
  return(list(par=b_eta_i, accp=accp))
  
} 
