
Bayesian_Semi=function(Z,X,G=NULL,K=5,nIter=1100,
                       burnIn=100,thin=5,return_Bden=T,useCov=T,
                       diag_Sigma_b_eta = F, true_para = NULL,
                       byStep=0.01, same_b_eta_variance = T, 
                       s_b_eta = 1.1, plotWhileRun = F,
                       check_para = NULL,zero_inflate = T, ref_b_spline = 1,
                       auto_check_convergence = F, max_Z = NULL, same_var = F, 
                       c = 0.01, d = 0.01, alpha_prior_var = 10, 
                       gamma_prior_var=2.5) 
{
  #Z: Vector of semi-continuous data values
  #X: A matrix of covariates, with rows as observations same length as vector Z
  #G: group index for the observations
  #K: Number of b-spline density basis, greater than 5
  #nIter: Total number of MCMC iterations
  #burnIn: number of burn in iterations
  #thin: MCMC sampling thin number
  #return_Bden: If True, returns the original B-spline density basis
  #useCov: equal to TRUE if use covariates (not only intercept)
  #diag_Sigma_b_eta=T, then the Sigma_b_eta is digonal with variance follows a inverse gamma prior
  #byStep=0.01, the step size for constructing the b densities, if Zs are integer, we could do stepsize  = 1
  #same_b_eta_variance = T, constraint all the b_eta on different b densities to have the same variance
  #s_b_eta：scalar for adpative MH algorithm
  #plotWhileRun=F：Plot out the convergence figures while running
  #check_para = NULL, then no computing rate for credible interval covering the true value. Or else
  #             this should be a list contain true value for Gamma (vector) and Alpha (matrix) with rows corresponding to covariates and columns to b-splines,
  #             Output would be a list contain two vector for Gamma (mean and whether CI covers the true value, 1 for yes), two matrix (mean for alpha and whether CI covers the true value)
  #zero_inflate = T, then we would have zero inflated part. If false, then only the b-spline part exists.
  #ref_b_spline = 1:this is the index for reference b-spline of Alpha coefficient
  #auto_check_convergence: if TRUE, auto checks convergence, if converged could end MCMC before the specified nIter
  #max_Z: the true maximum of Z, only validate in simulation, for B_density
  #same_var: if it is T then the variance for global intercept and local random intercept are the same
  #c, d: the hyperparameters for both sigma_b_eta and sigma_b_delta prior
  #alpha_prior_var: the sigma_b_eta prior variance, a scaler so that we assume every alpha_k are the same prior
  #gamma_prior_var: prior variance for gamma 
  
  
  #Load packages/functions
  library(tmvtnorm)
  library(mnormt)
  library(magic)
  library(arm)
  library(pscl)
  library(splines)
  library(plyr)
  library(doParallel)
  library(coda)
  
  # Creating a folder for plot pdf
  if(plotWhileRun == T){
    dir.create("temp_results", showWarnings = FALSE)
  }
  
  # Check if there is groups of subjects, change the grouping number to 1:length(unique(G))
  if(is.null(G)){
    n_sub = length(Z)
    G = 1:n_sub
    no_G = T
  }else{
    G_old = G
    G = as.numeric(as.factor(G))
    n_sub = length(unique(G)) # number of groups in the dataset
    no_G = F
  }
  
  # Check if there is missing cases
  if(sum(complete.cases(X))!=nrow(X) | sum(is.na(Z))>0){
    stop("missing cases in data");return;
  }
  
  # Check if there is at least 5 b-spline basis
  if(K<4 & K!= 3){stop("K needs to be >=5 or equal to 3");return; }
  
  N=length(Z)# total number of observations
  
  # see if subjects have the same number of observations
  same_length_G = F
  if(var(plyr::count(G)$freq)==0){
    same_length_G = T
  }
  # Whether if there is an intercept
  if (useCov==T){
    if(length(which(X[,1]==1))!= length(Z)) {X=cbind(Intcpt=1,X)}
    M=dim(X)[2]
  }else{
    X=matrix(1,nrow=N,ncol=1)
    
    M=1;
  }
  # Obtain the variable names if there is any
  if(!is.null(colnames(X))){
    col_names  = colnames(X)
  }else{
    col_names = NULL
  }
  # hyperparameters	
  Sigma_Gamma=matrix(0,nrow=M,ncol=M)
  diag(Sigma_Gamma)= gamma_prior_var; #Prior: P(Gamma) ~ N(0,Sigma_Gamma)
  
  
  
  # Parameter arrays	
  
  ALPHA_array=list()  #  M * K Matrix (here M is the M+1 in the paper)
  Tau_sq_array=list() # M vector 
  a_array=list() # M vector
  GAMMA_array=list() # M+1 vector
  sigma_delta_sq_array=list()
  Sigma_eta_array=list()
  Accp_Rate_array_a=list() 	#Alpha draw accept rate
  Accp_Rate_array_g=list()	#Gamma draw accept rate
  Accp_Rate_array_b_delta=list()
  Accp_Rate_array_b_eta=list()
  b_delta_array = list()
  b_eta_array = list()
  effect_mean_array = list()
  effect_median_array =list()
  effect_0.25_array = list()
  effect_0.75_array =list()
  effect_dist_array = list()
  array_ind=0
  
  
  # Initialize global parameters 	
  pi0=0.5 
  gamma0=log((1-pi0)/pi0) 
  Gamma=array(0,dim=c(M,1)); Gamma[1]=gamma0
  Tau_sq=array(1,dim=c(M,1)) ; 
  if(no_G==F){
    sigma_b_delta_sq = 1
    Sigma_b_eta = diag(1, nrow = K-1)
  }else{
    sigma_b_delta_sq = 0
    Sigma_b_eta = diag(0, nrow = K-1)
  }
  b_delta = rnorm(n = n_sub, mean = 0, sd = sqrt(sigma_b_delta_sq))
  b_eta = MASS::mvrnorm(n = n_sub, mu = rep(0, K-1), Sigma = Sigma_b_eta) # n_sub by K-1 matrix
  
  Gamma_mean=array(0,dim=c(M,1));
  Tau_sq_mean=array(0,dim=c(M,1));
  a_mean=array(0,dim=c(M,1));
  
  Phi=1-as.numeric(abs(Z)<= sort(abs(Z))[round(pi0*N)]);
  Phi_mean=0*Phi	
  PHI_match_rate=NULL
  
  Alpha=matrix(0,nrow=M,ncol=K);
  for (m in 1:M){
    Alpha[m,2:3]=mvrnorm(n = 1, c(0,0), diag(2))
    if(K >= 4){
      for (i in 4:K){
        Alpha[m,i]=rnorm(1,mean=0,sd=diag(2))	
      }
    }
  }
  Alpha_mean=matrix(0,nrow=M,ncol=K);
  Alpha_accp=rep(0,M); #draw alpha[m,] from m=1 to M;
  y = rep(0,N) # Latent variable for probit
  
  ## Adaptive MH alogrithm
  C_t = list() # For b_eta list contain the covariance matrix n_sub elements with (K-1)*(K-1)
  b_eta_mean = b_eta # The previous iter mean, n_sub*(K-1)
  
  # Getting the b-spline density basis
  #source("../functions/Bdensity.R")
  Gden<-matrix(0,nrow=N,ncol=K);
  Gden<-Bdensity(byStep,mu = 0,K,N,Z,Gden, ref = ref_b_spline, max_Z = max_Z) # Get the density for Z
  Gden_raw <- Bdensity_raw(byStep,mu = 0,K,N,Z,ref = ref_b_spline, max_Z = max_Z) # return the original b-spline density
  
  # Check if model is non-zero inflated or there is no random effect
  if(zero_inflate==F){
    Gamma=array(0,dim=c(M,1))
    Gamma[1]=100 # the intercept is very large, so that the non-null probability would be 1 when drawing b_eta
    b_delta = rep(0, n_sub)
    Phi = rep(1, N)
  }
  if(no_G == T){
    b_delta = rep(0, n_sub)
    b_eta = matrix(0L, nrow = dim(b_eta)[1], ncol = dim(b_eta)[2])  # n_sub by K-1 matrix
  }
  
  #beginning of the iterations
  for(iter in 1:nIter){
    print(iter)
    
    #Getting the non-zero Zs
    Z1=abs(Z[Phi==1])
    #Getting the covariates for non-zero Zs
    if(useCov==T){
      X1=X[Phi==1,]
      G1=G[Phi==1]
      Gden1 = Gden[Phi==1,]
    }else{
      X1=X[Phi==1]
      G1=G[Phi==1]
      Gden1 = Gden[Phi==1,]
    }
    
    # Index for non-zero observations
    nn=(1:N)[Phi==1]
    
    
    #get B-spline densities for Z1;
    gden<-matrix(0,nrow=length(Z1),ncol=K);
    gden<-Gden[nn,];
    
    
    # Draw Eta: local indicator
    Eta<-rep(0,N);
    XA=X1%*%Alpha; #length(nn)*M times M*K
    j=1;
    for (i in nn){
      g_indx = G[i] # the group number for ith observation
      p=exp(XA[j,] + c(0,b_eta[g_indx,]))*gden[j,];
      p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
      Eta[i]=try(sample(c(1:K),size=1,replace=T,prob=p),silent=TRUE);
      if(length(grep("Error", Eta[i], fixed=TRUE))>0){
        Eta[i]=sample(c(1:K),size=1,replace=T,prob=rep(1/K,K))
      }
      j=j+1;
    }
    Eta = as.numeric(Eta)
    
    
    # Draw Gamma:global coefficients
    objg=Draw_Gamma_log_M(y=y, X=X,G=G,SG = Sigma_Gamma, Phi =Phi,sigma_b_delta_sq=sigma_b_delta_sq, same_length_G = same_length_G)
    Gamma=objg$par;
    
    #browser()
    # Draw latent variable y
    y = c()
    i = 1
    if(same_length_G){
      if(no_G==F){# if there is grouping
        H = Matrix(solve(objg$V[[1]]), sparse = T)# Sparse matrix for faster computation
        for(g in unique(G)){
          indx = G==g
          Phi_i = Phi[indx]
          X_i = X[indx,]
          lb = ifelse(Phi_i==1, 0, -Inf)
          ub = ifelse(Phi_i==1, Inf, 0)
          y_temp = tmvtnorm::rtmvnorm.sparseMatrix(1, mean = c(X_i%*%Gamma), H = H, lower = lb, upper = ub,burn.in.samples = 100)
          #y_temp = tmvtnorm::rtmvnorm(1, mean = c(X_i%*%Gamma),sigma = objg$V[[i]], lower = lb, upper = ub,alogrithm="gibbs",burn.in.samples = 300)
          i = i + 1
          y = c(y, y_temp)
        }
      }else{# no random effect case
        for(g in unique(G)){
          indx = G==g
          Phi_i = Phi[indx]
          X_i = X[indx,]
          lb = ifelse(Phi_i==1, 0, -Inf)
          ub = ifelse(Phi_i==1, Inf, 0)
          #y_temp = tmvtnorm::rtmvnorm(1, mean = c(X_i%*%Gamma),sigma = objg$V[[i]], lower = lb, upper = ub, algorithm = "gibbs",burn.in.samples = 100)
          y_temp = truncnorm::rtruncnorm(1,a = lb, b= ub, mean = c(X_i%*%Gamma),sd = objg$V[[i]])
          #y_temp = tmvtnorm::rtmvnorm(1, mean = c(X_i%*%Gamma),sigma = objg$V[[i]], lower = lb, upper = ub)
          i = i + 1
          y = c(y, y_temp)
        }
      }
    }else{#no the same number of observation for each subject
      if(no_G==F){
        for(g in unique(G)){
          indx = G==g
          Phi_i = Phi[indx]
          X_i = X[indx,]
          lb = ifelse(Phi_i==1, 0, -Inf)
          ub = ifelse(Phi_i==1, Inf, 0)
          #y_temp = tmvtnorm::rtmvnorm.sparseMatrix(1, mean = c(X_i%*%Gamma), H = Matrix(solve(objg$V[[i]]), sparse = T), lower = lb, upper = ub,burn.in.samples = 100)
          y_temp = tmvtnorm::rtmvnorm(1, mean = c(X_i%*%Gamma),sigma = objg$V[[i]], lower = lb, upper = ub, algorithm = "gibbs")
          i = i + 1
          y = c(y, y_temp)
        }
      }else{
        for(g in unique(G)){
          indx = G==g
          Phi_i = Phi[indx]
          X_i = X[indx,]
          lb = ifelse(Phi_i==1, 0, -Inf)
          ub = ifelse(Phi_i==1, Inf, 0)
          #y_temp = tmvtnorm::rtmvnorm.sparseMatrix(1, mean = c(X_i%*%Gamma), H = Matrix(solve(objg$V[[i]]), sparse = T), lower = lb, upper = ub,burn.in.samples = 100)
          y_temp = truncnorm::rtruncnorm(1,a = lb, b= ub, mean = c(X_i%*%Gamma),sd = objg$V[[i]])
          i = i + 1 
          y = c(y, y_temp)
        }
      }
      
    }
    
    
    
    
    # Draw b_delta
    if(zero_inflate == T & no_G == F){#zero inflated and random effect
      if(is.null(true_para$b_delta) == F){#For model diagnostics, plug in true parameter
        b_delta = true_para$b_delta
        b_delta_accp = 1
      }else{
        b_delta_accp = rep(0, n_sub)
        for(i in 1:n_sub){
          b_delta_i = b_delta[i]
          obj_b_delta = Draw_b_delta_i(b_delta_i, i, G, X, Gamma, sigma_b_delta_sq, Phi, var = NULL,y=y)
          b_delta[i] = obj_b_delta$par
          b_delta_accp[i] = obj_b_delta$accp#here the acceptance is always 1
        }
        b_delta_accp = mean(b_delta_accp)
      }
    }
    
    # Draw sigma_b_delta_sq 
    if(zero_inflate==T & no_G == F){
      if(is.null(true_para$sigma_b_delta_sq)){
        sigma_b_delta_sq = rigamma(1, alpha = c + n_sub/2, beta = d + 0.5*sum(b_delta^2))
      }else{#For model diagnoistics, plug in true parameter
        sigma_b_delta_sq = true_para$sigma_b_delta_sq
      }
    }
    
    # Draw ALPHA (PG)
    if(is.null(true_para$Alpha)==F){#For model diagnostics, plug in true parameter
      Alpha = true_para$Alpha
      Alpha_accp[1:M]=1
    }else{
      for (k in 2:K){
        if(no_G==F){
          obj=Draw_Alpha_PG(Alpha[,k],Alpha[,-k],X1,b_eta,G1,Eta[!Eta ==0],k, Sigma_alpha = diag(rep(alpha_prior_var, M-1)))	#no intercept in the draw		
          Alpha[2:M,k]=obj$par;	
        }else{
          obj=Draw_Alpha_PG_full(Alpha[,k],Alpha[,-k],X1,b_eta,G1,Eta[!Eta ==0],k, Sigma_alpha = diag(c(10,rep(alpha_prior_var, M-1))))	#intercept in the draw
          Alpha[1:M,k]=obj$par;
        }
      }
    }
    Alpha_accp[1:M]=1;
    
    # Draw b_eta
    if(no_G==F){
      if(is.null(true_para$b_eta) == F){#For model diagnostics, plug in true parameter
        b_eta = true_para$b_eta
        b_eta_accp = 1
      }else{
        b_eta_accp = rep(0, n_sub)
        for(i in 1:n_sub){
          b_eta_i = b_eta[i,]
          #calculate the MH covariance
          if(iter == 1){
            C_t[[i]] = s_b_eta * Sigma_b_eta + s_b_eta*0.0001*diag(1, nrow = K-1)
          }else{
            t = iter - 1
            b_eta_mean_old_i = b_eta_mean[i,] #old mean b_eta
            b_eta_mean[i,] = (b_eta_mean_old_i * (iter-1) + b_eta_i)/iter # new mean
            C_t[[i]] = (t-1)*C_t[[i]]/t + s_b_eta*(t*b_eta_mean_old_i%*%t(b_eta_mean_old_i)-
                                                     (t+1)*b_eta_mean[i,]%*%t(b_eta_mean[i,])+
                                                     b_eta_i%*%t(b_eta_i)+
                                                     0.0001*diag(1, nrow = K-1))/t
          }
          #Adapative MH algorithm
          obj_b_eta = Draw_b_eta_i(b_eta_i, i, G, X, Gamma, Eta = Eta, Z,
                                   Alpha, Sigma_b_eta, K = K, iter = iter,burnIn = burnIn, 
                                   individual = F, cov = C_t[[i]])
          b_eta[i,] = obj_b_eta$par
          b_eta_accp[i] = obj_b_eta$accp
          
        }
        b_eta_accp = mean(b_eta_accp)
      }
    }
    
    
    # Draw Sigma_b_eta
    if(no_G == F){
      if(is.null(true_para$Sigma_b_eta)==F){#For model diagnostics, plug in true parameter
        Sigma_b_eta = true_para$Sigma_b_eta
      }else if(same_var ==T){# If the variance is same as the global random effect
        Sigma_b_eta = diag(sigma_b_delta_sq, nrow = K-1)
      }else{
        if(diag_Sigma_b_eta == T & same_b_eta_variance == F){#Diagonal covariance matrix with entries different
          G1_index = unique(G1)
          N1 = length(G1_index)
          for(k in 1:(K-1)){
            b_eta_k = b_eta[G1_index,k]
            Sigma_b_eta[k,k] = rigamma(1, alpha = c + N1/2, beta = d + 0.5*sum((b_eta_k)^2))
          }
        }else if(same_b_eta_variance == T){# Diagonal covariance matrix with entries same
          G1_index = unique(G1)
          N1 = length(G1_index)
          b_eta_G1 = c(b_eta[G1_index,])
          Sigma_b_eta = diag(rigamma(1, alpha = c + (K-1)*N1/2, beta = d + 0.5*sum(b_eta_G1^2)), nrow = K-1)
          # Sigma_b_eta = diag(rigamma(1, alpha = vu/2 + (K-1)*N1/2, beta = vu/a + 0.5*sum(b_eta_G1^2)), nrow = K-1)
          # a = rigamma(1, alpha = 1/2, beta = 1/A^2)
        }else{# Unstructured covariance matrix
          G1_indx = unique(G1)
          Sigma_b_eta = MCMCpack::riwish(length(G1_indx) + K-1, diag(1, nrow = K-1) +
                                           t(b_eta[G1_indx, ]) %*% b_eta[G1_indx, ]
          )
          # Sigma_b_eta = LaplacesDemon::rinvwishart(nu = length(G1_indx)+K-1, 
          #                                          S = diag(1, nrow = K-1) + 
          #                                            t(b_eta[G1_indx, ]) %*% b_eta[G1_indx, ])
        }
      }
    }
    
    #Draw \beta_0
    if(no_G==F){
      var_beta0 = Sigma_b_eta[1,1]/n_sub
      for(k in 2:K){
        mean_beta0 = Alpha[1,k]+sum(b_eta[,k-1])/n_sub# mean of beta0 of all the runs
        beta0_old = Alpha[1,k]#old beta0
        Alpha[1,k]=rnorm(1, mean = mean_beta0, sd = sqrt(var_beta0))# sample beta0
        b_eta[,k-1] = b_eta[,k-1]+beta0_old-Alpha[1,k]# adjust for the random intercept
      }
    }
    
    
    if(zero_inflate == T){# zero inflation then add delta indicator
      if(is.null(true_para$delta)){
        Phi = ifelse(Z>0,1,0)
        P_phi = rep(1,N)
      }else{
        Phi = true_para$delta
        P_phi = rep(1,N)
      }
    }
    
    
    
    # Save posterior sample in array once >burn in for each thin
    if(iter%%thin==0 & iter>=burnIn){
      array_ind=array_ind+1
      dist_effect = unlist(cal_quantity(NULL, x = rep(1, length(Gamma)), Alpha, Gamma, 
                                        Gden_raw$phi.mat, Gden_raw$grid, byStep))
      effect = coeff_effect(NULL, Alpha, Gamma, 
                            Gden_raw$phi.mat, Gden_raw$grid, byStep)
      effect_dist_array[[array_ind]] = dist_effect
      effect_mean_array[[array_ind]] = effect$effect_mean
      effect_median_array[[array_ind]] = effect$effect_median
      effect_0.25_array[[array_ind]] = effect$effect_0.25
      effect_0.75_array[[array_ind]] = effect$effect_0.75
      ALPHA_array[[array_ind]]=Alpha
      if(zero_inflate==T){
        GAMMA_array[[array_ind]]=Gamma
        if(no_G==F){
          b_delta_array[[array_ind]]=b_delta
          sigma_delta_sq_array[[array_ind]]=sigma_b_delta_sq
          Accp_Rate_array_b_delta[[array_ind]]=b_delta_accp
        }
        Accp_Rate_array_g[[array_ind]]=objg$accp;
      }
      Tau_sq_array[[array_ind]]=Tau_sq
      if(no_G==F){
        b_eta_array[[array_ind]]=b_eta
        Sigma_eta_array[[array_ind]]=Sigma_b_eta
        Accp_Rate_array_b_eta[[array_ind]]=b_eta_accp
      }
      Accp_Rate_array_a[[array_ind]]=Alpha_accp;
      for (m in 1:M){
        print(paste("Mean of Alpha m =",m));
        Alpha_mean[m,]=((array_ind-1)*Alpha_mean[m,]+ALPHA_array[[array_ind]][m,])/array_ind
        print(Alpha_mean[m,])
        print(paste("Average Multiple-try MH Accept Rate for Alpha m =",m,":",mean(sapply(Accp_Rate_array_a, '[', m))));
      }
      if(zero_inflate == T){
        print("Gamma mean:");
        print(Gamma_mean)
        print(paste("Multiple-try MH Accept Rate for Gamma (mean):",mean(as.numeric(Accp_Rate_array_g))));
      }
      
      if(no_G==F){
        print("random subject mean b_eta");
        print(b_eta_mean[sample(1:n_sub,2),])
        print(paste("Multiple-try MH Accept Rate for b_eta (mean):",mean(as.numeric(Accp_Rate_array_b_eta))));
      }
    }
    
    # Plotting code
    if(iter%%(thin*50)==0 & iter%%thin==0 & iter>=burnIn & plotWhileRun == T){
      if(plotWhileRun == T){
        pdf(file = "temp_results/Convergence plots.pdf")
        if(zero_inflate==T){
          mcmc_Gamma = mcmc(t(sapply(GAMMA_array, cbind)))
          plot(mcmc_Gamma, ylab = "Gamma")
        }
        
        for(i in 1:M){
          mcmc_Alpha = mcmc(t(sapply(ALPHA_array,function(x)x[i,2:K])))
          plot(mcmc_Alpha, ylab = paste("Alpha", i,"covariate",sep = ""))
        }
        if(zero_inflate==T & no_G==F){
          mcmc_sigma_delta= mcmc(unlist(sigma_delta_sq_array))
          plot(mcmc_sigma_delta, ylab = "sigma_delta_sq")
        }
        if(no_G==F){
          mcmc_Sigma_eta= mcmc(sapply(Sigma_eta_array, function(x)x[1,1]))
          plot(mcmc_Sigma_eta, ylab = "Sigma_eta_sq")
        }
        dev.off()
        
        results_tmp=list()
        results_tmp[["alpha_array"]]=ALPHA_array
        if(zero_inflate==T){
          results_tmp[["gamma_array"]]=GAMMA_array
          results_tmp[["accp_rate_array_g"]]=Accp_Rate_array_g
          if(no_G==F){
            results_tmp[["b_delta"]]=b_delta_array
            results_tmp[["sigma_b_delta_sq"]]=sigma_delta_sq_array
            results_tmp[["accp_rate_b_delta"]]=Accp_Rate_array_b_delta
          }
        }
        results_tmp[["tau_sq_array"]]=Tau_sq_array
        if(no_G==F){
          results_tmp[["b_eta"]]=b_eta_array
          results_tmp[["grouping_order"]]= sort(unique(G_old))
          results_tmp[["Sigma_b_eta"]]=Sigma_eta_array
          results_tmp[["accp_rate_b_eta"]]=Accp_Rate_array_b_eta
        }
        results_tmp[["accp_rate_array_a"]]=Accp_Rate_array_a
        results_tmp[["a_array"]] = a_array
        save(file="temp_results/results_tmp.R",results_tmp,X,Z,N,nIter,burnIn,thin,mu,K,useCov)
      }
    }
    
    # Check convergence automatically
    if(auto_check_convergence == T & iter%%(thin*100)==0 & iter>=(burnIn+100)){
      test_temp = c()
      if(zero_inflate==T){
        mcmc_Gamma = mcmc(t(sapply(tail(GAMMA_array, 1000), cbind)))
        test_Gamma = geweke.diag(mcmc_Gamma)$z
        test_temp = c(test_temp, 0.05 < pnorm(abs(test_Gamma), lower.tail = F))
      }
      for(i in 1:M){
        mcmc_Alpha = mcmc(t(sapply(tail(ALPHA_array,1000),function(x)x[i,2:K])))
        test_Alpha = geweke.diag(mcmc_Alpha)$z
        test_temp = c(test_temp, 0.05 < pnorm(abs(test_Alpha), lower.tail = F))
      }
      if(all(test_temp)){
        break
      }
    }
    
  }	
  
  #Calculate SD;
  for (m in 1:M){#Alpah
    print(paste("SD of Alpha m = ",m));
    for(k in 1:K){
      print(sd(sapply(ALPHA_array,'[',m,k)))
    }
  }
  if(zero_inflate == T){#Gamma
    print("Gamma SD:");
    for(i in 1:dim(X)[2]){
      print(sd(as.numeric(unlist(lapply(GAMMA_array, function(x) x[i])))))
    }
  }
  
  
  #return the posterior sampling if we don't have the true parameters to check for coverage and bias
  if(is.null(check_para)){
    results=list()
    results[["alpha_array"]]=ALPHA_array
    if(zero_inflate==T){
      results[["gamma_array"]]=GAMMA_array
      results[["b_delta"]]=b_delta_array
      results[["sigma_b_delta_sq"]]=sigma_delta_sq_array
      results[["accp_rate_b_delta"]]=Accp_Rate_array_b_delta
      results[["accp_rate_array_g"]]=Accp_Rate_array_g
      results[["gamma_mean"]]=Gamma_mean
      results[["p_phi"]]=P_phi
      results[["phi_mean"]]=Phi_mean
    }
    results[["tau_sq_array"]]=Tau_sq_array
    results[["b_eta"]]=b_eta_array
    results[["Sigma_b_eta"]]=Sigma_eta_array
    results[["accp_rate_b_eta"]]=Accp_Rate_array_b_eta
    if(no_G==F){
      results[["grouping_order"]]= sort(unique(G_old))
    }
    results[["alpha_mean"]]=Alpha_mean
    results[["accp_rate_array_a"]]=Accp_Rate_array_a
    results[["phi"]]=Phi #last iteration Phi value
    if (return_Bden ==T){
      results[["Bden"]]=Gden
    }
    results[["Bden_raw"]]=Gden_raw$phi.mat
    results[["grid"]]=Gden_raw$grid
    results[["eta"]] = Eta #last iteration Eta value
    results[["zero_inflate"]]=zero_inflate
    results[["grouped"]]=!no_G
    results[["col_names"]]=col_names
    results[["effect_dist"]] = effect_dist_array
    results[["effect_mean"]]=effect_mean_array 
    results[["effect_median"]]=effect_median_array
    results[["effect_0.25"]]=effect_0.25_array 
    results[["effect_0.75"]]=effect_0.75_array
    return(results)
  }else{#If there is true parameter vector, return coverage T/F and mean of posterior
    ##Alpha 
    if(!is.null(dim(check_para$Alpha))){
      AMN_lower = apply(simplify2array(ALPHA_array), 1:2, quantile, probs = 0.025)
      AMN_upper = apply(simplify2array(ALPHA_array), 1:2, quantile, probs = 0.975)
      cover_matrix_alpha = AMN_upper >= check_para$Alpha & check_para$Alpha >= AMN_lower
      AMN<-matrix(0,nrow=K,ncol=M);
      for (m in 1:M){		
        for(k in 1:K){   		 
          AMN[k,m]=mean(sapply(ALPHA_array,'[',m,k))
        }
      }
      ##If zero-inflated then calculate Gamma and sigma_delta mean and coverage
      if(zero_inflate==T){
        GMN_upper<-NULL
        for(i in 1:M){
          GMN_upper<-rbind(GMN_upper,quantile(as.numeric(unlist(lapply(GAMMA_array, function(x) x[i]))), probs = 0.975))
        }
        GMN_lower<-NULL
        for(i in 1:M){
          GMN_lower<-rbind(GMN_lower,quantile(as.numeric(unlist(lapply(GAMMA_array, function(x) x[i]))), probs = 0.025))
        }
        cover_vector_gamma = c(GMN_upper >= check_para$Gamma & check_para$Gamma >= GMN_lower)
        GMN<-NULL
        for(i in 1:M){
          GMN<-rbind(GMN,mean(as.numeric(unlist(lapply(GAMMA_array, function(x) x[i])))))
        }
        ## sigma_delta_sq
        sigma_delta_sq_mean = mean(unlist(sigma_delta_sq_array))
        sigma_delta_sq_upper = quantile(unlist(sigma_delta_sq_array), probs = 0.975)
        sigma_delta_sq_lower = quantile(unlist(sigma_delta_sq_array), probs = 0.025)
        if(no_G==F){
          cover_sigma_delta_sq = sigma_delta_sq_upper >= check_para$sigma_b_delta^2 &  check_para$sigma_b_delta^2 >= sigma_delta_sq_lower
        }else{
          cover_sigma_delta_sq = 1
        }
      }
      # Sigma_eta
      if(no_G==F){
        Sigma_eta_lower = apply(simplify2array(Sigma_eta_array), 1:2, quantile, probs = 0.025)
        Sigma_eta_upper = apply(simplify2array(Sigma_eta_array), 1:2, quantile, probs = 0.975)
        cover_matrix_Sigma_eta = Sigma_eta_upper >= check_para$Sigma_b_eta^2 & check_para$Sigma_b_eta^2 >= Sigma_eta_lower
        Sigma_eta_mean = apply(simplify2array(Sigma_eta_array), 1:2, mean)
      }else{
        Sigma_eta_lower = NA
        Sigma_eta_upper = NA
        cover_matrix_Sigma_eta = diag(1, K-1)
        Sigma_eta_mean = diag(0, K-1)
      }
    }
    
    ## Distribution effects mean 
    effect_mean_upper<-NULL
    for(i in 1:M){
      effect_mean_upper<-rbind(effect_mean_upper,quantile(as.numeric(unlist(lapply(effect_mean_array, function(x) x[i]))), probs = 0.975))
    }
    effect_mean_lower<-NULL
    for(i in 1:M){
      effect_mean_lower<-rbind(effect_mean_lower,quantile(as.numeric(unlist(lapply(effect_mean_array, function(x) x[i]))), probs = 0.025))
    }
    cover_vector_mean = c(effect_mean_upper >= check_para$effect_mean & check_para$effect_mean >= effect_mean_lower)
    effect_mean<-NULL
    for(i in 1:M){
      effect_mean<-rbind(effect_mean,mean(as.numeric(unlist(lapply(effect_mean_array, function(x) x[i])))))
    }
    ## Distribution effects median 
    effect_median_upper<-NULL
    for(i in 1:M){
      effect_median_upper<-rbind(effect_median_upper,quantile(as.numeric(unlist(lapply(effect_median_array, function(x) x[i]))), probs = 0.975))
    }
    effect_median_lower<-NULL
    for(i in 1:M){
      effect_median_lower<-rbind(effect_median_lower,quantile(as.numeric(unlist(lapply(effect_median_array, function(x) x[i]))), probs = 0.025))
    }
    cover_vector_median = c(effect_median_upper >= check_para$effect_median & check_para$effect_median >= effect_median_lower)
    effect_median<-NULL
    for(i in 1:M){
      effect_median<-rbind(effect_median,mean(as.numeric(unlist(lapply(effect_median_array, function(x) x[i])))))
    }
    ## Distribution effects 0.25
    effect_0.25_upper<-NULL
    for(i in 1:M){
      effect_0.25_upper<-rbind(effect_0.25_upper,quantile(as.numeric(unlist(lapply(effect_0.25_array, function(x) x[i]))), probs = 0.975))
    }
    effect_0.25_lower<-NULL
    for(i in 1:M){
      effect_0.25_lower<-rbind(effect_0.25_lower,quantile(as.numeric(unlist(lapply(effect_0.25_array, function(x) x[i]))), probs = 0.025))
    }
    cover_vector_0.25 = c(effect_0.25_upper >= check_para$effect_0.25 & check_para$effect_0.25 >= effect_0.25_lower)
    effect_0.25<-NULL
    for(i in 1:M){
      effect_0.25<-rbind(effect_0.25,mean(as.numeric(unlist(lapply(effect_0.25_array, function(x) x[i])))))
    }
    ## Distribution effects 0.75
    effect_0.75_upper<-NULL
    for(i in 1:M){
      effect_0.75_upper<-rbind(effect_0.75_upper,quantile(as.numeric(unlist(lapply(effect_0.75_array, function(x) x[i]))), probs = 0.975))
    }
    effect_0.75_lower<-NULL
    for(i in 1:M){
      effect_0.75_lower<-rbind(effect_0.75_lower,quantile(as.numeric(unlist(lapply(effect_0.75_array, function(x) x[i]))), probs = 0.025))
    }
    cover_vector_0.75 = c(effect_0.75_upper >= check_para$effect_0.75 & check_para$effect_0.75 >= effect_0.75_lower)
    effect_0.75<-NULL
    for(i in 1:M){
      effect_0.75<-rbind(effect_0.75,mean(as.numeric(unlist(lapply(effect_0.75_array, function(x) x[i])))))
    }
    effect_dist_upper<-NULL
    for(i in 1:4){
      effect_dist_upper<-rbind(effect_dist_upper,quantile(as.numeric(unlist(lapply(effect_dist_array, function(x) x[i]))), probs = 0.975))
    }
    effect_dist_lower<-NULL
    for(i in 1:4){
      effect_dist_lower<-rbind(effect_dist_lower,quantile(as.numeric(unlist(lapply(effect_dist_array, function(x) x[i]))), probs = 0.025))
    }
    cover_vector_dist = c(effect_dist_upper >= check_para$effect_dist & check_para$effect_dist >= effect_dist_lower)
    effect_dist<-NULL
    for(i in 1:4){
      effect_dist<-rbind(effect_dist,mean(as.numeric(unlist(lapply(effect_dist_array, function(x) x[i])))))
    }
    ## return the results for sim
    if(!is.null(dim(check_para$Alpha))){
      return(list(GMN = GMN, AMN = AMN, cover_gamma = cover_vector_gamma, cover_alpha = cover_matrix_alpha, 
                  sigma_delta_sq_mean = sigma_delta_sq_mean, cover_sigma_delta_sq = cover_sigma_delta_sq, 
                  Sigma_eta_mean = Sigma_eta_mean, cover_Sigma_eta = cover_matrix_Sigma_eta,
                  effect_mean = effect_mean, effect_median = effect_median, effect_0.25 = effect_0.25, 
                  effect_0.75 = effect_0.75, effect_dist = effect_dist,cover_mean = cover_vector_mean,
                  cover_median = cover_vector_median, cover_0.25 = cover_vector_0.25, cover_0.75 = cover_vector_0.75, 
                  cover_dist = cover_vector_dist
      ))
    }else{
      return(list(effect_mean = effect_mean, effect_median = effect_median, effect_0.25 = effect_0.25, 
                  effect_0.75 = effect_0.75, effect_dist = effect_dist,cover_mean = cover_vector_mean,
                  cover_median = cover_vector_median, cover_0.25 = cover_vector_0.25, cover_0.75 = cover_vector_0.75, 
                  cover_dist = cover_vector_dist
      ))
    }
    
  }
}

