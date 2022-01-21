

Bayesian_Semi_probit=function(Z,X,G=NULL,K=5,nIter=1100,
                              burnIn=100,thin=5,return_Bden=T,useCov=T,
                              true_para = NULL,
                              byStep=0.01, plotWhileRun = F,
                              check_para = NULL,zero_inflate = T, ref_b_spline = 1,
                              auto_check_convergence = F, max_Z = NULL,
                              c = 0.01, d = 0.01, alpha_prior_var = 10, 
                              gamma_prior_var=2.5,step_thres = 0.2, repara_lambda = F,
                              specify_knot = F, s_theta = 0.4, internal_knots, intercept = T) 
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
  #byStep=0.01, the step size for constructing the b densities, if Zs are integer, we could do stepsize  = 1
  #plotWhileRun=F：Plot out the convergence figures while running
  #check_para = NULL, then no computing rate for credible interval covering the true value. Or else
  #             this should be a list contain true value for Gamma (vector) and Alpha (matrix) with rows corresponding to covariates and columns to b-splines,
  #             Output would be a list contain two vector for Gamma (mean and whether CI covers the true value, 1 for yes), two matrix (mean for alpha and whether CI covers the true value)
  #zero_inflate = T, then we would have zero inflated part. If false, then only the b-spline part exists.
  #ref_b_spline = 1:this is the index for reference b-spline of Alpha coefficient
  #auto_check_convergence: if TRUE, auto checks convergence, if converged could end MCMC before the specified nIter
  #max_Z: the true maximum of Z, only validate in simulation, for B_density
  #c, d: the hyperparameters for both sigma_b_eta and sigma_b_delta prior
  #alpha_prior_var: the sigma_b_eta prior variance, a scaler so that we assume every alpha_k are the same prior
  #gamma_prior_var: prior variance for gamma 
  # step_thres: MH step size for threhold
  
  
  #Load packages/functions
  #library(tmvtnorm)
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
  thres_array = list()
  array_ind=0
  thres_accp = c()
  
  # Initialize global parameters 	
  pi0=0.5 
  gamma0=log((1-pi0)/pi0) 
  Gamma=array(0,dim=c(M,1)); Gamma[1]=gamma0
  Tau_sq=array(1,dim=c(M,1)) ; 
  if(no_G==F){
    sigma_b_delta_sq = 1
    Sigma_b_eta = 1
    Sigma_b = diag(c(1,1))
  }else{
    sigma_b_delta_sq = 0
    Sigma_b_eta = 0
  }
  b_delta = rnorm(n = n_sub, mean = 0, sd = sqrt(sigma_b_delta_sq))
  b_eta = rep(0, n_sub)
  
  
  Alpha=mvrnorm(n = 1, rep(0,M), diag(rep(1, M)))
  Alpha_mean=c();
  y = rep(0,N) # Latent variable for probit part 1
  thres = rep(0, K-1)
  thres[2:(K-1)] = 1:(K-2)/K
  #reparaterized thres
  theta = -(1:(K-2))/(K/4)
  #MH adaptive for theta
  ## Adaptive MH alogrithm
  C_t = list() # For b_eta list contain the covariance matrix n_sub elements with (K-1)*(K-1)
  theta_mean = theta # The previous iter mean, n_sub*(K-1)
  
  
  # Getting the b-spline density basis
  #browser()
  Gden<-matrix(0,nrow=N,ncol=K);
  if(specify_knot){
    if(is.null(internal_knots )){
      knots = quantile(Z[Z>0], probs = seq(0.1,0.999,length.out = K-2))
    }else{
      knots = internal_knots
    }
    
    knots = c(byStep, knots, max(abs(Z))+byStep)
  }else{
    knots = NULL
  }
  #browser()
  Gden<-Bdensity(byStep,mu = 0,K,N,Z,Gden, ref = ref_b_spline, max_Z = max_Z, knots = knots, intercept = intercept) # Get the density for Z
  Gden_raw <- Bdensity_raw(byStep,mu = 0,K,N,Z,ref = ref_b_spline, max_Z = max_Z, knots = knots, intercept= intercept) # return the original b-spline density
  
  # Check if model is non-zero inflated or there is no random effect
  if(zero_inflate == T){# zero inflation then add delta indicator
    Phi = ifelse(Z>0,1,0)
    P_phi = rep(1,N)
  }
  if(zero_inflate==F){
    Gamma=array(0,dim=c(M,1))
    Gamma[1]=100 # the intercept is very large, so that the non-null probability would be 1 when drawing b_eta
    b_delta = rep(0, n_sub)
    Phi = rep(1, N)
  }
  if(no_G == T){
    b_delta = rep(0, n_sub)
    b_eta = rep(0, n_sub)  # n_sub by K-1 matrix
  }
  
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
  l1 = rep(0,sum(Phi==1)) # Latent variable for probit part 2

  #get B-spline densities for Z1;
  gden<-matrix(0,nrow=length(Z1),ncol=K);
  gden<-Gden[nn,];
  
  
  #beginning of the iterations
  for(iter in 1:nIter){
    print(iter)
    
    if(is.null(true_para$Eta)==F){
      Eta = true_para$Eta
    }else{
      Eta<-rep(0,N);
      XA=X1%*%Alpha; #length(nn)*M times M*1
      j=1;
      for (i in nn){
        g_indx = G[i] # the group number for ith observation
        p_cum=c(pnorm(thres - XA[j,] - b_eta[g_indx]),1)
        p = c(p_cum[1], diff(p_cum)) * gden[j,]
        p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
        p_final = p/sum(p)
        Eta[i]=try(sample(c(1:K),size=1,replace=T,prob=p),silent=TRUE);
        if(length(grep("Error", Eta[i], fixed=TRUE))>0){
          Eta[i]=sample(c(1:K),size=1,replace=T,prob=rep(1/K,K))
        }
        j=j+1;
      }
      Eta = as.numeric(Eta)
    }
    
    
    # Draw latent variable y and l
    #browser()
    Eta1 = Eta[Phi==1]
    thres_aug = c(-Inf,thres, Inf)
    yl = list()
    V = list()
    X_list = list()#list for X that contain both gamma and alpha X
    para_join = c(Gamma, Alpha)
    if(no_G==F){# if there is grouping
      #H = Matrix(solve(objg$V[[1]]), sparse = T)# Sparse matrix for faster computation
      for(i in unique(G)){
        indx = G==i
        J = sum(indx)
        Phi_i = Phi[indx]
        X_i = X[indx,]
        lb = ifelse(Phi_i==1, 0, -Inf)
        ub = ifelse(Phi_i==1, Inf, 0)
        indx1 = G1==i
        J1 = sum(indx1)
        #browser()
        if(J1==0){
          if(J ==1){
            X_list[[i]] = t(as.matrix(c(X_i,rep(0,M))))
          }else{
            matrix_0 =  matrix(0, nrow = dim(X_i)[1], ncol = dim(X_i)[2])
            X_list[[i]] = cbind(X_i,matrix_0)
          }
          V[[i]] = diag(J)+rep(1,J)%*%t(rep(1,J))*Sigma_b[1,1]
          yl[[i]] = TruncatedNormal::rtmvnorm(1,lb = lb, ub= ub, mu = c(X_list[[i]]%*%para_join), sigma = V[[i]])
        }else{
          Eta_i = Eta1[indx1]
          X_1i = X1[indx1,]
          lb1 = thres_aug[Eta_i]
          ub1 = thres_aug[(Eta_i+1)]
          w = matrix(
            cbind(c(rep(1,J), rep(0,J1)), c(rep(0,J), rep(1,J1)) ), ncol = 2
                     )
          if(J1 ==1 & J>1){
            matrix_0 =  matrix(0, nrow = dim(X_i)[1], ncol = dim(X_i)[2])
            X_list[[i]] = as.matrix(rbind(cbind(X_i,matrix_0), c(rep(0,M), X_1i)))
          }else if(J1>1 & J>1){
            matrix_0 =  matrix(0, nrow = dim(X_i)[1], ncol = dim(X_i)[2])
            matrix_01 =  matrix(0, nrow = dim(X_1i)[1], ncol = dim(X_1i)[2])
            X_list[[i]] = as.matrix(rbind(cbind(X_i,matrix_0), cbind(matrix_01, X_1i)))
          }else if(J1==1 & J==1){
            X_list[[i]] = as.matrix(rbind(c(X_i,rep(0,M)), c(rep(0,M), X_1i)))
          }
          
          V[[i]] = diag(J+J1)+w%*%Sigma_b%*%t(w)
          yl[[i]] = TruncatedNormal::rtmvnorm(1,lb = c(lb,lb1), ub= c(ub,ub1), mu = c(X_list[[i]]%*%para_join), sigma = V[[i]])
        }
      }
    }else{# no random effect case
      for(g in unique(G)){
        indx = G==g
        Phi_i = Phi[indx]
        X_i = X[indx,]
        lb = ifelse(Phi_i==1, 0, -Inf)
        ub = ifelse(Phi_i==1, Inf, 0)
        #y_temp = tmvtnorm::rtmvnorm(1, mean = c(X_i%*%Gamma),sigma = objg$V[[i]], lower = lb, upper = ub, algorithm = "gibbs")
        y_temp = truncnorm::rtruncnorm(1,a = lb, b= ub, mean = c(X_i%*%Gamma),sd = 1)#no random effect the variance is always 1
        #y_temp = tmvtnorm::rtmvnorm(1, mean = c(X_i%*%Gamma),sigma = objg$V[[i]], lower = lb, upper = ub,alogrithm="gibbs",burn.in.samples = 300)
        y[indx] = y_temp
      }
      l1 = rep(0, length(G1))
      for(g in unique(G1)){
        indx = G1==g
        Eta_i = Eta1[indx]
        X_i = X1[indx,]
        lb = thres_aug[Eta_i]
        ub = thres_aug[(Eta_i+1)]
        #y_temp = tmvtnorm::rtmvnorm.sparseMatrix(1, mean = c(X_i%*%Gamma), H = Matrix(solve(objg$V[[i]]), sparse = T), lower = lb, upper = ub,burn.in.samples = 100)
        l_temp = truncnorm::rtruncnorm(1,a = lb, b= ub, mean = c(X_i%*%Alpha),sd = 1)
        l1[indx] = l_temp
      }
    }

    
    # Draw Gamma and Alpha
    if(no_G == F){
      T_1 = list()
      T_2 = list()
      Sigma_gamma_alpha = diag(c(rep(gamma_prior_var, length(Gamma)),
                                 rep(alpha_prior_var,length(Alpha))))
      for(i in unique(G)){
        V_i_inv = solve(V[[i]])
        T_1[[i]] = t(X_list[[i]])%*%V_i_inv%*%X_list[[i]]
        T_2[[i]] = t(X_list[[i]])%*%V_i_inv%*%(yl[[i]])
      }
      var = solve(solve(Sigma_gamma_alpha)+Reduce("+",T_1))
      mean = var%*%Reduce("+",T_2)
      Gamma_Alpha = MASS::mvrnorm(n=1,mu = mean,Sigma=var)
      Gamma = Gamma_Alpha[1:length(Gamma)]
      Alpha = Gamma_Alpha[(length(Gamma)+1):length(Gamma_Alpha)]
    }else{
      objg=Draw_Gamma_log_M(y=y, X=X,G=G,SG = Sigma_Gamma, Phi =Phi,b_delta=b_delta, same_length_G = same_length_G)
      Gamma=objg$par;
      if(is.null(true_para$Alpha)==F){#For model diagnostics, plug in true parameter
        Alpha = true_para$Alpha
      }else{
        obj=Draw_Alpha_probit(l1 = l1, X1 = X1, G1 = G1, 
                              Sigma_alpha = diag(rep(alpha_prior_var, M)), 
                              b_eta = b_eta)		
        Alpha=obj$par;	
      }
    }

    
    #browser()
    # Draw b_eta and b_delta
    #browser()
    if(no_G == F){#zero inflated and random effect
      if(is.null(true_para$b_delta) == F){#For model diagnostics, plug in true parameter
        b_delta = true_para$b_delta
        b_eta = true_para$b_eta
      }else{
        for(i in 1:n_sub){
          obj_b = Draw_b_delta_eta_i_probit_wishart(i, G, G1, X, X1, Gamma,Alpha, Sigma_b, Phi, yl_i = yl[[i]])
          b_delta[i] = obj_b$par[1]
          b_eta[i] = obj_b$par[2]
        }
      }
    }
    
    
    # Draw Sigma_b
    if(no_G == F){
      if(is.null(true_para$Sigma_b)){
        Sigma_b = MCMCpack::riwish(n_sub + 2, diag(1, 2)+t(matrix(cbind(b_delta, b_eta), ncol = 2))%*%matrix(cbind(b_delta, b_eta), ncol = 2))
        sigma_b_delta_sq = Sigma_b[1,1]
        Sigma_b_eta = Sigma_b[2,2]
      }else{
        Sigma_b = true_para$Sigma_b
        sigma_b_delta_sq = Sigma_b[1,1]
        Sigma_b_eta = Sigma_b[2,2]
      }
    }
    
   
    if(repara_lambda){
      #browser()
      theta_accp = rep(0, K-2)
      #Sigma_theta = diag(c(rep(0.005, ceiling((K-2)*0.2)), rep(0.01, K-2-ceiling((K-2)*0.2))))
      Sigma_theta = diag(0.01, nrow = K-2)
      #calculate the MH covariance
      if(iter == 1){
        C_t = s_theta * Sigma_theta + s_theta*0.0001*diag(1, nrow = K-2)
      }else{
        t = iter - 1
        theta_mean_old = theta_mean #old mean b_eta
        theta_mean = (theta_mean_old * (iter-1) + theta)/iter # new mean
        C_t = (t-1)*C_t/t + s_theta*(t*theta_mean_old%*%t(theta_mean_old)-
                                       (t+1)*theta_mean%*%t(theta_mean)+
                                       theta%*%t(theta)+
                                       0.0001*diag(1, nrow = K-2))/t
      }
      #browser()
      # multiple = 2
      # log_p_theta_star=rep(0,multiple)
      # log_p_theta_2star=rep(0,multiple)
      #theta_opt <-tryCatch({optim(theta, Eta1 = Eta1, G1 = G1, X1 = X1, Alpha = Alpha, b_eta = b_eta,K = K, Sigma_theta = Sigma_theta,
      #                            log_thres_p_repar,method="Nelder-Mead",
      #                            hessian=TRUE,control=list(maxit=10,fnscale=-1))},error=function(e) NULL)
      #browser()
      # if(!is.null(theta_opt)){
      #   theta = theta_opt$par
      #   C_t = solve(-theta_opt$hessian)
      #   theta_star = mvtnorm::rmvt(n = multiple, delta = theta,sigma = C_t*0.7, df = 5)
      # }else{
      #   theta_star = MASS::mvrnorm(n = multiple, mu = theta,Sigma = C_t)
      # }
      # 
      # theta_star = MASS::mvrnorm(n = multiple, mu = theta,Sigma = C_t)
      # for (j in 1:multiple){
      #   log_p_theta_star[j]=log_thres_p_repar(theta_star[j,], Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
      # }
      # d_p_theta_star = mvtnorm::dmvnorm(theta_star, mean = theta, sigma = C_t, log = T)
      # if(any(is.na(log_p_theta_star))){
      #   print("overflow in theta density")
      # }
      # # browser()
      # max_p_theta = median(log_p_theta_star)
      # p=exp(log_p_theta_star-d_p_theta_star-max_p_theta)/sum(exp(log_p_theta_star - d_p_theta_star-max_p_theta))
      # #in case there is still overfloat;
      # p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
      # v=sample(c(1:multiple),1,prob=p);
      # 
      # theta_2star = MASS::mvrnorm(n = multiple-1, mu = theta_star[v,],Sigma = C_t)
      # theta_2star <-rbind(theta_2star,theta)
      # 
      # for (j in 1:multiple){
      #   log_p_theta_2star[j]=log_thres_p_repar(theta_2star[j,], Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
      # }
      # d_p_theta_2star = mvtnorm::dmvnorm(theta_2star, mean = theta_star[v,], sigma = C_t, log = T)
      # if(any(is.na(log_p_theta_2star))){
      #   print("overflow in theta density")
      # }
      # 
      # #control overfloat
      # num=sum(exp(log_p_theta_star -d_p_theta_star-max_p_theta))
      # den=sum(exp(log_p_theta_2star -d_p_theta_2star-max_p_theta))
      # 
      # rho=min(1,num/den)
      # 
      # #in case overfloat again
      # if(is.na(rho)) {
      #   print("theta rho is na")
      #   rho=0.01}
      # 
      # accp = 0
      # u = log(runif(1, min = 0, max = 1))
      # if(u < rho){
      #   theta = theta_star[v,]
      #   accp = 1
      #   for(j in 2:(K-1)){
      #     thres[j] = sum(exp(theta[1:(j-1)]))
      #   }
      # }
      
      for(k in 1:(K-2)){
        log_p_old = log_thres_p_repar(theta, Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
        theta_star = theta
        theta_star[k] =rnorm(n = 1, mean = theta[k],sd = sqrt(C_t[k,k]))
        log_p_new = log_thres_p_repar(theta_star, Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
        rho = min(0, log_p_new - log_p_old)
        accp = 0
        u = log(runif(1, min = 0, max = 1))
        if(is.na(rho)){
          rho = 0.5
          print("Thres probability NA")
        }
        if(u < rho){
          theta = theta_star
          accp = 1
        }
      }
      for(j in 2:(K-1)){
        thres[j] = sum(exp(theta[1:(j-1)]))
      }
      theta_accp = c(theta_accp, mean(accp))
      #theta_opt <-tryCatch({optim(theta, Eta1 = Eta1, G1 = G1, X1 = X1, Alpha = Alpha, b_eta = b_eta,K = K, Sigma_theta = Sigma_theta,
                                  # log_thres_p_repar,method="Nelder-Mead",
                                  # hessian=TRUE,control=list(maxit=10,fnscale=-1))},error=function(e) NULL)
      # if(!is.null(theta_opt)){
      #   theta = theta_opt$par
      #   #C_t = solve(-theta_opt$hessian)
      #   theta_star = c(mvtnorm::rmvt(n = 1, delta = theta,sigma = C_t*0.7, df = 5))
      # }else{
      #   theta_star = MASS::mvrnorm(n = 1, mu = theta,Sigma = C_t)
      # }
      # log_p_old = log_thres_p_repar(theta, Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
      # log_p_new = log_thres_p_repar(theta_star, Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
      # rho = min(0, log_p_new - log_p_old)
      # 
      # #in case overfloat again
      # if(is.na(rho)) {
      #   print("theta rho is na")
      #   rho=0.05}
      # 
      # accp = 0
      # u = log(runif(1, min = 0, max = 1))
      # if(u < rho){
      #   theta = theta_star[]
      #   accp = 1
      #   for(j in 2:(K-1)){
      #     thres[j] = sum(exp(theta[1:(j-1)]))
      #   }
      # }
    }else{
      accep_v = c()
      for(k in 2:(K-1)){
        thres_aug = c(-Inf,thres, Inf)
        r = step_thres
        log_p_old = log_thres_p(thres_aug, Eta1, G1, X1, Alpha, b_eta)
        thres_star = runif(1, min = max(thres_aug[k], thres_aug[k+1]-r),
                           max = min(thres_aug[k+2], thres_aug[k+1]+r))
        thres_aug_star = thres_aug
        thres_aug_star[k+1] =  thres_star
        log_p_new = log_thres_p(thres_aug_star, Eta1, G1, X1, Alpha, b_eta)
        rho = min(0, log_p_new - log_p_old)
        accp = 0
        u = log(runif(1, min = 0, max = 1))
        if(is.na(rho)){
          rho = 0.5
          print("Thres probability NA")
        }
        if(u < rho){
          thres[k] = thres_star
          accp = 1
        }
        accep_v = c(accep_v, accp)
      }
      thres_accp = c(thres_accp, mean(accep_v))
      # for(k in 2:(K-1)){
      #   low = max(0,max(l1[Eta1 == k]), thres[k-1],-20)
      #   up = min(min(l1[Eta1 == (k+1)]), thres[(k+1)],20, na.rm = T)#na.rm for the last thres[(k+1)], out of bound, gives NA
      #   thres[k] = runif(1, min = low, max = up)
      # }
    }
    
    
  
    # browser()
    # multiple = 2
    # log_p_theta_star=rep(0,multiple)
    # log_p_theta_2star=rep(0,multiple)
    # 
    # theta_star = MASS::mvrnorm(n = multiple, mu = theta,Sigma = C_t)
    # 
    # for (j in 1:multiple){
    #   log_p_theta_star[j]=log_thres_p_repar(theta_star[j,], Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
    # }
    # d_p_theta_star = mvtnorm::dmvnorm(theta_star, mean = theta, sigma = C_t, log = T)
    # if(any(is.na(log_p_theta_star))){
    #   print("overflow in theta density")
    # }
    # # browser()
    # max_p_theta = median(log_p_theta_star)
    # p=exp(log_p_theta_star-d_p_theta_star-max_p_theta)/sum(exp(log_p_theta_star - d_p_theta_star-max_p_theta))
    # #in case there is still overfloat;
    # p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
    # v=sample(c(1:multiple),1,prob=p);
    # 
    # theta_2star = MASS::mvrnorm(n = multiple-1, mu = theta_star[v,],Sigma = C_t)
    # theta_2star <-rbind(theta_2star,theta)
    # 
    # for (j in 1:multiple){
    #   log_p_theta_2star[j]=log_thres_p_repar(theta_2star[j,], Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
    # }
    # d_p_theta_2star = mvtnorm::dmvnorm(theta_2star, mean = theta_star[v,], sigma = C_t, log = T)
    # if(any(is.na(log_p_theta_2star))){
    #   print("overflow in theta density")
    # }
    # 
    # #control overfloat
    # num=sum(exp(log_p_theta_star -d_p_theta_star-max_p_theta))
    # den=sum(exp(log_p_theta_2star -d_p_theta_2star-max_p_theta))
    # 
    # rho=min(1,num/den)
    # 
    # #in case overfloat again
    # if(is.na(rho)) {
    #   print("theta rho is na")
    #   rho=0.01}
    # 
    # accp = 0
    # u = log(runif(1, min = 0, max = 1))
    # if(u < rho){
    #   theta = theta_star[v,]
    #   accp = 1
    #   for(j in 2:(K-1)){
    #     thres[j] = sum(exp(theta[1:(j-1)]))
    #   }
    # }
    # 
    #browser()
    
    #browser()
    # for(k in 1:(K-2)){
    #   log_p_old = log_thres_p_repar(theta, Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
    #   theta_star = theta
    #   theta_star[k] =rnorm(n = 1, mean = theta[k],sd = sqrt(C_t[k,k]))
    #   log_p_new = log_thres_p_repar(theta_star, Eta1, G1, X1, Alpha, b_eta,K = K, Sigma_theta = Sigma_theta)
    #   rho = min(0, log_p_new - log_p_old)
    #   accp = 0
    #   u = log(runif(1, min = 0, max = 1))
    #   if(is.na(rho)){
    #     rho = 0.5
    #     print("Thres probability NA")
    #   }
    #   if(u < rho){
    #     theta = theta_star
    #     accp = 1
    #   }
    # }
    # for(j in 2:(K-1)){
    #   thres[j] = sum(exp(theta[1:(j-1)]))
    # }
    # theta_accp = c(theta_accp, mean(accp))
    
    
    
    
    
    
    
    # Save posterior sample in array once >burn in for each thin
    if(iter%%thin==0 & iter>=burnIn){
      array_ind=array_ind+1
      dist_effect = unlist(cal_quantity_probit(NULL, x = rep(1, length(Gamma)), Alpha, Gamma, 
                                               Gden_raw$phi.mat, Gden_raw$grid, byStep, K = K, thres = thres))
      effect = coeff_effect_probit(NULL, Alpha, Gamma, 
                                   Gden_raw$phi.mat, Gden_raw$grid, byStep, K = K, thres = thres)
      effect_dist_array[[array_ind]] = dist_effect
      effect_mean_array[[array_ind]] = effect$effect_mean
      effect_median_array[[array_ind]] = effect$effect_median
      effect_0.25_array[[array_ind]] = effect$effect_0.25
      effect_0.75_array[[array_ind]] = effect$effect_0.75
      ALPHA_array[[array_ind]]=Alpha
      thres_array[[array_ind]]=thres
      if(zero_inflate==T){
        GAMMA_array[[array_ind]]=Gamma
        if(no_G==F){
          b_delta_array[[array_ind]]=b_delta
          sigma_delta_sq_array[[array_ind]]=sigma_b_delta_sq
        }
      }
      Tau_sq_array[[array_ind]]=Tau_sq
      if(no_G==F){
        b_eta_array[[array_ind]]=b_eta
        Sigma_eta_array[[array_ind]]=Sigma_b_eta
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
  
  
  #return the posterior sampling if we don't have the true parameters to check for coverage and bias
  if(is.null(check_para)){
    results=list()
    results[["alpha_array"]]=ALPHA_array
    results[["thres_array"]]=thres_array
    if(zero_inflate==T){
      results[["gamma_array"]]=GAMMA_array
      results[["b_delta"]]=b_delta_array
      results[["sigma_b_delta_sq"]]=sigma_delta_sq_array
      results[["accp_rate_b_delta"]]=Accp_Rate_array_b_delta
      results[["accp_rate_array_g"]]=Accp_Rate_array_g
      results[["p_phi"]]=P_phi
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
    results[["knots"]]=Gden_raw$knots
    results[["eta"]] = Eta #last iteration Eta value
    results[["zero_inflate"]]=zero_inflate
    results[["grouped"]]=!no_G
    results[["col_names"]]=col_names
    results[["effect_dist"]] = effect_dist_array
    results[["effect_mean"]]=effect_mean_array 
    results[["effect_median"]]=effect_median_array
    results[["effect_0.25"]]=effect_0.25_array 
    results[["effect_0.75"]]=effect_0.75_array
    results[["thres_accp"]] = thres_accp
    return(results)
  }else{#If there is true parameter vector, return coverage T/F and mean of posterior
    ##Alpha 
    AMN_upper<-NULL
    for(i in 1:M){
      AMN_upper<-rbind(AMN_upper,quantile(as.numeric(unlist(lapply(ALPHA_array, function(x) x[i]))), probs = 0.975))
    }
    AMN_lower<-NULL
    for(i in 1:M){
      AMN_lower<-rbind(AMN_lower,quantile(as.numeric(unlist(lapply(ALPHA_array, function(x) x[i]))), probs = 0.025))
    }
    cover_vector_alpha = c(AMN_upper >= check_para$Alpha & check_para$Alpha >= AMN_lower)
    AMN<-NULL
    for(i in 1:M){
      AMN<-rbind(AMN,mean(as.numeric(unlist(lapply(ALPHA_array, function(x) x[i])))))
    }
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
    # Sigma_eta
    if(no_G==F){
      Sigma_eta_lower = quantile(unlist(Sigma_eta_array), probs = 0.025)
      Sigma_eta_upper = quantile(unlist(Sigma_eta_array), probs = 0.975)
      cover_matrix_Sigma_eta = Sigma_eta_upper >= check_para$Sigma_b_eta^2 & check_para$Sigma_b_eta^2 >= Sigma_eta_lower
      Sigma_eta_mean = mean(unlist(Sigma_eta_array))
    }else{
      Sigma_eta_lower = NA
      Sigma_eta_upper = NA
      cover_matrix_Sigma_eta = diag(1, K-1)
      Sigma_eta_mean = diag(0, K-1)
    }
    ##If zero-inflated then calculate Gamma and sigma_delta mean and coverage
    if(zero_inflate==T){
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
    if(!is.null(check_para$thres)){
      thres_upper<-NULL
      for(i in 1:(K-1)){
        thres_upper<-rbind(thres_upper,quantile(as.numeric(unlist(lapply(thres_array, function(x) x[i]))), probs = 0.975))
      }
      thres_lower<-NULL
      for(i in 1:(K-1)){
        thres_lower<-rbind(thres_lower,quantile(as.numeric(unlist(lapply(thres_array, function(x) x[i]))), probs = 0.025))
      }
      cover_vector_thres = c(thres_upper >= check_para$thres & check_para$thres >= thres_lower)
      thres_mean <- NULL
      for(i in 1:(K-1)){
        thres_mean<-rbind(thres_mean,mean(as.numeric(unlist(lapply(thres_array, function(x) x[i])))))
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
    if(!is.null(check_para$thres)){
      return(list(GMN = GMN, AMN = AMN, thres_mean = thres_mean, cover_gamma = cover_vector_gamma, cover_alpha = cover_vector_alpha, 
                  cover_thres = cover_vector_thres, sigma_delta_sq_mean = sigma_delta_sq_mean, cover_sigma_delta_sq = cover_sigma_delta_sq, 
                  Sigma_eta_mean = Sigma_eta_mean, cover_Sigma_eta = cover_matrix_Sigma_eta,
                  effect_mean = effect_mean, effect_median = effect_median, effect_0.25 = effect_0.25, 
                  effect_0.75 = effect_0.75, effect_dist = effect_dist,cover_mean = cover_vector_mean,
                  cover_median = cover_vector_median, cover_0.25 = cover_vector_0.25, cover_0.75 = cover_vector_0.75, 
                  cover_dist = cover_vector_dist
      ))
    }else{
      return(list(GMN = GMN,cover_gamma = cover_vector_gamma, AMN = AMN,cover_alpha = cover_vector_alpha, 
                  sigma_delta_sq_mean = sigma_delta_sq_mean, cover_sigma_delta_sq = cover_sigma_delta_sq, 
                  Sigma_eta_mean = Sigma_eta_mean, cover_Sigma_eta = cover_matrix_Sigma_eta,
                  effect_mean = effect_mean, effect_median = effect_median, effect_0.25 = effect_0.25, 
                  effect_0.75 = effect_0.75, effect_dist = effect_dist,cover_mean = cover_vector_mean,
                  cover_median = cover_vector_median, cover_0.25 = cover_vector_0.25, cover_0.75 = cover_vector_0.75, 
                  cover_dist = cover_vector_dist
      ))
    }
  }
}
