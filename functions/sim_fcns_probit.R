# Simulation functions
library(matrixStats)
# Generate K b-spline density basis
Bdensity_sim = function(byStep=0.01, K=5, end = 150, degree = 3, intercept = T){
  grid=seq(byStep,end+byStep,by=byStep)
  if(K>4){
    knots=seq(byStep,end+byStep,len=K) #equal space;This includes starting and ending point
    Basis.mat=bs(grid,knots=knots[2:(length(knots)-1)],degree=degree,intercept=intercept,
                 Boundary.knots=c(knots[1],knots[length(knots)]))
    phi.mat=Basis.mat[,1:K]# Get rid of the last two basis since usually this density would be low at high value of z
    
    #normalization to turn b-spline basis into b-spline density
    for(k in 1:K){
      phi.mat[,k]= phi.mat[,k]/(sum(phi.mat[,k])*byStep)
    }
  }else{
    #no knot, only 4 b spline, minus 1, use only 3
    Basis.mat=bs(grid,degree=degree,intercept=T)
    phi.mat=Basis.mat[,1:3]
    phi.mat[,1] = phi.mat[,1]/(end/4)
    phi.mat[,2] = phi.mat[,2]/(end/4)
    phi.mat[,3] = phi.mat[,3]/(end/4)
  }
  
  # phi.mat[,3:(K-2)]=phi.mat[,3:(K-2)]/(knots[2]-knots[1])
  # phi.mat[,1]=phi.mat[,1]*2/(knots[2]-knots[1])
  # phi.mat[,2]=phi.mat[,2]*4/(3*(knots[2]-knots[1]))
  # phi.mat[,(K-1)]=phi.mat[,(K-1)]*4/(3*(knots[2]-knots[1]))
  # phi.mat[,K]=phi.mat[,K]*2/(knots[2]-knots[1]);
  return(phi.mat)
}

# Function for visualize the B spline plot
draw_b_spline = function(phi.mat, byStep=0.01,
                         end = 150, draw_mix = F,prior=c(0.2,0.1,0.1,0.1,0.5)){
  grid=seq(byStep,end+byStep,by=byStep)
  plot(phi.mat[,1]~grid, ylim=c(0,max(phi.mat)), type='l', lwd=2, col=1, 
       xlab="Cubic B-spline basis", ylab="")
  for (j in 2:ncol(phi.mat)) lines(phi.mat[,j]~grid, lwd=2, col=j)
  #Drawing the combined density
  if(draw_mix == T){
    for(i in 1:ncol(phi.mat)){
      phi.mat[,i] = phi.mat[,i] * prior[i]
    }
    lines(rowSums(phi.mat)~grid, lwd=4, col=ncol(phi.mat)+1)
  }
  
}



#prior = c(0.1,0.2,0.3,0.3,0.1)
#Generate observations from an empirical b-spline density distribution 
gen_z_f1_probit = function(u, eta, byStep=0.01, phi.mat, prior,end = 150){
  #u is vector of numbers generated from U(0,1)
  #phi.mat is a empirical b-spline basis
  #prior is a vector of probablities of observation coming from different b-spline densities
  grid=seq(byStep,end+byStep,by=byStep)
  # multiple the prior to each basis
  # for(i in 1:ncol(phi.mat)){
  #   phi.mat[,i] = phi.mat[,i] * prior[i]
  # }
  # CDF at each point
  # gen_num = rep(0, length(u))
  # cdf = cumsum(rowSums(phi.mat*byStep))
  # for(j in 1:length(u)){
  #   gen_num[j] = grid[which.min(abs(u[j] - cdf))]
  # }
  # which density to generate from
  #Getting the cdf for each density
  for(i in 1:ncol(phi.mat)){
    phi.mat[,i] = cumsum(phi.mat[,i]*byStep)
  }
  cdf = phi.mat
  # start to generate number
  gen_num = rep(0, length(u))
  for(j in 1:length(u)){
    ind_density = eta# density index
    gen_num[j] = grid[which.min(abs(u[j] - cdf[,ind_density]))]
  }
  return(gen_num)
}

#test
# u = runif(10) 
# gen_z_f1(u = u, phi.mat = phi.mat, prior = prior, Z = Z)

# Function to compare the true b-splines and the model output b-splines
b_spline_compare = function(phi.mat.true, phi.mat.model,byStep=0.01,
                            end = 150, true_prior, model_prior, G_indx = NULL, 
                            legend = T,xlab_T = T, ylab_T = T, title_T = T,multiple = T, CI = NULL){
  grid=seq(byStep,end+byStep,by=byStep)
  
  # the true density for b-splines with the prior prob
  for(i in 1:ncol(phi.mat.true)){
    phi.mat.true[,i] = phi.mat.true[,i] * true_prior[i]
  }
  
  # the model density for b-splines with the model prob
  for(i in 1:ncol(phi.mat.model)){
    phi.mat.model[,i] = phi.mat.model[,i] * model_prior[i]
  }
  
  if(multiple == T){
    plot(rowSums(phi.mat.true)~grid, ylim=c(0,max(phi.mat)), type='l', lty=1, col="red",
         xlab = NULL, ylab = NULL
    )
    if(xlab_T ==T){
      title(xlab="Observation value", ylab="Density")
    }
    if(ylab_T == T){
      title(ylab="Density")
    }
    if(title_T == T){
      title(main = G_indx)
    }
    lines(rowSums(phi.mat.model)~grid, lty=2, col="blue")
    if(legend == T){
      legend("topright", legend=c("True density", "Model density"),
             col=c("red", "blue"), lty=1:2, cex=0.8, ncol = 2)
    }
  }else{
    g = ggplot() + geom_line(aes(y = rowSums(phi.mat.true), x = grid, color = "True")) + 
      geom_line(aes(y = rowSums(phi.mat.model), x = grid, color = "Model"))+
      geom_ribbon(aes(x = grid,ymin = CI[,1], ymax = CI[,2]), alpha = 0.3)+
      theme_bw() 
      
    if(legend == T){
      g = g + theme(legend.position = "bottom")+labs(x = "Observation values",
          y = "Density",
          title = "Population distribution for positive part comparing model to true",
          color = "")
    }else{
      g= g+  theme(legend.position="none",
                     axis.title = element_blank())
    }
    return(g)
  }
  
}

# CI for the positive density
cal_CI_B_density_probit = function(MCMCfit, x, last_n = 1000, g = NULL){
  alpha_array = tail(MCMCfit$alpha_array, last_n)
  list_density = list()
  phi.mat = MCMCfit$Bden_raw
  K = dim(phi.mat)[2]
  thres_array = tail(MCMCfit$thres_array, last_n)
  if(!is.null(g)){
    b_eta_array = tail(MCMCfit$b_eta, last_n)
  }
  for(j in 1:last_n){
    c_vec = c()
    est_Alpha = alpha_array[[j]]
    est_thres = thres_array[[j]]
    if(!is.null(g)){
      b_eta = b_eta_array[[j]][g]
      for(k in 1:K){
        if(k==1){
          c_vec = c(c_vec, pnorm(est_thres[k] - x %*% est_Alpha - b_eta))
        }else if(k==K){
          c_vec = c(c_vec, 1 - pnorm(est_thres[k-1] - x %*% est_Alpha - b_eta))
        }else{
          c_vec = c(c_vec, pnorm(est_thres[k] - x %*% est_Alpha - b_eta) - sum(c_vec) )
        }
      }
    }else{
      for(k in 1:K){
        if(k==1){
          c_vec = c(c_vec, pnorm(est_thres[k] - x %*% est_Alpha ))
        }else if(k==K){
          c_vec = c(c_vec, 1 - pnorm(est_thres[k-1] - x %*% est_Alpha ))
        }else{
          c_vec = c(c_vec, pnorm(est_thres[k] - x %*% est_Alpha ) - sum(c_vec) )
        }
      }
    }
    
    prior_est = c_vec
    phi.mat.model = phi.mat
    for(i in 1:ncol(phi.mat)){
      phi.mat.model[,i] = phi.mat[,i] * prior_est[i]
    }
    list_density[[j]] = rowSums(phi.mat.model)
  }
  temp_df = do.call(rbind,list_density)
  return(colQuantiles(temp_df, probs = c(0.025, 0.975)))
  
}


# Function to plot the density on the positive part comparing estimate and true
positive_sim_plot_probit = function(x, true_Alpha, est_Alpha, phi.mat.true, phi.mat.model,byStep = 0.01, 
                             G = NULL, true_b_eta = NULL, MCMCfit, X = NULL,
                             G_indx = NULL, last_n = 1000, true_thres, est_thres){
  #x is the vector for the evaluating x vector
  # true_Alpha: a matrix of true alpha, row number equal to number of covariates, column number of b-spline
  #phi.mat: the b-spline generating basis 
  #byStep: 
  #G: the group index for each observation, if it is not null, then we draw 9 plots for 9 random subject based on their ture b_eta and model b_eta
  #true_b_eta: a matrix of #subject*(K-1)
  #MCMCfit: model fit
  #X: the matrix of covariates we used to fit the model
  #G_indx: the group index we are going to plot, recommend 4 or 9 indexes 
  K_true = dim(phi.mat.true)[2]
  K_model = dim(phi.mat.model)[2]
  if(is.null(G)){
    CI = cal_CI_B_density_probit(MCMCfit, x, last_n = last_n)
    c_vec = c()
    for(k in 1:K_true){
      if(k==1){
        c_vec = c(c_vec, pnorm(true_thres[k] - x %*% true_Alpha ))
      }else if(k==K_true){
        c_vec = c(c_vec, 1 - pnorm(true_thres[k-1] - x %*% true_Alpha ))
      }else{
        c_vec = c(c_vec, pnorm(true_thres[k] - x %*% true_Alpha) - sum(c_vec) )
      }
    }
    prior_true = c_vec
    
    
    c_vec = c()
    for(k in 1:K_model){
      if(k==1){
        c_vec = c(c_vec, pnorm(est_thres[k] - x %*% est_Alpha ))
      }else if(k==K_model){
        c_vec = c(c_vec, 1 - pnorm(est_thres[k-1] - x %*% est_Alpha ))
      }else{
        c_vec = c(c_vec, pnorm(est_thres[k] - x %*% est_Alpha) - sum(c_vec) )
      }
    }
    prior_est = c_vec
    if(length(prior_est)==length(prior_true)){
      d_result = data.frame(parameters = paste("b_prior", 1:length(c_vec)), 
                          True = c(prior_true),
                          Model = c(prior_est ))
    }else{
      d_result = NULL
    }
    
    return(list(b_spline_compare(phi.mat.true, phi.mat.model, byStep=byStep,
                                 end = 150, prior_true, prior_est, multiple = F, CI = CI), d_result) )
  }else{
    if(is.null(G_indx)){
      i_sample = sample(unique(as.numeric(as.factor(G))), 9)
    }else{
      i_sample = G_indx
    }
    par(mfrow = c(3,3), xpd = TRUE)
    p_list = list()
    p_list = lapply(1:length(i_sample), function(j){
      i = i_sample[j]
      b_eta_true = true_b_eta[i]
      x_test = X[G == sort(unique(G))[i],]
      x_test = x_test[1,]
      c_vec = c()
      for(k in 1:K_true){
        if(k==1){
          c_vec = c(c_vec, pnorm(true_thres[k] - x_test %*% true_Alpha - b_eta_true ))
        }else if(k==K){
          c_vec = c(c_vec, 1 - pnorm(true_thres[k-1] - x_test %*% true_Alpha - b_eta_true ))
        }else{
          c_vec = c(c_vec, pnorm(true_thres[k] - x_test %*% true_Alpha - b_eta_true) - sum(c_vec) )
        }
      }
      prior_true = c_vec
      
      
      # Mean of b_eta model
      b_eta_model = mean(sapply(tail(MCMCfit$b_eta, last_n), 
                                         function(x) x[i])) 
      c_vec = c()
      for(k in 1:K_model){
        if(k==1){
          c_vec = c(c_vec, pnorm(est_thres[k] - x_test %*% est_Alpha - b_eta_model ))
        }else if(k==K){
          c_vec = c(c_vec, 1 - pnorm(est_thres[k-1] - x_test %*% est_Alpha - b_eta_model ))
        }else{
          c_vec = c(c_vec, pnorm(est_thres[k] - x_test %*% est_Alpha - b_eta_model) - sum(c_vec) )
        }
      }
      prior_est = c_vec
      CI = cal_CI_B_density_probit(MCMCfit,x_test, last_n = last_n, g = i)
      b_spline_compare(phi.mat.true, phi.mat.model,byStep=byStep,
                        end = 150, prior_true, prior_est, G_indx = i, legend = F,
                        xlab_T = F, ylab_T = T, title_T = T, multiple = F, CI = CI)
      
    })
    return(grid.arrange(grobs = p_list, ncol = 3))
  }
}



#function for simulating data from a probit ordinal regression
sim_data_probit = function(n, J, gamma, Alpha, sigma_b, phi.mat, x_type = NULL, byStep = 0.01,
                           K, thres,random_effect_equal = F){
  # simulate data
  # x_type: a vector of either "numeric" or "binary", same length as length(gamma)-1
  # thres: threhold for the ordinal category of splines, should be K-1
  # K: number of b-spline density basis
  # Global indicator random effect 
  b =  mvrnorm(n, mu = rep(0,2), Sigma = sigma_b) 
  b_delta = b[,2]
  # Local indicator covraince matrix (generate the vector first regardless of the delta=1 or not)
  if(random_effect_equal){
    b_eta = b_delta#equal to b_delta
  }else{
    b_eta = b[,1]
  }
  # n by K-1 matrix
  # Predictor
  x_length = length(gamma)
  X = list()
  for(i in 1:n){
    if(is.null(x_type)){
      X[[i]] = cbind(rep(1,J), mvrnorm(J, mu = rep(0, x_length-1), Sigma = diag(x_length-1)),i)#, runif(n, min = -1, max = 5)
    }else{
      c = rep(1,J)
      for(type in x_type){
        if(type=="numeric"){
          c = cbind(c, rnorm(J))
        }else if(type == "binary"){
          c = cbind(c, rbinom(J, 1,0.5))
        }
      }
      X[[i]] = cbind(c,i)
    } 
  }
  
  
  X = do.call(rbind, X) #combine all the different subjects into one matrix
  thres_aug = c(-Inf, thres, Inf)
  # Generate data Using doparalle, around 21secs
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  Z = foreach(i=1:n, .export = "gen_z_f1_probit") %dopar% {
    Z_i = rep(0,J)
    delta_i = rep(0, J)
    l_i = rep(0,J)
    p_i = rep(0,J)
    eta_i = rep(0,J)
    for(j in 1:J){
      ind_i = X[,length(gamma)+1]==i
      if(J==1){
        x = X[ind_i, 1:length(gamma)]
      }else{
        x = X[ind_i, 1:length(gamma)][j,]
      }
      #Bernoulli probability
      #p = exp(x %*% gamma + b_delta[i])/(1+exp(x %*% gamma + b_delta[i]))
      p = pnorm(x %*% gamma + b_delta[i])
      delta = rbinom(1, size = 1, p)
      eta = 0
      if(delta == 1){
        c_vec = c()
        linear = as.numeric(x %*% Alpha + b_eta[i] + rnorm(1))
        for(k in 1:K){
          if(thres_aug[k] < linear & thres_aug[k+1] > linear ){
            eta = k
          }
          # if(k==1){
          #   c_vec = c(c_vec, pnorm(thres[k] - x %*% Alpha - b_eta[i])) 
          # }else if(k==K){
          #   c_vec = c(c_vec, 1 - pnorm(thres[k-1] - x %*% Alpha - b_eta[i])) 
          # }else{
          #   c_vec = c(c_vec, pnorm(thres[k] - x %*% Alpha - b_eta[i]) - sum(c_vec) ) 
          # }
        }
        z_gen = gen_z_f1_probit(runif(1), eta = eta,byStep = byStep,phi.mat = phi.mat, prior = prior, end = 150)
        l_i[j] = linear
      }else{
        z_gen= 0
        l_i[j] = NA
      }
      eta_i[j] = eta
      Z_i[j] = z_gen
      delta_i[j] = delta
      p_i[j]=p
    }
    list(Z_i, delta_i, p_i, l_i,eta_i)
  }
  stopCluster(cl)
  Z1 = lapply(Z, function(x) do.call(cbind, x))
  frame_temp = do.call(rbind, Z1)
  Z = frame_temp[,1]
  delta = frame_temp[,2]
  p = frame_temp[,3]
  l = frame_temp[,4]
  Eta = frame_temp[,5]
  if(J==1){
    G=NULL
  }else{
    G = X[,length(gamma)+1]
  }
  X = X[,-(length(gamma)+1)] # make the X matrix with only data
  return(list(Z = Z,X = X,G = G,b_delta = b_delta, b_eta = b_eta, delta = delta, p_delta = p, l =l,
              Eta = Eta))
}

cal_prob_probit = function(x, Alpha, gamma, K,thres){
  # Function to calculate the global probability and the local probability for given x
  # x: a vector with the same length as gamma
  # Alpha: a matrix with columns for each covariate
  # gamma: a vector with the same length as x
  p_glob = c(pnorm(x %*% gamma)) # prob for non-zero, single value
  c_vec = c()
  for(k in 1:K){
    if(k==1){
      c_vec = c(c_vec, pnorm(thres[k] - x %*% Alpha )) 
    }else if(k==K){
      c_vec = c(c_vec, 1 - pnorm(thres[k-1] - x %*% Alpha )) 
    }else{
      c_vec = c(c_vec, pnorm(thres[k] - x %*% Alpha) - sum(c_vec)) 
    }
  }
  p_local = c_vec
  
  return(list(p_glob = p_glob, p_local = p_local))
}

cal_quantity_probit = function(quantity = "mean", x, Alpha, gamma, phi.mat, grid, byStep, K, thres = thres){
  # Function to calculate a quantity of the whole and local distribution
  # quantity: "mean", or certain percentage 0.5(median)
  # x: a vector with the same length as gamma
  # Alpha: a vector
  # gamma: a vector with the same length as x
  # phi.mat: output from Bdensity_raw$phi.mat or Bdensity_sim
  # grid: output from Bdensity_raw$grid
  
  #Calculate global and local probabablity
  prob_list = cal_prob_probit(x, Alpha, gamma, K, thres = thres)
  #Calculate the b-spline combined distribution
  phi.mat.model = phi.mat
  for(i in 1:ncol(phi.mat.model)){
    phi.mat.model[,i] = phi.mat[,i] * prob_list$p_local[i]
  }
  grid_prob = rowSums(phi.mat.model)
  #Calculate quantity
  if(is.null(quantity)){
    quantity_mean = sum(grid * grid_prob * byStep) * prob_list$p_glob
    if(0.5 > (1-prob_list$p_glob)){
      quantity_glob = 0.5 - (1-prob_list$p_glob)
      quantity_median = grid[which.min(abs(
        quantity_glob - cumsum(grid_prob*byStep*prob_list$p_glob)))] 
    }else{
      quantity_median = 0
    }
    if(0.25 > (1-prob_list$p_glob)){
      quantity_glob = 0.25 - (1-prob_list$p_glob)
      quantity_0.25 = grid[which.min(abs(
        quantity_glob - cumsum(grid_prob*byStep*prob_list$p_glob)))] 
    }else{
      quantity_0.25 = 0
    }
    if(0.75 > (1-prob_list$p_glob)){
      quantity_glob = 0.75 - (1-prob_list$p_glob)
      quantity_0.75 = grid[which.min(abs(
        quantity_glob - cumsum(grid_prob*byStep*prob_list$p_glob)))] 
    }else{
      quantity_0.75 = 0
    }
    return(list(mean = quantity_mean, median = quantity_median, 
                perc.25 = quantity_0.25, perc.75 = quantity_0.75))
     
  }else{
    if(quantity == "mean"){
      local_quantity = sum(grid * grid_prob * byStep)# mean for the b-spline part
      global_quantity = local_quantity * prob_list$p_glob# global mean
    }else if(is.numeric(quantity) & quantity<1){
      if(quantity > (1 - prob_list$p_glob)){
        local_quantity = grid[which.min(abs(quantity - cumsum(grid_prob*byStep) ))]
        quantity_glob = quantity - (1-prob_list$p_glob)
        global_quantity = grid[which.min(abs(
          quantity_glob - cumsum(grid_prob*byStep*prob_list$p_glob) 
        ))]
      }else{
        global_quantity = 0
        local_quantity = grid[which.min(abs(quantity - cumsum(grid_prob*byStep) ))]
      }
    }else{
      stop("quantity must be less than 1 or 'mean'")
    }
    return(list(local = local_quantity, global = global_quantity))
  }
}


coeff_effect_probit = function(quantity = "mean", Alpha, gamma, phi.mat, grid, byStep, local = F, K, thres){
  # Function to calculate the coefficient size for one unit increase in x
  # quantity: "mean", or certain percentage 0.5(median)
  # Alpha: a vector with each element for different covariate
  # gamma: a vector with the same length as x
  # phi.mat: output from Bdensity_raw$phi.mat or Bdensity_sim
  # grid: output from Bdensity_raw$grid
  M = length(gamma)
  if(is.null(quantity)){
    base_quantity = cal_quantity_probit(NULL, x = c(1, rep(0,M-1)),
                                 Alpha, gamma, phi.mat, grid, byStep, K = K, thres = thres)
    quantity_mean = c()
    quantity_median = c()
    quantity_0.25 = c()
    quantity_0.75 = c()
    for(i in 2:M){
      x = rep(0,M)
      x[1] = 1
      x[i] = 1
      q = cal_quantity_probit(NULL, x, Alpha, gamma, phi.mat, grid, byStep, K = K, thres = thres)
      quantity_mean = c(quantity_mean, q$mean)
      quantity_median = c(quantity_median, q$median)
      quantity_0.25 = c( quantity_0.25,q$perc.25)
      quantity_0.75 = c( quantity_0.75,q$perc.75)
    }
    effect_mean  = quantity_mean - base_quantity$mean
    effect_median = quantity_median - base_quantity$median
    effect_0.25 = quantity_0.25 - base_quantity$perc.25
    effect_0.75 = quantity_0.75 - base_quantity$perc.75
    return(list(effect_mean = c(base_quantity$mean, effect_mean),
                effect_median = c(base_quantity$median, effect_median),
                effect_0.25 = c(base_quantity$perc.25, effect_0.25),
                effect_0.75 = c(base_quantity$perc.75, effect_0.75)
                ))
  }else{
    base_quantity = cal_quantity_probit(quantity, x = c(1, rep(0,M-1)),
                                 Alpha, gamma, phi.mat, grid, byStep,K, thres = thres)
    quantity_global = c()
    quantity_local = c()
    for(i in 2:M){
      x = rep(0,M)
      x[1] = 1
      x[i] = 1
      q = cal_quantity_probit(quantity, x, Alpha, gamma, phi.mat, grid, byStep,K = K, thres = thres)
      quantity_global = c(quantity_global, q$global)
      quantity_local = c(quantity_local, q$local)
    }
    effect_global  = quantity_global - base_quantity$global
    effect_local = quantity_local - base_quantity$local
    return(list(
      effect_global = c(base_quantity$global, effect_global),
      effect_local = c(base_quantity$local, effect_local)
    ))
  }
  
}


gamma_b_spline_compare = function( theta, kappa, delta, phi.mat.model,byStep=0.01,
                                   grid = NULL,  model_prior, G_indx = NULL, 
                                   legend = T,xlab_T = T, ylab_T = T, title_T = T,multiple = T, CI = NULL){
  
  
  # the model density for b-splines with the model prob
  for(i in 1:ncol(phi.mat.model)){
    phi.mat.model[,i] = phi.mat.model[,i] * model_prior[i]
  }
  
  d_gamma = dggamma(grid, theta=theta, kappa=kappa, delta = delta, log = F)
  g = ggplot() + geom_line(aes(y = d_gamma, x = grid, color = "True")) + 
    geom_line(aes(y = rowSums(phi.mat.model), x = grid, color = "Model"))+
    geom_ribbon(aes(x = grid,ymin = CI[,1], ymax = CI[,2]), alpha = 0.3)+
    theme_bw() 
  
  if(legend == T){
    g = g + theme(legend.position = "bottom")+labs(x = "Observation values",
                                                   y = "Density",
                                                   title = "Population distribution for positive part comparing model to true",
                                                   color = "")
  }else{
    g= g+  theme(legend.position="none",
                 axis.title = element_blank())
  }
  return(g)
  
}



gamma_b_spline_compare_plot = function(x, alpha, delta_coef, k, est_Alpha, phi.mat.model,byStep = 0.01, 
                                       G = NULL, true_b_eta = NULL, MCMCfit, X = NULL,
                                       G_indx = NULL, last_n = 1000, est_thres){
  #x is the vector for the evaluating x vector
  # true_Alpha: a matrix of true alpha, row number equal to number of covariates, column number of b-spline
  #phi.mat: the b-spline generating basis 
  #byStep: 
  #G: the group index for each observation, if it is not null, then we draw 9 plots for 9 random subject based on their ture b_eta and model b_eta
  #true_b_eta: a matrix of #subject*(K-1)
  #MCMCfit: model fit
  #X: the matrix of covariates we used to fit the model
  #G_indx: the group index we are going to plot, recommend 4 or 9 indexes 
  k = k
  mu = x %*% as.vector(alpha)
  sigma = sqrt(exp(x %*% as.vector(delta_coef)))
  delta = k/sigma
  theta = exp(mu)/((1/k^2)^(sigma/k))
  kappa = 1/(sigma*k)
  
  K_model = dim(phi.mat.model)[2]
  CI = cal_CI_B_density_probit(MCMCfit, x, last_n = last_n)
  
  c_vec = c()
  for(k in 1:K_model){
    if(k==1){
      c_vec = c(c_vec, pnorm(est_thres[k] - x %*% est_Alpha ))
    }else if(k==K_model){
      c_vec = c(c_vec, 1 - pnorm(est_thres[k-1] - x %*% est_Alpha ))
    }else{
      c_vec = c(c_vec, pnorm(est_thres[k] - x %*% est_Alpha) - sum(c_vec) )
    }
  }
  prior_est = c_vec
  d_result = data.frame(parameters = paste("b_prior", 1:length(c_vec)), 
                        Model = c(prior_est ))
  return(list(gamma_b_spline_compare(theta, kappa, delta, phi.mat.model, byStep=byStep,
                               grid = MCMCfit$grid, prior_est, multiple = F, CI = CI), d_result) )
}