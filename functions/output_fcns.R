#############################################################################

# Functions for output tidy results from zero-inflated semi-parametric model

#############################################################################

# Function for output ALpha, Gamma, sigma_delta and sigma_eta posterior mean and credibility interval
meanCI <- function(model, level = 0.05, last_n = 1000, true_para = NULL, return_all = T){
  # model: a zero-inflated b-spline model returned result
  #level: the type 1 error rate
  #last_n=1000: take the last_n number of iteration for mean and credibility interval
  # return_all = T: return all the nuisance parameters as well
  #A list containing the true parameters: gamma, 
  zero_inflate = model$zero_inflate
  if(is.null(model$col_names)){
    covariate_names = c("intercept")
    for(m in 1:(length(GMN)-1)){
      covariate_names= c(covariate_names, paste("covariate", m))
    }
  }else{
    covariate_names = model$col_names
  }
  if(return_all == T){
    #Gamma
    if(zero_inflate){
      GMN = apply(simplify2array(tail(model$gamma_array, last_n)), 1, mean)
      GMN_upper = apply(simplify2array(tail(model$gamma_array, last_n)), 1, quantile, probs = 1 - level/2)
      GMN_lower = apply(simplify2array(tail(model$gamma_array, last_n)), 1, quantile, probs = level/2)
      Gamma_frame = data.frame(variable = covariate_names, Gamma_mean = GMN, 
                               Gamma_up = GMN_upper, Gamma_low = GMN_lower)
    }else{
      Gamma_frame = NULL
    }
    
    # Random effect variance
    if(model$grouped & zero_inflate){
      # sigma_delta
      sigma_b_delta_sq_mean = mean(unlist(model$sigma_b_delta_sq))
      sigma_b_delta_sq_upper = quantile(unlist(model$sigma_b_delta_sq), 1 - level/2)
      sigma_b_delta_sq_lower = quantile(unlist(model$sigma_b_delta_sq), level/2)
      # sigma_eta
      sigma_b_eta_sq_mean = apply(simplify2array(tail(model$Sigma_b_eta, last_n)), 1:2, mean)[1,1]
      sigma_b_eta_sq_upper = apply(simplify2array(tail(model$Sigma_b_eta, last_n)), 1:2, quantile, probs = 1- level/2)[1,1]
      sigma_b_eta_sq_lower = apply(simplify2array(tail(model$Sigma_b_eta, last_n)), 1:2, quantile, probs = level/2)[1,1]
      sigma_frame = data.frame(variable = c("sigma_b_delta_sq", "sigma_b_eta_sq"),
                               mean = c(sigma_b_delta_sq_mean, sigma_b_eta_sq_mean), 
                               up = c(sigma_b_delta_sq_upper, sigma_b_eta_sq_upper),
                               low = c(sigma_b_delta_sq_lower, sigma_b_eta_sq_lower), row.names = NULL)
    }else if(model$grouped & !zero_inflate){
      # sigma_eta
      sigma_b_eta_sq_mean = apply(simplify2array(tail(model$Sigma_b_eta, last_n)), 1:2, mean)[1,1]
      sigma_b_eta_sq_upper = apply(simplify2array(tail(model$Sigma_b_eta, last_n)), 1:2, quantile, probs = 1- level/2)[1,1]
      sigma_b_eta_sq_lower = apply(simplify2array(tail(model$Sigma_b_eta, last_n)), 1:2, quantile, probs = level/2)[1,1]
      sigma_frame = data.frame(variable = c( "sigma_b_eta_sq"),
                               mean = c(sigma_b_eta_sq_mean), 
                               up = c(sigma_b_eta_sq_upper),
                               low = c(sigma_b_eta_sq_lower), row.names = NULL)
    }else{
      sigma_frame = NULL
    }
    
    # ALpha
    AMN = apply(simplify2array(tail(model$alpha_array, last_n)), 1:2, mean)
    AMN_upper = apply(simplify2array(tail(model$alpha_array, last_n)), 1:2, quantile, probs = 1 - level/2)
    AMN_lower = apply(simplify2array(tail(model$alpha_array, last_n)), 1:2, quantile, probs = level/2)
    
    Alpha_frame = data.frame(matrix(nrow = dim(AMN)[1] * (dim(AMN)[2]-1), ncol = 4))
    names(Alpha_frame) = c("covariate_names",  "Alpha mean", "Alpha_up", "Alpha_low")
    if(is.null(model$col_names)){
      for(m in 1:dim(AMN)[1]){
        for(k in 2:dim(AMN)[2]){
          if(m==1){
            Alpha_frame[(m-1)*(dim(AMN)[2]-1)+(k-1),1] = paste("intercept", "_", "spline", k, sep = "")
          }else{
            Alpha_frame[(m-1)*(dim(AMN)[2]-1)+(k-1),1] = paste("covarate",m,"_","spline", k, sep = "")
          }
          Alpha_frame[(m-1)*(dim(AMN)[2]-1)+(k-1),2] = AMN[m,k]
          Alpha_frame[(m-1)*(dim(AMN)[2]-1)+(k-1),3] = AMN_upper[m,k]
          Alpha_frame[(m-1)*(dim(AMN)[2]-1)+(k-1),4] = AMN_lower[m,k]
        }
      }
    }else{
      for(m in 1:dim(AMN)[1]){
        for(k in 2:dim(AMN)[2]){
          Alpha_frame[(m-1)*(dim(AMN)[2]-1)+(k-1),1] = paste(model$col_names[m],"_","spline", k, sep = "")
          Alpha_frame[(m-1)*(dim(AMN)[2]-1)+(k-1),2] = AMN[m,k]
          Alpha_frame[(m-1)*(dim(AMN)[2]-1)+(k-1),3] = AMN_upper[m,k]
          Alpha_frame[(m-1)*(dim(AMN)[2]-1)+(k-1),4] = AMN_lower[m,k]
        }
      }
    }
    
    mean_effect = apply(simplify2array(tail(model$effect_mean, last_n)), 1, mean)
    mean_effect_upper = apply(simplify2array(tail(model$effect_mean, last_n)), 1, quantile, probs = 1 - level/2)
    mean_effect_lower = apply(simplify2array(tail(model$effect_mean, last_n)), 1, quantile, probs = level/2)
    mean_effect_frame = data.frame(variable = covariate_names, 
                                   mean_effect = mean_effect, 
                                   mean_up = mean_effect_upper, 
                                   mean_low = mean_effect_lower)
    
    
    median_effect = apply(simplify2array(tail(model$effect_median, last_n)), 1, mean)
    median_effect_upper = apply(simplify2array(tail(model$effect_median, last_n)), 1, quantile, probs = 1 - level/2)
    median_effect_lower = apply(simplify2array(tail(model$effect_median, last_n)), 1, quantile, probs = level/2)
    median_effect_frame = data.frame(variable = covariate_names, 
                                     median_effect = median_effect, 
                                     median_up = median_effect_upper, 
                                     median_low = median_effect_lower)
    
    perc.25_effect = apply(simplify2array(tail(model$effect_0.25, last_n)), 1, mean)
    perc.25_effect_upper = apply(simplify2array(tail(model$effect_0.25, last_n)), 1, quantile, probs = 1 - level/2)
    perc.25_effect_lower = apply(simplify2array(tail(model$effect_0.25, last_n)), 1, quantile, probs = level/2)
    perc.25_effect_frame = data.frame(variable = covariate_names, 
                                      perc.25_effect = perc.25_effect, 
                                      perc.25_up = perc.25_effect_upper, 
                                      perc.25_low = perc.25_effect_lower)
    
    
    perc.75_effect = apply(simplify2array(tail(model$effect_0.75, last_n)), 1, mean)
    perc.75_effect_upper = apply(simplify2array(tail(model$effect_0.75, last_n)), 1, quantile, probs = 1 - level/2)
    perc.75_effect_lower = apply(simplify2array(tail(model$effect_0.75, last_n)), 1, quantile, probs = level/2)
    perc.75_effect_frame = data.frame(variable = covariate_names, 
                                      perc.75_effect = perc.75_effect, 
                                      perc.75_up = perc.75_effect_upper, 
                                      perc.75_low = perc.75_effect_lower)
    
    dist_effect = apply(simplify2array(tail(model$effect_dist, last_n)), 1, mean)
    dist_effect_upper = apply(simplify2array(tail(model$effect_dist, last_n)), 1, quantile, probs = 1 - level/2)
    dist_effect_lower = apply(simplify2array(tail(model$effect_dist, last_n)), 1, quantile, probs = level/2)
    dist_effect_frame = data.frame(variable = c("mean","median", "perc.25", "perc.75"), 
                                   dist_effect = dist_effect, 
                                   dist_up = dist_effect_upper, 
                                   dist_low = dist_effect_lower)
    
    if(!is.null(true_para)){
      Gamma_frame = cbind(Gamma_frame, true = true_para$gamma)
      mean_effect_frame = cbind(mean_effect_frame, true= true_para$effect_mean)
      median_effect_frame = cbind(median_effect_frame, true = true_para$effect_median)
      mean_effect_frame$contain = mean_effect_frame$true <= mean_effect_frame$mean_up & 
        mean_effect_frame$true >= mean_effect_frame$mean_low
      median_effect_frame$contain = median_effect_frame$true <= median_effect_frame$median_up & 
        median_effect_frame$true >= median_effect_frame$median_low
      
      perc.25_effect_frame = cbind(perc.25_effect_frame, true= true_para$effect_0.25)
      perc.75_effect_frame = cbind(perc.75_effect_frame, true = true_para$effect_0.75)
      perc.25_effect_frame$contain = perc.25_effect_frame$true <= perc.25_effect_frame$perc.25_up & 
        perc.25_effect_frame$true >= perc.25_effect_frame$perc.25_low
      perc.75_effect_frame$contain = perc.75_effect_frame$true <= perc.75_effect_frame$perc.75_up & 
        perc.75_effect_frame$true >= perc.75_effect_frame$perc.75_low
      
      dist_effect_frame = cbind(dist_effect_frame, true = true_para$effect_dist)
      dist_effect_frame$contain = dist_effect_frame$true <= dist_effect_frame$dist_up & 
        dist_effect_frame$true >= dist_effect_frame$dist_low
      
      true_alpha = c()
      for(m in 1: length(GMN)){
        true_alpha = c(true_alpha, true_para$alpha[m,-1])
      }
      if(nrow(Alpha_frame)==length(true_alpha)){
        Alpha_frame = cbind(Alpha_frame, true = true_alpha)
        Alpha_frame$contain = Alpha_frame$true <= Alpha_frame$Alpha_up & 
          Alpha_frame$true >= Alpha_frame$Alpha_low
      }
      if(model$grouped){
        sigma_frame = cbind(sigma_frame, true = c(true_para$sigma_b_delta^2, true_para$sigma_b_eta^2))
        sigma_frame$contain = sigma_frame$true <= sigma_frame$up & 
          sigma_frame$true >= sigma_frame$low
      }
      
      Gamma_frame$contain = Gamma_frame$true <= Gamma_frame$Gamma_up & 
        Gamma_frame$true >= Gamma_frame$Gamma_low
      
      
    }
    
    return(list(gamma_frame = Gamma_frame, alpha_frame = Alpha_frame, sigma_frame = sigma_frame,
                mean_effect_frame = mean_effect_frame, median_effect_frame=median_effect_frame,
                perc.25_effect_frame=perc.25_effect_frame, perc.75_effect_frame=perc.75_effect_frame,
                dist_effect_frame = dist_effect_frame))
  }else{# only return the distribution effects
    mean_effect = apply(simplify2array(tail(model$effect_mean, last_n)), 1, mean)
    mean_effect_upper = apply(simplify2array(tail(model$effect_mean, last_n)), 1, quantile, probs = 1 - level/2)
    mean_effect_lower = apply(simplify2array(tail(model$effect_mean, last_n)), 1, quantile, probs = level/2)
    mean_effect_frame = data.frame(variable = covariate_names, 
                                   mean_effect = mean_effect, 
                                   mean_up = mean_effect_upper, 
                                   mean_low = mean_effect_lower)
    
    
    median_effect = apply(simplify2array(tail(model$effect_median, last_n)), 1, mean)
    median_effect_upper = apply(simplify2array(tail(model$effect_median, last_n)), 1, quantile, probs = 1 - level/2)
    median_effect_lower = apply(simplify2array(tail(model$effect_median, last_n)), 1, quantile, probs = level/2)
    median_effect_frame = data.frame(variable = covariate_names, 
                                     median_effect = median_effect, 
                                     median_up = median_effect_upper, 
                                     median_low = median_effect_lower)
    
    perc.25_effect = apply(simplify2array(tail(model$effect_0.25, last_n)), 1, mean)
    perc.25_effect_upper = apply(simplify2array(tail(model$effect_0.25, last_n)), 1, quantile, probs = 1 - level/2)
    perc.25_effect_lower = apply(simplify2array(tail(model$effect_0.25, last_n)), 1, quantile, probs = level/2)
    perc.25_effect_frame = data.frame(variable = covariate_names, 
                                      perc.25_effect = perc.25_effect, 
                                      perc.25_up = perc.25_effect_upper, 
                                      perc.25_low = perc.25_effect_lower)
    
    
    perc.75_effect = apply(simplify2array(tail(model$effect_0.75, last_n)), 1, mean)
    perc.75_effect_upper = apply(simplify2array(tail(model$effect_0.75, last_n)), 1, quantile, probs = 1 - level/2)
    perc.75_effect_lower = apply(simplify2array(tail(model$effect_0.75, last_n)), 1, quantile, probs = level/2)
    perc.75_effect_frame = data.frame(variable = covariate_names, 
                                      perc.75_effect = perc.75_effect, 
                                      perc.75_up = perc.75_effect_upper, 
                                      perc.75_low = perc.75_effect_lower)
    
    dist_effect = apply(simplify2array(tail(model$effect_dist, last_n)), 1, mean)
    dist_effect_upper = apply(simplify2array(tail(model$effect_dist, last_n)), 1, quantile, probs = 1 - level/2)
    dist_effect_lower = apply(simplify2array(tail(model$effect_dist, last_n)), 1, quantile, probs = level/2)
    dist_effect_frame = data.frame(variable = c("mean","median", "perc.25", "perc.75"), 
                                   dist_effect = dist_effect, 
                                   dist_up = dist_effect_upper, 
                                   dist_low = dist_effect_lower)
    if(!is.null(true_para)){
      mean_effect_frame = cbind(mean_effect_frame, true= true_para$effect_mean)
      median_effect_frame = cbind(median_effect_frame, true = true_para$effect_median)
      mean_effect_frame$contain = mean_effect_frame$true <= mean_effect_frame$mean_up & 
        mean_effect_frame$true >= mean_effect_frame$mean_low
      median_effect_frame$contain = median_effect_frame$true <= median_effect_frame$median_up & 
        median_effect_frame$true >= median_effect_frame$median_low
      
      
      perc.75_effect_frame = cbind(perc.75_effect_frame, true = true_para$effect_0.75)
      perc.75_effect_frame$contain = perc.75_effect_frame$true <= perc.75_effect_frame$perc.75_up & 
        perc.75_effect_frame$true >= perc.75_effect_frame$perc.75_low
    
    }
    return(list(mean_effect_frame = mean_effect_frame, median_effect_frame=median_effect_frame, perc.75_effect_frame=perc.75_effect_frame))
  }
  
}


convergence_plot <- function(model , save = T){
  zero_inflate = model$zero_inflate
  M = dim(model$alpha_array[[1]])[1]
  K = dim(model$alpha_array[[1]])[2]
  if(save==T){
    pdf(file = "Convergence plots.pdf")
  }
  if(model$grouped){
    g_indx = sample(as.numeric(as.factor(model$grouping_order)),2) #the subjects for visualizing mean random effect
  }
  
  if(zero_inflate==T){
    if(is.null(model$col_names)){
      mcmc_Gamma = mcmc(t(sapply(model$gamma_array, cbind)))
      plot(mcmc_Gamma, ylab = "Gamma")
    }else{
      mcmc_Gamma = mcmc(t(sapply(model$gamma_array, cbind)))
      plot(mcmc_Gamma, ylab = "Gamma", main = paste(model$col_names))
    }
  }
  
  for(i in 1:M){
    mcmc_Alpha = mcmc(t(sapply(model$alpha_array,function(x)x[i,2:K])))
    plot(mcmc_Alpha, ylab = paste("Alpha", model$col_names[i],sep = " "))
  }
  if(zero_inflate==T & model$grouped){
    mcmc_sigma_delta= mcmc(unlist(model$sigma_b_delta_sq))
    plot(mcmc_sigma_delta, ylab = "sigma_delta_sq")
    for(i in g_indx){
      mcmc_b_delta = mcmc(sapply(model$b_delta, function(x)x[i]))
      plot(mcmc_b_delta, main = paste("b_delta",i), ylab = paste("b_delta","subject", i))
    }
  }
  if(model$grouped){
    mcmc_Sigma_eta= mcmc(sapply(model$Sigma_b_eta, function(x)x[1,1]))
    plot(mcmc_Sigma_eta, ylab = "Sigma_eta_sq")
    for(i in g_indx){
      mcmc_b_eta = mcmc(t(sapply(model$b_eta, function(x)x[i,])))
      plot(mcmc_b_eta, main = paste("b_eta",i), ylab = paste("b_eta","subject", i))
    }
  }
  
  if(is.null(model$col_names)){
    mcmc_mean = mcmc(t(sapply(model$effect_mean, cbind)))
    plot(mcmc_global, ylab = "Mean effect")
  }else{
    mcmc_mean = mcmc(t(sapply(model$effect_mean, cbind)))
    plot(mcmc_mean, ylab = "Mean effect", main = paste(model$col_names))
  }
  if(is.null(model$col_names)){
    mcmc_median = mcmc(t(sapply(model$effect_median, cbind)))
    plot(mcmc_median, ylab = "Median effect")
  }else{
    mcmc_median = mcmc(t(sapply(model$effect_median, cbind)))
    plot(mcmc_median, ylab = "Median effect", main = paste(model$col_names))
  }
  if(is.null(model$col_names)){
    mcmc_perc0.25 = mcmc(t(sapply(model$effect_0.25, cbind)))
    plot(mcmc_perc0.25, ylab = "perc 0.25 effect")
  }else{
    mcmc_perc0.25 = mcmc(t(sapply(model$effect_0.25, cbind)))
    plot(mcmc_perc0.25, ylab = "perc 0.25 effect", main = paste(model$col_names))
  }
  if(is.null(model$col_names)){
    mcmc_perc0.75 = mcmc(t(sapply(model$effect_0.75, cbind)))
    plot(mcmc_perc0.75, ylab = "perc 0.75 effect")
  }else{
    mcmc_perc0.75 = mcmc(t(sapply(model$effect_0.75, cbind)))
    plot(mcmc_perc0.75, ylab = "perc 0.75 effect", main = paste(model$col_names))
  }
  
  mcmc_dist = mcmc(t(sapply(model$effect_dist, cbind)))
  plot(mcmc_dist, ylab = "Distribution quantity")
  
  if(save == T){
    dev.off()
    return("Plots saved to working directory")
  }
  
}


# Draw the corresponding b-spline basis from the algorithm
draw_b_spline = function(model, draw_mix = F,prior=c(0.5,0.1,0.1,0.1,0.2)){
  require(ggplot2)
  grid=model$grid
  n_spline = ncol(model$Bden_raw)
  if(draw_mix == F){
    phi.mat = data.frame(cbind(model$Bden_raw, grid))
    names(phi.mat) = c(paste("spline", 1:n_spline, sep = ""), "grid")
    local_frame = NULL
    for(c in paste("spline", 1:n_spline, sep = "")){
      min_g = min(grid[phi.mat[,c]>0])
      max_g = max(grid[phi.mat[,c]>0])
      local_frame = rbind(local_frame, data.frame(spline = c, start = min_g, end = max_g))
    }
    phi.mat_long = phi.mat %>% pivot_longer(1:all_of(n_spline), names_to = "b_spline", values_to = "value")
    # Plotting of individual spline density
    p = phi.mat_long %>% 
    ggplot(aes(x = grid, y = value, color = b_spline)) + 
      geom_line() + labs(y = "Density", title = "B-Spline Density", x = "Support") +
      theme_bw()+theme(legend.position="bottom")
      
  }else{
    phi.mat = model$Bden_raw
    for(i in 1:ncol(phi.mat)){
      phi.mat[,i] = phi.mat[,i] * prior[i]
    }
    plot_frame = data.frame(cbind(density = rowSums(phi.mat), grid = grid))
    p = plot_frame %>% ggplot(aes(x = grid, y = density)) + geom_line() + 
      labs(title = "Combined b-spline density")+theme_bw()+
      theme(legend.position="bottom")
  }
  return(list(local_frame, p))
}

plot_compare_2x = function(x1, x2, model, last_n = 1000,
                           include_zero = F){
  # generate phi.mat
  grid=model$grid
  phi.mat=model$Bden_raw
  
  # Placeholder for filtering the data 
  
  # Calculate mean Gamma and Alpha
  GMN = apply(simplify2array(tail(model$gamma_array, last_n)), 1, mean)
  AMN = apply(simplify2array(tail(model$alpha_array, last_n)), 1:2, mean)
  # Calculate global p and local ps
  denom_m1 = 0
  for(col in 1:ncol(AMN)){
    denom_m1 = denom_m1 + exp(x1 %*% AMN[,col])
  }
  model_b_prior1 = c()
  for(col in 1:ncol(AMN)){
    model_b_prior1 = c(model_b_prior1, exp(x1 %*% AMN[,col])/denom_m1)
  }
  
  denom_m2 = 0
  for(col in 1:ncol(AMN)){
    denom_m2 = denom_m2 + exp(x2 %*% AMN[,col])
  }
  model_b_prior2 = c()
  for(col in 1:ncol(AMN)){
    model_b_prior2 = c(model_b_prior2, exp(x2 %*% AMN[,col])/denom_m2)
  }
  # p0 = 1/(exp(x %*% GMN) + 1)
  # model_glob_prior = c(p0,1-p0)
  
  # the model density for b-splines with the model prob
  phi.mat.model1 = phi.mat
  for(i in 1:ncol(phi.mat.model1)){
    phi.mat.model1[,i] = phi.mat[,i] * model_b_prior1[i]
  }
  phi.mat.model2 = phi.mat
  for(i in 1:ncol(phi.mat.model2)){
    phi.mat.model2[,i] = phi.mat[,i] * model_b_prior2[i]
  }
  
  # getting the model density
  CI1 = cal_CI_B_density(model, x1, last_n = last_n, g = NULL)
  CI2 = cal_CI_B_density(model, x2, last_n = last_n, g = NULL)
  density1 = rowSums(phi.mat.model1) 
  density2 = rowSums(phi.mat.model2) 
  p =  ggplot() + geom_line(aes(y = density1, x = grid, color = "Baseline")) + 
    geom_ribbon(aes(x = grid,ymin = CI1[,1], ymax = CI1[,2]), alpha = 0.3)+
    geom_line(aes(y = density2, x = grid, color = "age standardized = 3")) + 
    geom_ribbon(aes(x = grid,ymin = CI2[,1], ymax = CI2[,2]), alpha = 0.3)+
    labs(x = "Observation values",
         y = "Density",
         color = "Model") +
    theme_bw()+theme(legend.position="bottom")+
    scale_color_manual(values=c("#CC6666", "#9999CC"))
    #scale_linetype_identity()
    #scale_linetype_manual(values = c("Baseline" = "dashed","age standardized = 3" = "solid"))
  
  return(p)
}


plot_checkfit_average = function(Z, X, x = NULL, model, last_n = 1000,
                                 include_zero = F){
  # generate phi.mat
  grid=model$grid
  phi.mat=model$Bden_raw
  
  # Placeholder for filtering the data 
  
  #If x is null, then we use the whole data and use x = mean(x) as fitting point
  if(is.null(x)){
    x = c()
    for(j in 1:ncol(X)){
      x = c(x, mean(X[,j]))
    }
    x = round(x, 3)
  }
  
  # Calculate mean Gamma and Alpha
  GMN = apply(simplify2array(tail(model$gamma_array, last_n)), 1, mean)
  AMN = apply(simplify2array(tail(model$alpha_array, last_n)), 1:2, mean)
  # Calculate global p and local ps
  denom_m = 0
  for(col in 1:ncol(AMN)){
    denom_m = denom_m + exp(x %*% AMN[,col])
  }
  model_b_prior = c()
  for(col in 1:ncol(AMN)){
    model_b_prior = c(model_b_prior, exp(x %*% AMN[,col])/denom_m)
  }
  p0 = 1/(exp(x %*% GMN) + 1)
  model_glob_prior = c(p0,1-p0)
  
  # the model density for b-splines with the model prob
  phi.mat.model = phi.mat
  for(i in 1:ncol(phi.mat.model)){
    phi.mat.model[,i] = phi.mat[,i] * model_b_prior[i]
  }
  if(include_zero == T){
    grid = c(0,grid)
    # getting the model density
    bin = round(max(grid),1)
    density = rowSums(phi.mat.model) * model_glob_prior[2] #b density 
    density =  c(model_glob_prior[1], density) # model density after adding the 0 probability mass
    p =  ggplot() + geom_line(aes(y = density[-1], x = grid[-1], color = "positive")) + 
      geom_point(aes(y = density[1], x = grid[1], color = "zero"))+
      geom_histogram(aes(x = Z, y = ..density.., fill = "data"), bins = bin, alpha = 0.8)+
      scale_fill_manual(values = c("data"= "gray"))+
      labs(x = "Observation values",
           y = "Density",
           color = "Model",
           fill = "") +theme_classic()+
      theme_bw()+theme(legend.position="bottom")
  }else{
    # getting the model density
    CI = cal_CI_B_density(model, x, last_n = last_n, g = NULL)
    bin = round(max(grid),1)/3
    density = rowSums(phi.mat.model) 
    Z = Z[Z!=0]
    p =  ggplot() + geom_line(aes(y = density, x = grid, color = "positive")) + 
      geom_histogram(aes(x = Z, y = ..density.., fill = "data"), bins = bin, alpha = 0.3)+
      scale_fill_manual(values = c("data"= "yellow"))+
      geom_ribbon(aes(x = grid,ymin = CI[,1], ymax = CI[,2]), alpha = 0.3)+
      labs(x = "Observation values",
           y = "Density",
           color = "Model",
           fill = "") +theme_classic()+
      theme_bw()+theme(legend.position="bottom")
  }
  
  return(list(p, spline_prob = model_b_prior, global_prob = model_glob_prior))
}



plot_checkfit_individual = function(Z, X, G, x = NULL, model, last_n = 1000, G_indx = NULL){
  require(gridExtra)
  # generate phi.mat
  grid=model$grid
  phi.mat=model$Bden_raw
  
  # Placeholder for filtering the data 
  
  #If x is null, then we use the whole data and use x = mean(x) as fitting point
  x_list = list() # a list for the mean x for each subject
  z_list = list() # a list for the count for each subject
  if(is.null(G_indx)){
    G_indx = sample(unique(G[duplicated(G)]),4)
    for(i in 1:length(G_indx)){
      x = colMeans(X[G==G_indx[i],]) 
      x_list[[i]] =  x 
      z = Z[G == G_indx[i]]
      z_list[[i]]= z
    }
    G_indx = c(1:length(unique(G)))[sort(unique(G)) %in% G_indx] 
  }else{ 
    for(i in 1:length(G_indx)){
      x = colMeans(X[G==G_indx[i],]) 
      x_list[[i]] =  x 
      z = Z[G == G_indx[i]]
      z_list[[i]]= z
    }
    G_indx = c(1:length(unique(G)))[sort(unique(G)) %in% G_indx]
  }
  # Calculate mean b_eta and b_delta
  b_eta_list = list()
  b_delta_list = list()
  for(i in 1:length(G_indx)){
    b_eta_model = rowMeans(sapply(tail(model$b_eta, last_n), function(x) x[G_indx[i],])) 
    b_delta_model = mean(sapply(tail(model$b_delta, last_n), function(x) x[G_indx[i]]))
    b_eta_list[[i]] =  b_eta_model
    b_delta_list[[i]] =  b_delta_model
  }
  # Calculate mean Gamma and Alpha
  GMN = apply(simplify2array(tail(model$gamma_array, last_n)), 1, mean)
  AMN = apply(simplify2array(tail(model$alpha_array, last_n)), 1:2, mean)
  # Calculate global p and local ps
  plot_list = list()
  model_b_prior_list = list()
  model_glob_prior_list = list()
  grid = c(0, grid)
  plot_list = lapply(1:length(x_list), function(i){
    x = x_list[[i]]
    b_eta = c(0,b_eta_list[[i]])
    b_delta = b_delta_list[[i]]
    z = z_list[[i]]
    denom_m = 0
    for(col in 1:ncol(AMN)){
      denom_m = denom_m + exp(x %*% AMN[,col] + b_eta[col])
    }
    model_b_prior = c()
    for(col in 1:ncol(AMN)){
      model_b_prior = c(model_b_prior, exp(x %*% AMN[,col] + b_eta[col])/denom_m)
    }
    model_b_prior_list[[i]]<<- model_b_prior
    p0 = 1/(exp(x %*% GMN + b_delta) + 1)
    model_glob_prior = c(p0,1-p0)
    model_glob_prior_list[[i]]<<- model_glob_prior
    # the model density for b-splines with the model prob
    phi.mat.model = phi.mat
    for(j in 1:ncol(phi.mat.model)){
      phi.mat.model[,j] = phi.mat[,j] * model_b_prior[j]
    }
    # getting the model density
    density = rowSums(phi.mat.model) * model_glob_prior[2] #b density 
    density =  c(model_glob_prior[1], density) # model density after adding the 0 probability mass
    bin = round(max(grid),1)
    ggplot() + geom_line(aes(y = density[-1], x = grid[-1], color = "positive")) + 
      geom_point(aes(y = density[1], x = grid[1], color = "zero"))+
      geom_histogram(aes(x = z, y = ..density.., fill = "data"), bins = bin, alpha = 0.8)+
      scale_fill_manual(values = c("data"= "gray"))+
      labs(title = paste("Subject", paste(G_indx[i])),
           fill = "",
           colour = "") +
      theme_bw()+
      theme(legend.position="none",
            axis.title = element_blank())
      
  
  })
  p1 = arrangeGrob(grobs = plot_list[1:4], ncol = 2)
  return(p1)
  #return(list(p1, spline_prob = model_b_prior_list, global_prob = model_glob_prior_list, z = z_list))
}

plot_b_spline = function(phi.mat){
  plot(phi.mat[,1]~grid, ylim=c(0,max(phi.mat)), type='l', lwd=2, col=1, 
     xlab="Cubic B-spline basis", ylab="")
for (j in 2:ncol(phi.mat)) lines(phi.mat[,j]~grid, lwd=2, col=j)
}
