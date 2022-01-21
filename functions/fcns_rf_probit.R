rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}


Draw_Gamma_log_M = function(y, X, G, SG, Phi,b_delta,same_length_G){
  T_1 = list()
  T_2 = list()
  V = list()
  i = 1
  for(g in unique(G)){
    indx = G == g
    n = sum(indx)
    X_i = X[indx,]
    y_i = y[indx]
    #V_i = diag(n)+rep(1,n)%*%t(rep(1,n))*sigma_b_delta_sq
    V_i = diag(n)
    V[[i]]=V_i
    V_i_inv = solve(V_i)
    if(n>1){
      T_1[[i]] = t(X_i)%*%V_i_inv%*%X_i
      T_2[[i]] = t(X_i)%*%V_i_inv%*%(y_i-b_delta[g])
      
    }else{
      T_1[[i]] = X_i%*%V_i_inv%*%t(X_i)
      T_2[[i]] = X_i%*%V_i_inv%*%(y_i-b_delta[g])
    }
    i = i+1
  }
  
  var = solve(solve(SG)+Reduce("+",T_1))
  mean = var%*%Reduce("+",T_2)
  Gamma = MASS::mvrnorm(n=1,mu = mean,Sigma=var)
  return(list(par=Gamma, accp=1, V = V))
}


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




Draw_Alpha_probit = function(l1, X1, G1,Sigma_alpha, b_eta){
  T_1 = list()
  T_2 = list()
  V = list()
  i = 1
  
  for(g in unique(G1)){
    indx = G1 == g
    n = sum(indx)
    X_i = X1[indx,]
    l_i = l1[indx]
    #V_i = diag(n)+t(t(rep(1,n)))%*%Sigma_b_eta%*%t(rep(1,n))
    V_i = diag(n)
    V[[i]]=V_i
    V_i_inv = solve(V_i)
    if(n>1){
      T_1[[i]] = t(X_i)%*%V_i_inv%*%X_i
      T_2[[i]] = t(X_i)%*%V_i_inv%*%(l_i-b_eta[g])
    }else{
      T_1[[i]] = X_i%*%V_i_inv%*%t(X_i)
      T_2[[i]] = X_i%*%V_i_inv%*%(l_i-b_eta[g])
    }
    i = i+1
  }

  var = solve(solve(Sigma_alpha)+Reduce("+",T_1))
  mean = var%*%Reduce("+",T_2)
  Alpha = MASS::mvrnorm(n=1,mu = mean,Sigma=var)
  return(list(par=Alpha, accp=1, V = V))
}


Draw_b_eta_i_probit = function(g, l1, G1, X1, Alpha, Sigma_b_eta)
{
  indx = G1==g
  J = sum(indx)
  l_i = l1[indx]
  X_i = X1[indx,]
  v_square = Sigma_b_eta/(1+J*Sigma_b_eta)
  mu = v_square * sum(l_i - X_i%*%Alpha)
  b_eta_i_new = rnorm(1, mean = mu, sd = sqrt(v_square))
  accp = 1
  return(list(par=b_eta_i_new, accp=accp))
} 

log_thres_p = function(thres_aug, Eta1, G1, X1, Alpha, b_eta){
  XA = X1%*%Alpha
  log_p = sum(log(pnorm(thres_aug[Eta1+1] - XA - b_eta[G1])-
                    pnorm(thres_aug[Eta1] - XA - b_eta[G1])))
  return(log_p)
}
log_thres_p_repar = function(theta, Eta1, G1, X1, Alpha, b_eta,K,Sigma_theta){
  # theta is one to one to thres hold, K-2 length since thres[1] is 0 fixed.
  thres = rep(0, K-1)
  for(j in 2:(K-1)){
    thres[j] = sum(exp(theta[1:(j-1)]))
  }
  thres_aug = c(-Inf, thres, Inf)
  XA = X1%*%Alpha
  temp =  (pnorm(thres_aug[Eta1+1] - XA - b_eta[G1])-
             pnorm(thres_aug[Eta1] - XA - b_eta[G1]))
  log_p = sum(log( temp
                       ), -1/2*t(theta) %*% Sigma_theta %*%theta )
  return(log_p)
}

log_p_b_eta_probit = function(b_eta_i, i, G, X_i, thres, Eta_i,
                              Alpha, Sigma_b_eta){
  thres_aug = c(-Inf,thres, Inf)
  thres_up = thres_aug[Eta_i+1]
  thres_low = thres_aug[Eta_i]
  log_p = -0.5 * b_eta_i^2/Sigma_b_eta  + 
    sum(
      log(pnorm(thres_up - X_i%*%Alpha-b_eta_i)-pnorm(thres_low - X_i%*%Alpha-b_eta_i))
    )
  return(log_p)
}

Draw_b_eta_i_MH_probit = function(b_eta_i, i, G, X, thres, Eta = Eta,Z,
                        Alpha, Sigma_b_eta, K = 5, iter,burnIn, individual = F,
                        cov = diag(0.5, nrow = K-1)){
  #browser()
  X_i = X[G==i & Z > 0,]
  Eta_i = Eta[G==i & Z>0]
  if(length(Eta_i)==0){
    return(list(par=rnorm(n = 1, mean = 0,sd = sqrt(Sigma_b_eta)), accp=1))
  }

  b_eta_i_star = rnorm(n = 1, mean = b_eta_i,sd = sqrt(cov))
  
  
  log_p_b_eta_i_star = log_p_b_eta_probit(b_eta_i_star, i, G, X_i, thres,Eta_i,
                                   Alpha, Sigma_b_eta)
  log_p_b_eta_i_old = log_p_b_eta_probit(b_eta_i, i, G, X_i, thres,Eta_i,
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


Draw_b_delta_eta_i_probit_wishart = function(i, G,G1, X,X1, Gamma, Alpha, Sigma_b, Phi, yl_i)
{
  indx = G==i
  J = sum(indx)
  indx1 = G1==i
  J1 = sum(indx1)
  y_i = yl_i[1:J]
  X_i = X[indx,]
  X_ib = matrix(cbind(rep(1,J),rep(0, J)), ncol = 2)
  if(J1<1){
    V = solve(solve(Sigma_b)+t(X_ib)%*%X_ib)
    mu = V%*%(t(X_ib)%*%y_i-t(X_ib)%*%X_i%*%Gamma)
    MASS::mvrnorm(n=1,mu = mu,Sigma=V)
    b_new = MASS::mvrnorm(n=1,mu = mu,Sigma=V)
  }else{
    l1_i = yl_i[(J+1):length(yl_i)]
    X_1i = X1[indx1,]
    X_1ib = matrix(cbind(rep(0, J1), rep(1,J1)), ncol = 2)
    V = solve(solve(Sigma_b)+t(X_ib)%*%X_ib+t(X_1ib)%*%X_1ib)
    mu = V%*%(t(X_ib)%*%y_i-t(X_ib)%*%X_i%*%Gamma+t(X_1ib)%*%l1_i-t(X_1ib)%*%X_1i%*%Alpha)
    MASS::mvrnorm(n=1,mu = mu,Sigma=V)
    b_new = MASS::mvrnorm(n=1,mu = mu,Sigma=V)
  }
  
  accp = 1
  return(list(par=b_new, accp=accp))
} 


Draw_b_delta_eta_i_probit = function(i, G,G1, X,X1, Gamma, Alpha, sigma_b_delta_sq, Phi, y=y, l1 =l1 )
{
  indx = G==i
  J = sum(indx)
  indx1 = G1==i
  J1 = sum(indx1)
  l1_i = l1[indx1]
  X_1i = X1[indx1,]
  y_i = y[indx]
  X_i = X[indx,]
  if(J1<2){
    v_square = sigma_b_delta_sq/(1+(J)*sigma_b_delta_sq)
    mu = v_square * sum(y_i - X_i%*%gamma)
    b_delta_i_new = rnorm(1, mean = mu, sd = sqrt(v_square))
  }else{
    v_square = sigma_b_delta_sq/(1+(J+J1)*sigma_b_delta_sq)
    mu = v_square * (sum(y_i - X_i%*%gamma)+sum(l1_i-X_1i%*%Alpha))
    b_delta_i_new = rnorm(1, mean = mu, sd = sqrt(v_square))
  }
  
  accp = 1
  return(list(par=b_delta_i_new, accp=accp))
} 
