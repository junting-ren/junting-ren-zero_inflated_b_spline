library(splines)

#K-2 is the number of knots
# N is the number of observation Z
# Z is the observations: Binge drinking
# Gden:matrix(0,nrow=N,ncol=K), empty matrix with elements 0
#ref=1: this is the index for reference b-spline of Alpha coefficient
Bdensity<-function(byStep=0.001,mu=0,K,N,Z,Gden, ref=1, max_Z = NULL, knots = NULL, intercept = T){
  #max_Z is the true maximum of Z
  if(is.null(max_Z)){
    max_Z=max(abs(Z))
  }
  grid=seq(byStep,max_Z+byStep,by=byStep)
  if(K>=3){
    if(is.null(knots)){
      knots=seq(byStep,max_Z+byStep,len=K) #equal space;This includes starting and ending point
    }
    
    Basis.mat=matrix(0,nrow=length(grid),ncol=K);
    
    Basis.mat=bs(grid,knots=knots[2:(length(knots)-1)],degree=3,intercept=intercept,
                 Boundary.knots=c(knots[1],knots[length(knots)]))
    
    phi.mat=Basis.mat[,1:K]# Get rid of the last two basis since usually this density would be low at high value of z
    
    
    
    # With intercept
    for(k in 1:K){
      phi.mat[,k]= phi.mat[,k]/(sum(phi.mat[,k])*byStep)
    }
    phi.mat = phi.mat[,c(ref, c(1:K)[-ref] )]
    # phi.mat[,3:(K-2)]=phi.mat[,3:(K-2)]*4/(3*(knots[2]-knots[1]))
    # phi.mat[,1]=phi.mat[,1]*4/(knots[2]-knots[1])
    # phi.mat[,2]=phi.mat[,2]*2/(knots[2]-knots[1])
    # phi.mat[,(K-1)]=phi.mat[,(K-1)]/(knots[2]-knots[1])
    # phi.mat[,K]=phi.mat[,K]*4/(3*(knots[2]-knots[1]));
    #interpolate(get the density of z at different b-spline);
    for(i in 1:K){				
      Gden[,i]=approx(grid,phi.mat[,i],abs(Z),rule=2,ties="ordered")$y
    }
  }else{
    #no knot, only 4 b spline
    Basis.mat=bs(grid,degree=3,intercept=intercept)
    if(K!=4){#When K<4
      phi.mat=Basis.mat[,-4]
    }
    for(k in 1:ncol(phi.mat)){
      phi.mat[,k]= phi.mat[,k]/(sum(phi.mat[,k])*byStep)
    }
    for(i in 1:ncol(phi.mat)){				
      Gden[,i]=approx(grid,phi.mat[,i],abs(Z),rule=2,ties="ordered")$y
    }
  }
  return (Gden);
}


Bdensity_raw<-function(byStep=0.001,mu=0,K,N,Z,ref=1, max_Z = NULL, knots = NULL, intercept = T){
  if(is.null(max_Z)){
    max_Z=max(abs(Z))
  }
  grid=seq(byStep,max_Z+byStep,by=byStep)
  if(K>=3){
    if(is.null(knots)){
      knots=seq(byStep,max_Z+byStep,len=K) #equal space;This includes starting and ending point
    }
    
    
    Basis.mat=matrix(0,nrow=length(grid),ncol=K);
    
    Basis.mat=bs(grid,knots=knots[2:(length(knots)-1)],degree=3,intercept=intercept,
                 Boundary.knots=c(knots[1],knots[length(knots)]))
    
    phi.mat=Basis.mat[,1:K]# Get rid of the last two basis since usually this density would be low at high value of z
    # With intercept
    for(k in 1:K){
      phi.mat[,k]= phi.mat[,k]/(sum(phi.mat[,k])*byStep)
    }
    phi.mat = phi.mat[,c(ref, c(1:K)[-ref] )]
    # phi.mat[,3:(K-2)]=phi.mat[,3:(K-2)]*4/(3*(knots[2]-knots[1]))
    # phi.mat[,1]=phi.mat[,1]*4/(knots[2]-knots[1])
    # phi.mat[,2]=phi.mat[,2]*2/(knots[2]-knots[1])
    # phi.mat[,(K-1)]=phi.mat[,(K-1)]/(knots[2]-knots[1])
    # phi.mat[,K]=phi.mat[,K]*4/(3*(knots[2]-knots[1]));
    #interpolate(get the density of z at different b-spline);
  }else{
    #no knot, only 3 b spline
    Basis.mat=bs(grid,degree=3,intercept=intercept)
    if(K!=4){
      phi.mat=Basis.mat[,-4]
    }
    for(k in 1:ncol(phi.mat)){
      phi.mat[,k]= phi.mat[,k]/(sum(phi.mat[,k])*byStep)
    }
  }
  return(list(phi.mat = phi.mat, grid = grid, knots = knots));
}
