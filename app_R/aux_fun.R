times <- function(lambda,df) {   
  time1=log(df$psd0[1]/df$psd0[2]*df$ct[2]/df$ct[1])/-(lambda[1]-lambda[2])
  time2=log(df$psd0[2]/df$psd0[3]*df$ct[3]/df$ct[2])/-(lambda[2]-lambda[3])
  time3=log(df$psd0[1]/df$psd0[3]*df$ct[3]/df$ct[1])/-(lambda[1]-lambda[3])
  return(c(time1,time2,time3))
}

SSt <- function(lambda,df) {   
  time = times(lambda,df)
  res=(time-mean(time))^2
  if (mean(time)>0)
    return(sum(res))
  else 
    return(sum(res)+1e2)
}

## backward simulation
psd0 <- function(lambda,time,cinput){
  c0=rep(NaN,length(lambda))
  for(i in 1:length(lambda)){
    c0[i]=cinput[i]*exp((-lambda[i])*time)
  }
  csum=sum(c0)
  return(c(c0/csum,csum))
}


## single backward simulation
single_pred <- function(tmax,dt,lambda,filt){
  time = -seq(0,tmax,dt)
  t=tibble()
  for (i in 1:length(filt)){
    ln_norm_c = lambda[2*i-1] + (min(lambda[2*i],0))*time
    # lncpred2 = lambda[3] + (lambda[4])*time
    min(lambda[2*i],0)
    norm_c = exp(ln_norm_c)
    # cpred2 = exp(lncpred2)
    Filter= rep(filt[i],length(time))
    t = rbind(t,tibble(time,ln_norm_c,norm_c,Filter))
  }
  aa = t%>%select(time,norm_c,Filter)%>%
    pivot_wider(names_from = Filter, values_from = norm_c)%>%
    mutate(csum = rowSums(select(.,-time)))%>%
    mutate(psd = select(.,-c(time,csum))/csum)
  # pivot_longer(select(.,psd),everything())
  # mutate(psd2 = select(.,"1.2")/csum)  
  t=cbind(t,aa$psd%>%pivot_longer(everything(),values_to = 'frac')%>%
            arrange(name)%>%select(frac))
  
  return(t)
}

## estimating lambda 3 based on 2 samples
lambdac <- function(cs2,cs1,Nsizes,Nsim,lambda,dsamples,vstream){
  ts = dsamples/vstream
  lbd1 = lambda[1,]
  lbd2 = lambda[2,]
  ct = sum(cs1)
  
  c1=rep(NaN,Nsim)
  c2=rep(NaN,Nsim)
  c3=rep(NaN,Nsim)
  
  for (i in 1:Nsim){
    c1[i] = cs2[1]*exp(-lbd1[i]*ts)
    c2[i] = cs2[2]*exp(-lbd2[i]*ts)
    c3[i] = ct - c1[i] - c2[i]
  }
  
  lambdac = -1/ts * log(c3/cs2[3])
  return(lambdac)
}


## Mode of a distribution
mode<- function(x){
 bin_size <- (max(x) - min(x))/15  
 histogram <- hist(x, breaks = seq(min(x), max(x), by = bin_size), include.lowest = TRUE,plot = FALSE)
 max_bin <- which.max(histogram$counts)
 mode <- mean(histogram$breaks[max_bin:(max_bin+1)])
return(mode)
}


MC_new <- function(table,Nsim,dt,tmax,Nsizes){
  set.seed(1)
  #names of the filters
  filt =  unique(table$Filter)
  # tmax = tmax
  # dt = dt
  # 
  dfMC=tibble()
  # withProgress(message = 'Calculation in progress', {
  for (i in 1:Nsim){ 
  lambda =  runif(2*Nsizes,min=table$conf.low,max=table$conf.high)
  t = single_pred(tmax,dt,lambda,filt)
  t$Nsim = i
  dfMC = rbind(dfMC,t)
  incProgress(1/Nsim)
  }
  # })  
  return(dfMC)
}

## Monte Carlo simulation 
MC <- function(table,Nsim,dt,tmax,Nsizes,csample){
  set.seed(1)
  df = table%>%
    select(Filter,estimate,conf.low,conf.high)
  # Nsim = input$nsim 
  # dt = input$dt
  time = seq(0,tmax,dt) 
  lambda_low = df$conf.low
  lambda_high = df$conf.high
  
  lambda=matrix(NaN,nrow=Nsizes,ncol=Nsim)
  f_names=rep(NaN,Nsizes)
  lbd_names=rep(NaN,Nsizes)
  psd0_st=matrix(NaN,nrow=Nsim*length(time),ncol=2*Nsizes+3)
  
  for(i in 1:Nsizes){
    lambda[i,]= runif(Nsim,min=df$conf.low[i],max=df$conf.high[i])
    f_names[i] = paste0("f0", i)
    lbd_names[i] = paste0("lbd", i)
  }
  
  for (j in 1:length(time)){
    for (i in 1:Nsim){
      psd0_st[i+(j-1)*Nsim,1] = time[j]
      psd0_st[i+(j-1)*Nsim,2:(Nsizes+1)] = lambda[1:Nsizes,i]
      psd0_st[i+(j-1)*Nsim,(Nsizes+2):(2*Nsizes+2)] = 
        psd0(lambda[,i],time[j],csample)
      psd0_st[i+(j-1)*Nsim,(2*Nsizes+3)] = i
    }
  }
  
  dfMC = as_tibble(psd0_st)
  names(dfMC) = c('time',lbd_names,f_names,'C0','Nsim')
  return(dfMC)
}

## Monte Carlo simulation with partition
MC_p <- function(table,Nsim,dt,tmax,Nsizes,csample){
  set.seed(1)
  df = table%>%
    select(Filter,estimate,conf.low,conf.high)

  time = seq(0,tmax,dt) 
  lambda_low = df$conf.low
  lambda_high = df$conf.high
  
  lambda=matrix(NaN,nrow=Nsizes,ncol=Nsim)
  partit=matrix(NaN,nrow=Nsizes,ncol=Nsim)
  f_names=rep(NaN,Nsizes)
  f2_names=rep(NaN,Nsizes)
  lbd_names=rep(NaN,Nsizes)
  p_names=rep(NaN,Nsizes)
  psd0_st=matrix(NaN,nrow=Nsim*length(time),ncol=4*Nsizes+4)
  
  for(i in 1:Nsizes){
    partit[i,]= runif(Nsim,min=0,max=1)
    lambda[i,]= runif(Nsim,min=df$conf.low[i],max=df$conf.high[i])
    f_names[i] = paste0("f0", i)
    f2_names[i] = paste0("f20", i)
    lbd_names[i] = paste0("lbd", i)
    p_names[i]= paste0("p", i)
  }
  
  for (j in 1:length(time)){
    for (i in 1:Nsim){
      psd0_st[i+(j-1)*Nsim,1] = time[j]
      psd0_st[i+(j-1)*Nsim,2:(Nsizes+1)] = lambda[1:Nsizes,i]
      psd0_st[i+(j-1)*Nsim,(Nsizes+2):(2*Nsizes+1)] = partit[1:Nsizes,i]
      psd0_st[i+(j-1)*Nsim,(2*Nsizes+2):(3*Nsizes+2)] = 
        psd0(lambda[,i],time[j],csample*partit[,i])
      
      psd0_st[i+(j-1)*Nsim,(3*Nsizes+3):(4*Nsizes+3)] = 
        psd0(lambda[,i],time[j],csample*(1-partit[,i]))
      
      psd0_st[i+(j-1)*Nsim,(4*Nsizes+4)] = i
    }
  }
  
  dfMC = as_tibble(psd0_st)
  names(dfMC) = c('time',lbd_names,p_names,f_names,'C0',f2_names,'C20','Nsim')
  return(dfMC)
}


## Monte Carlo simulation with lambda 3 comlementary
MC_compl <- function(table,Nsim,dt,tmax,Nsizes,csample,csample1,dsamples,vstream){
  set.seed(1)
  df = table%>%
    select(Filter,estimate,conf.low,conf.high)
  time = seq(1,tmax,dt) 
  # lambda_low = df$conf.low
  # lambda_high = df$conf.high
  
  lambda=matrix(NaN,nrow=Nsizes,ncol=Nsim)
  f_names=rep(NaN,Nsizes)
  lbd_names=rep(NaN,Nsizes)
  psd0_st=matrix(NaN,nrow=Nsim*length(time),ncol=2*Nsizes+3)
  
  for(i in 1:(Nsizes-1)){
    lambda[i,]= runif(Nsim,min=df$conf.low[i],max=df$conf.high[i])
  }
  lambda[Nsizes,] = lambdac(csample,csample1,Nsizes,Nsim,lambda,dsamples,vstream)
  
  for(i in 1:Nsizes){
    # lambda[i,]= runif(Nsim,min=df$conf.low[i],max=df$conf.high[i])
    f_names[i] = paste0("f0", i)
    lbd_names[i] = paste0("lbd", i)
  }
  
  for (j in 1:length(time)){
    for (i in 1:Nsim){
      psd0_st[i+(j-1)*Nsim,1] = time[j]
      psd0_st[i+(j-1)*Nsim,2:(Nsizes+1)] = lambda[1:Nsizes,i]
      psd0_st[i+(j-1)*Nsim,(Nsizes+2):(2*Nsizes+2)] =
        psd0(lambda[,i],time[j],csample)
      psd0_st[i+(j-1)*Nsim,(2*Nsizes+3)] = i
    }
  }
  
  dfMC = as_tibble(psd0_st)
  names(dfMC) = c('time',lbd_names,f_names,'C0','Nsim')
  return(dfMC)
}



# SMIM functions ----------------------------------------------------------

k_func = function(rs,rh,U,D,lambda, beta,logT1,logT2)
{
  t1 = 10^logT1
  t2 = 10^logT2
  Theta = rs +
    lambda*(1-exp(rh*t1)*gammainc(-beta,t1*(1/t2+rh))*(1+t2*rh)^beta/
              gammainc(-beta,t1/t2))
  k = (-U+sqrt(U^2+4*D*Theta))/(2*D);
  return (k)
}

C_func = function(xlim,C0,rs,rh,U,D,lambda, beta,logT1,logT2)
{  
  x=c(seq(0,xlim,length.out = 1000))
  C = C0*exp(-k_func(rs,rh,U,D,lambda, beta,logT1,logT2)*x)
  df = tibble(x,C)
  return (df)
}

C_loc = function(xloc,C0,rs,rh,U,D,lambda, beta,logT1,logT2)
{  
  C = C0*exp(-k_func(rs,rh,U,D,lambda, beta,logT1,logT2)*xloc)
  return (C)
}
