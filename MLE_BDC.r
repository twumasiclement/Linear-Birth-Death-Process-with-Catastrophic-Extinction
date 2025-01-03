
library(gmp) #Dealing with larger factorial limits with arbitrary precision
library(compiler)# byte code compilation
library(data.table)
library(parallel)
library(JuliaCall) #for julia commands & scripts

library("maxLik")#for maximum likelihood estimation/optimization
#library('latex2exp') #For adding LaTeX symbols to R plots

setwd("/home/clement/Clement_R_Results_folder")

#options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

# Calculate the number of cores for parLapply
#no_cores <- min(16, detectCores())

# Initiate cluster
#cluster <- makeCluster(no_cores)
#cluster
nc=min(16, detectCores())
nc



beta_n_j<-function(nmax){
    #calculate beta^n_j=gamma^n_j/n! for n=1,2 ,...,nmax and j=1,2,...,n
    bnj<-matrix(0,nmax,nmax)
    bnj[,1]<-1
    if(nmax==1) return(bnj)
    for(j in 2:nmax){
        for(n in j:nmax){
            bnj[n,j]<-bnj[n-1,j-1]/n +(n+j-1)*bnj[n-1,j]/n
        }
    }
    return(bnj[nmax,])
}

beta_n_j_compiler<-cmpfun(beta_n_j)




BDCconsts <- function(lambda, mu, rho,t) {
# Constants used in calculating distribution of BDC process at time t
  
  roots <- sort(Re(polyroot(c(mu, -(lambda+mu+rho), lambda))))
  v0 <- roots[1]
  v1 <- roots[2]
  sigma <- exp(-lambda*(v1 - v0)*t)
  k1 <- v0*v1*(1 - sigma)/(v1 - sigma*v0)
  k2 <- (v1*sigma - v0)/(v1 - sigma*v0)
  k3 <- (1 - sigma)/(v1 - sigma*v0)
  return(list(k1=k1, k2=k2, k3=k3, sigma=sigma, v0=v0, v1=v1))
}



#Function for the probability generating function G(z,t)
PGF_z=function(lambda,mu,rho,t,z,m){
   #v0<-((lambda+mu+rho)-sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
   #v1<-((lambda+mu+rho)+sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
    constants=BDCconsts(lambda,mu,rho,t)
    v0<- constants$v0
    v1<- constants$v1
    sigma<- constants$sigma
    num<-(v0*v1*(1-sigma))+(z*(v1*sigma-v0))
    den<- v1-(sigma*v0)-(z*(1-sigma))
    return( (num/den)^m)
}
PGF_z_compiler=cmpfun(PGF_z)
#PGF_z_compiler(lambda=0.513,mu=0.35,rho=0.003,t=1,z=1,m=2)


#Estimating C(t)=P(catastrophe resulting in 0 population)
Constant_added=function(lambda,mu,rho,t,z=1,m){
constant=1-PGF_z_compiler(lambda=lambda,mu=mu,rho=rho,t=t,z=z,m=m)
return(constant)
    }



Pn_m_t<-function(n,m,lambda, mu, rho,t){ 
 constants=BDCconsts(lambda=lambda, mu=mu, rho=rho,t=t)
    C_t=Constant_added(lambda=lambda, mu=mu, rho=rho,t=t,m=m)
    k1=constants$k1;k2=constants$k2;k3=constants$k3
    #if(n==Inf|m==Inf){return(0)}
    
    if(n== -1  && m>=1){return(C_t)} #P(X(t) = 0 and fish dead | X(0) = m) = C(t): 0 due to catastrophic event
    if(n== -1  && m==-1){return(1)} #P(X(t) = 0 and fish dead | X(0) = -1) = 1: 0 due to catastrophic event
    if(n== -1  && m==0){return(0)} #given parasites die out, the prob of catastrophe is 0.
    if(n>0  && m==-1){return(0)} #given catastrophe occurs, the prob of parasites>0 at any time t#given parasites die out, prob of catastrophe is 0. is 0.
    
    if(n==0 && m>=1) {return(k1^m)} #P( X(t) = 0 and fish alive | X(0) = m) = k1^m: 0 due to natural death
 
    if(n==0 && m==0){return(1)} #given 0 parasites the probability of having 0 parasites at any time t is 1
   
    if(n>0 && m==0){return(0)}  #given 0 parasites the probability of having  parasites n>0 at any time t is 0
   #Gamma is the coefficent for each term for the sums for any n
    ##############   FOR m<=n  #####################
    else if (n>0 && m<=n) {
    
   #sum_initializer=0
   # for (j in 1:m){
   # terms=beta_n_j(nmax=n)[1:j][j] *(k3^(n-j)) *  (k1^(m-j)) * ((k2+(k1*k3))^(j)) * (factorialZ(m)/(factorialZ(m-j)))
   # sum_initializer=sum_initializer+terms
     #   }
        
        #To detect the function global variables: "beta_n_j" & "factorialZ"
    #clusterExport(cluster,list("beta_n_j","factorialZ"),envir=globalenv())  
   
 terms<- lapply(1:m,function(j)       
 beta_n_j_compiler(nmax=n)[1:j][j] *(k3^(n-j)) *  (k1^(m-j)) * ((k2+(k1*k3))^(j)) * 
                  (factorialZ(m)/(factorialZ(m-j))))
       
    
    return(as.numeric(do.call("sum",terms)))
         
    }  
    
     ##############   FOR m>=n  #####################
    else if (n>0 && m>=n){
        
    #sum_initializer=0
    #for (j in 1:n){
    #terms=beta_n_j(nmax=n)[j] * (k3^(n-j)) *  (k1^(m-j)) * ((k2+(k1*k3))^(j)) * (factorialZ(m)/factorialZ(m-j))
    #sum_initializer=sum_initializer+terms
        #}
     
        terms<- lapply(1:n,function(j)
           
beta_n_j_compiler(nmax=n)[j]* (k3^(n-j)) *  (k1^(m-j)) * ((k2+(k1*k3))^(j)) * 
                        (factorialZ(m)/factorialZ(m-j)))
        return(as.numeric(do.call("sum",terms)))
        
   
    }
    
}



Pn_m_t_compiler=cmpfun(Pn_m_t)

#Pn_m_t_compiler(n=20,m=10,lambda=0.513, mu=0.35, rho=0.003,t=2)





# Maximum Likelikood Function (R)
#Obtaining MLE estimates given data across the 9 observed time points
MLE_BDC<-function(Parasite_numbers,X0=2,t=c(1,3,5,7,9,11,13, 15,17)){
   
    #LogLikelihood function to maximize
    BDC_Loglik=function(params){
    
   total_fish=dim(Parasite_numbers)[2] #total number of fish
        #The three parameters to be optimized
        lambda<-params[1]
        mu<-params[2]
        rho<-params[3]
   
    #n_t=colSums(Parasite_numbers)
    n_t<-Parasite_numbers
    n0=X0 #initial parasite number 
        
     log_prob<- lapply(1:total_fish,function(k)   
    c(log(Pn_m_t_compiler(n=n_t[[k]][1],m=n0,lambda=lambda, mu=mu, rho=rho,t=t[1])),
    log(Pn_m_t_compiler(n=n_t[[k]][2],m=n_t[[k]][1],lambda=lambda, mu=mu, rho=rho,t=t[2])),
   log(Pn_m_t_compiler(n=n_t[[k]][3],m=n_t[[k]][2],lambda=lambda, mu=mu, rho=rho,t=t[3])),
    log(Pn_m_t_compiler(n=n_t[[k]][4],m=n_t[[k]][3],lambda=lambda, mu=mu, rho=rho,t=t[4])),
    log(Pn_m_t_compiler(n=n_t[[k]][5],m=n_t[[k]][4],lambda=lambda, mu=mu, rho=rho,t=t[5])),
   log(Pn_m_t_compiler(n=n_t[[k]][6],m=n_t[[k]][5],lambda=lambda, mu=mu, rho=rho,t=t[6])),
   log(Pn_m_t_compiler(n=n_t[[k]][7],m=n_t[[k]][6],lambda=lambda, mu=mu, rho=rho,t=t[7])),
    log(Pn_m_t_compiler(n=n_t[[k]][8],m=n_t[[k]][7],lambda=lambda, mu=mu, rho=rho,t=t[8])),
    log(Pn_m_t_compiler(n=n_t[[k]][9],m=n_t[[k]][8],lambda=lambda, mu=mu, rho=rho,t=t[9]))
      )) 
        
         return(do.call("sum",log_prob))           
    }

BDC_Loglik_compiler=cmpfun(BDC_Loglik)
            
## Inequality constraints: lambda>0, mu>0, rho>0
    
A<-matrix(c(1, 0, 0,
              0, 1, 0,
              0, 0, 1), 3, 3, byrow=TRUE)
#B<-c(0,0,0.0001) #each >0 works
B<-c(0,0,0) #each >0
ineqCon <- list(ineqA=A, ineqB=B)

estimates<-maxLik(BDC_Loglik_compiler,start=c(lambda=2, mu=1,rho= 0.001),
                  constraints=ineqCon,method="NM") 

#return(print(summary(estimates))) #to print estimates & its standard errors and p-values 
      return(as.vector(estimates$estimate))#returning a vector of the estimates: birth, death & catastrophe
          }

MLE_BDC_compiler<- cmpfun(MLE_BDC) 





# MLE based on Julia's Log-Likelihood function
#Loading LogLikehood function from Julia
julia_command('BDCdir = "/home/clement/Clement_R_Results_folder"') # location of BDCfit.jl
julia_command("push!(LOAD_PATH, BDCdir)")
julia_command("using BDCfit")

MLE_BDC_Julia<-function(Parasite_numbers){
  
  #LogLikelihood function to maximize
  BDC_Loglik_Julia=function(params){
    
    #The three parameters to be optimized
    lambda<-params[1]
    mu<-params[2]
    rho<-params[3]
    
log_like_Julia<- julia_call("logL", lambda, mu, rho, t(Parasite_numbers)) 
    
  }
  
  BDC_Loglik_Julia_compiler=cmpfun(BDC_Loglik_Julia)
  
  ## Inequality constraints: lambda>0, mu>0, rho>0
  
  A<-matrix(c(1, 0, 0,
              0, 1, 0,
              0, 0, 1), 3, 3, byrow=TRUE)
  #B<-c(0,0,0.0001) #each >0 works
  B<-c(0,0,0) #each >0
  ineqCon <- list(ineqA=A, ineqB=B)
  
  estimates<-maxLik(BDC_Loglik_Julia_compiler,start=c(lambda=2, mu=1,rho= 0.001),
                    constraints=ineqCon,method="NM") 
  
  #return(print(summary(estimates))) #to print estimates & its standard errors and p-values 
  return(as.vector(estimates$estimate))#returning a vector of the estimates: birth, death & catastrophe
}


MLE_BDC_Julia_compiler<- cmpfun(MLE_BDC_Julia)





# MLE fitting using simulated data
#Function to determine & update events
SSA_update_event=function(X,b, d, c){
    
    fish_status <- 1 # fish starts out alive
    
    rate<-numeric(3)
    #probability of birth
    rate[1]<- b*X
    #probability of death
    rate[2]<- d*X
    #probability of catastrophe
    rate[3]<- c*X
    
    total_rate<- abs(rate[1]+rate[2]+rate[3])

     if (total_rate == 0) {
    return(list(X = X, t_incr = Inf)) # zero population
  }

    #Determine event occurence from single draw
    u<-runif(1,0, total_rate)
        if (u<abs(rate[1])){
             #birth of parasites
           X<-X+1   
        } else if(u<abs(rate[1]+rate[2])){
              #death of parasites
           X<-X-1  
       }  else {
           #catastrophe or death of fish
            X<-0
            fish_status<-2
    }
        t_incr <- rexp(1, total_rate) # time increment
   
  return(list(X = X, t_incr = t_incr,fish_status=fish_status))
    
}
SSA_update_event_compiler=cmpfun(SSA_update_event)

#Function for simulation based on the update function

#Parameters- (Birth, death and catastrophic events)
#parameters=c(0.512, 0.35, 0.003) 



#Function for exact simulation
Exact_BDC<-function(X0,b,d,c,ti=0,tmax=30){

pop_max <- 10000 # stop the simulation if total population exceeds this limit

#Time variable
#ti<- 0; tmax=30;
#Main loop
     save_ti <- as.vector(1:tmax)    #Times to simulate
    
    save_TF <- rep(FALSE, length(save_ti))
    pop <-matrix(-1,1,length(save_ti)) # parasite pop at each location (rows) and time point (cols) 
    alive <- rep(2, length(save_ti)) # host fish status at each time point
    alive_ti <- 1 #fish starts out alive
     X<-X0
    pop_ti <- sum(X)
    
    
    while(sum(save_TF) < length(save_ti)){# Calculate rate of events
  if (sum(pop_ti) > pop_max) {
      cat("pop max exceeded","\n")
    break
  } 
  if (alive_ti == 2) {
      break
    }
        
 
      
    output <-SSA_update_event_compiler(X,b, d, c)
    #Update time to next event
    ti <- ti + output$t_incr
  
    
    #if (X < 0) break #break if there is negative population
   
    # Events to occur
    save_new <- which((ti >= save_ti) & !save_TF)
    for (i in save_new){
    pop[,i]<- pop_ti
    alive[i] <- alive_ti
    }
      save_TF <- (ti >= save_ti)
      X<- output$X
     pop_ti <- sum(X)
        
    alive_ti <- output$fish_status
    
      
  
}
    
    return(list(pop=pop,alive = alive))
}

Exact_BDC_compiler=cmpfun(Exact_BDC)

#Print output from the Simulation Model
#Exact_BDC_compiler(X0=2,b=parameters[1],d=parameters[2],c=parameters[3],ti=0,tmax=17)





ExactSSA_sim <- function(X0,b,d,c,ti=0,tmax=30, nrep=1e1, nc=min(16, detectCores())) {


  f <- function(x) Exact_BDC_compiler(X0=X0,b=b,d=d,c=c,ti=ti,tmax=tmax)
  x <- simplify2array(mclapply(1:nrep, f, mc.cores=nc))
  #print(x[,,]) #output combined as a matrix with dimension simulation output length X nrep 
    pop=x["pop" ,]
    #alive=x["alive" ,]
    #return(list(pop=pop,alive=alive))
      return(pop)
}




#Log_likelihood based on R codes
BDC_Loglik=function(params,Parasite_numbers,X0=2,t=c(1,3,5,7,9,11,13, 15,17)){
    
   total_fish=dim(Parasite_numbers)[2] #total number of fish
        #The three parameters to be optimized
        lambda<-params[1]
        mu<-params[2]
        rho<-params[3]
   
    #n_t=colSums(Parasite_numbers)
    n_t<-Parasite_numbers
    n0=X0 #initial parasite number 
        
     log_prob<- mclapply(1:total_fish,function(k)   
    c(log(Pn_m_t_compiler(n=n_t[[k]][1],m=n0,lambda=lambda, mu=mu, rho=rho,t=t[1])),
    log(Pn_m_t_compiler(n=n_t[[k]][2],m=n_t[[k]][1],lambda=lambda, mu=mu, rho=rho,t=t[2])),
   log(Pn_m_t_compiler(n=n_t[[k]][3],m=n_t[[k]][2],lambda=lambda, mu=mu, rho=rho,t=t[3])),
    log(Pn_m_t_compiler(n=n_t[[k]][4],m=n_t[[k]][3],lambda=lambda, mu=mu, rho=rho,t=t[4])),
    log(Pn_m_t_compiler(n=n_t[[k]][5],m=n_t[[k]][4],lambda=lambda, mu=mu, rho=rho,t=t[5])),
   log(Pn_m_t_compiler(n=n_t[[k]][6],m=n_t[[k]][5],lambda=lambda, mu=mu, rho=rho,t=t[6])),
   log(Pn_m_t_compiler(n=n_t[[k]][7],m=n_t[[k]][6],lambda=lambda, mu=mu, rho=rho,t=t[7])),
    log(Pn_m_t_compiler(n=n_t[[k]][8],m=n_t[[k]][7],lambda=lambda, mu=mu, rho=rho,t=t[8])),
    log(Pn_m_t_compiler(n=n_t[[k]][9],m=n_t[[k]][8],lambda=lambda, mu=mu, rho=rho,t=t[9]))
      ),mc.cores=nc
       ) 
        
         return(do.call("sum",log_prob))           
    }
BDC_Loglik_compiler=cmpfun(BDC_Loglik)

#d1=c(2,3,5,6,7,8,2,5,6)
#d2=c(2,10,10,6,9,8,2,2,8)
#data=data.frame(d1,d2)

#est=BDC_Loglik_compiler(params=c(0.5,.3,.1),Parasite_numbers=data)
#est


#Simulated 100 realization of data of 50 fish with MLE
#Using both R and Julia log-likelihood functions

Simul_data_nruns_MLE<- function(X0,b,d,c,ti=0,tmax=30, nrep=1e1,nruns=2, nc=min(16, detectCores())) {
  
  Sim_data<-array(dim = c(nruns,9,nrep))
  Simulation_func <- function(x) ExactSSA_sim(X0=X0,b=b,d=d,c=c,ti=ti,tmax=tmax,nrep=nrep)
  x1_out <- simplify2array(mclapply(1:nruns, Simulation_func, mc.cores=nc))#convert output as an array
  
  
  for (k in 1:nruns){
    for(i in 1:nrep){
      # saving results for each 100 realizations and 
      Sim_data[k, ,i]<- x1_out[,k ][[i]][seq(1,17,by=2)]      
    }
  }
  
  #Computing MLE based on R function
  MLE_R <- function(k) MLE_BDC(Parasite_numbers=as.data.frame(Sim_data[k, ,]))
  MLE_estimates <- mclapply(1:nruns, MLE_R, mc.cores=nc)  #as lists
  
  LogLike_MLE<-as.data.frame(do.call("rbind",MLE_estimates)) 
  
  #Computing MLE based on Julia function 
  MLE_estimates_Julia <- mclapply(1:nruns, 
                                  function(k) MLE_BDC_Julia_compiler(Parasite_numbers=as.data.frame(Sim_data[k, ,])), 
                                  mc.cores=nc)  #as lists
  LogLike_MLE_Julia <- as.data.frame(do.call("rbind",MLE_estimates_Julia)) 
  
  #Computing the Log-likelihood based on MLE estimates (R)
  Loglik<-mclapply(1:nruns, function(k)
    BDC_Loglik_compiler(params=c(LogLike_MLE[[1]][k],LogLike_MLE[[2]][k],LogLike_MLE[[3]][k]),
                        Parasite_numbers=as.data.frame(Sim_data[k, ,])), mc.cores=nc)
  
  #Computing the Log-likelihood based on true parameters values (R)   
  LogLik_trueValue<- mclapply(1:nruns, function(k)
    BDC_Loglik_compiler(params=c(b,d,c),Parasite_numbers=as.data.frame(Sim_data[k, ,])), mc.cores=nc)
  
  
  #Computing the Log-likelihood based on true parameters values (Julia)   
  LogLik_trueValue_Julia<- mclapply(1:nruns, function(k)
    julia_call("logL", b,d,c, t(as.data.frame(Sim_data[k, ,]))), mc.cores=nc)
  
  
  
  #Computing the Log-likelihood based on MLE estimates from Julia function                           
  Loglike_Julia<-mclapply(1:nruns, function(k)
    julia_call("logL",  LogLike_MLE_Julia[[1]][k], LogLike_MLE_Julia[[2]][k],LogLike_MLE_Julia[[3]][k], 
               t(as.data.frame(Sim_data[k, ,]))), mc.cores=nc)                           
  
  
  #Computing means of mle of birth,death and catastrophe across all replicates 
  #Computing the loglike at these means
  MLE_means<- as.vector(apply(LogLike_MLE[ ,1:3],2,mean))
  
  
  Loglik_MLE_means<-mclapply(1:nruns, function(k)
    BDC_Loglik_compiler(params=c(MLE_means[1],MLE_means[2],MLE_means[3]),
                        Parasite_numbers=as.data.frame(Sim_data[k, ,])), mc.cores=nc)
  
  #Computing means of mle of birth,death and catastrophe across all replicates from Julia function
  #Computing the loglike at these means from Julia
  MLE_means_Julia<- as.vector(apply(LogLike_MLE_Julia[ ,1:3],2,mean))                         
  
  Loglik_MLE_means_Julia<-mclapply(1:nruns, function(k)
    julia_call("logL", MLE_means_Julia[1],MLE_means_Julia[2],MLE_means_Julia[3], 
               t(as.data.frame(Sim_data[k, ,]))), mc.cores=nc)                            
  
  
  #Returning final output for R (First part)
  LogLike_MLE$birth_MLE_mean<-rep(MLE_means[1],length=nruns)
  LogLike_MLE$death_MLE_mean<-rep(MLE_means[2],length=nruns)
  LogLike_MLE$catastrophe_MLE_mean<-rep(MLE_means[3],length=nruns)                           
  LogLike_MLE$Loglik_MLE<-unlist(Loglik)
  LogLike_MLE$LogLik_MLE_means<- unlist(Loglik_MLE_means) 
  LogLike_MLE$Loglik_True_values <-unlist(LogLik_trueValue)
  
  
  names(LogLike_MLE)<-c("b_MLE_R","d_MLE_R","c_MLE_R","b_MLE_mean_R","d_MLE_mean_R",
                        "c_MLE_mean_R","Loglik_MLE_R","LogLik_MLE_means_R","Loglik_TrueValues_R")
  
  #Returning final output for Julia (second part)
  LogLike_MLE_Julia$birth_MLE_mean_Julia<-rep(MLE_means_Julia[1],length=nruns)
  LogLike_MLE_Julia$death_MLE_mean_Julia<-rep(MLE_means_Julia[2],length=nruns)
  LogLike_MLE_Julia$catastrophe_MLE_mean_Julia<-rep(MLE_means_Julia[3],length=nruns)                           
  LogLike_MLE_Julia$Loglik_MLE_Julia<-unlist(Loglike_Julia) 
  LogLike_MLE_Julia$LogLik_MLE_means_Julia<- unlist(Loglik_MLE_means_Julia) 
  LogLike_MLE_Julia$Loglik_True_values_Julia <-unlist(LogLik_trueValue_Julia)
  
  
  names(LogLike_MLE_Julia)<-c("b_MLE_Julia","d_MLE_Julia","c_MLE",
                              "b_MLE_mean_Julia","d_MLE_mean_Julia",
                              "c_MLE_mean_Julia","Loglik_MLE_Julia","LogLik_MLE_means_Julia","Loglik_TrueValues_Julia")
  
  
  #returns MLE estimates and its logLik for each realization from R & Julia  
  return(cbind(LogLike_MLE,LogLike_MLE_Julia))  
  
  
  
}

Output=Simul_data_nruns_MLE(X0=2,b=0.5,d=0.3,c=0.001,ti=0,tmax=17,nrep=50,nruns=100)

write.csv(Output,"Loglike_MLE_updated.csv")

