# Loading relevant packages
library(data.table)
library(parallel)
library('latex2exp') #For adding LaTeX symbols to R plots

library(gmp) #Dealing with larger factorial limits
library(compiler)# byte code compilation


# setting working directory
setwd("/Users/clementtwumasi/Desktop/BDC_Paper_Journal_of_Applied_Probability /")

#A function to convert NAs to 0
na.zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
}

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
  #pop <-matrix(-1,1,length(save_ti)) # parasite pop at each location (rows) and time point (cols)
  pop <-matrix(NA,1,length(save_ti)) # parasite pop at observed time point (cols) 
  
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
Exact_BDC_compiler(X0=2,b=.5,d=.3,c=.001,ti=0,tmax=17)

# Tau-leaping algorithm proposed by Gillespie 2001
Hybrid_TauLeap_Gillespie2001=function(X0,b,d,c,error,ti=0,tmax=30){
  
  #ti<-0 #initial time
  #tmax<-30  #final time
  
  rate=numeric(3)
  
  #Nsteps=1 #initializin the number of steps for exact SSA (Nsteps<=100)
  
  
  save_ti <- as.vector(1:tmax)    #Times to simulate
  # leapsize=matrix(NA,1,length(save_ti)) #store leap sizes
  alive <- rep(2, length(save_ti)) # host fish status at each time point
  alive_ti <- 1 #fish starts out alive
  
  save_TF <- rep(FALSE, length(save_ti))
  pop <-matrix(NA,1,length(save_ti)) # parasite pop at observed time point (cols)     
  X<-X0
  
  
  pop_ti <- sum(X)
  
  while(ti<tmax){ ### beginning of while loop
    
    
    ###Computing list of propensity functions/event rates#####
    #probability of birth
    rate[1]<- b*X
    #probability of death
    rate[2]<- d*X
    #probability of catastrophe
    rate[3]<- c*X
    
    
    total_rate<- rate[1]+rate[2]+rate[3] #representing a0(x)
    
    tau= (error*(b+d))/(abs(b-d)*max(b,d))  #estimated leap size for B-D-C process
    
    
    
    #Running Tau-leaping
    if(total_rate*tau >2){ #Execute tau-leaping
      
      ti <- ti +tau 
      if (runif(1) < rate[3]*tau) { # catastrophe
        
        X <- 0
        alive_ti<-2
      } else{ # births and deaths
        X <- X + rpois(1,  rate[1]*tau) - rpois(1, rate[2]*tau)   
      }
      
      
      #print(paste("Pop=",X,"","time=",ti,":Results from Tauleaping"))
      
    } #end of tau-leaping
    
    
    
    else {#Execute the exact SSA
      
      out=SSA_update_event_compiler(X=X,b=b,d=d,c=c)
      X=out$X
      
      ti <- ti +out$t_incr
      
    } #end of exact SSA
    
    
    
    if (X < 0) break #break if there is negative population
    
    if (alive_ti == 2) {
      break
    }
    
    # Events to occur
    save_new <- which((ti >= save_ti) & !save_TF)
    
    for (i in save_new){
      pop[,i]<- pop_ti
      alive[i] <- alive_ti
      #leapsize[,i]=tau
    }
    save_TF <- (ti >= save_ti)
    X<- X
    pop_ti<- sum(X)
    alive_ti <- alive_ti
    
  }
  return(list(pop=pop,alive=alive))
  
}

Hybrid_TauLeap_Gillespie2001_compiler=cmpfun(Hybrid_TauLeap_Gillespie2001)

#Print output from the Simulation Model
Hybrid_TauLeap_Gillespie2001_compiler(X0=2,b=0.5,d=.3,c=.001,error=0.03,ti=0,tmax=17)


# Tau-leaping algorithm proposed by Gillespie & Petzold 2003

#Simulations
Hybrid_TauLeap_Gillespie2003=function(X0,b,d,c,error,ti=0,tmax=30){
  
  
  #ti<-0 #initial time
  #tmax<-30  #final time
  
  
  rate<-numeric(3) #store rate or propensity functions or rates
  
  Nc=0.1#multiplier
  #epsilon=0.03 #recommended error value/bound
  
  Nsteps=1 #initializing the number of steps for exact SSA (Nsteps<=100)
  
  
  Leap_sizes=vector("list",length=4)
  
  save_ti <- as.vector(1:tmax)    #Times to simulate
  
  #leapsize=matrix(NA,1,length(save_ti)) #store leap sizes
  
  alive <- rep(2, length(save_ti)) # host fish status at each time point
  alive_ti <- 1 #fish starts out alive
  
  save_TF <- rep(FALSE, length(save_ti))
  pop <-matrix(NA,1,length(save_ti)) # parasite pop at observed time point (cols)     
  X<-X0
  
  
  pop_ti <- sum(X)
  
  while(ti<tmax){ ### beginning of while loop
    
    
    
    
    
    v=c(+1,-1) # state-change vector for tau-leaping
    ###Computing list of propensity functions/event rates#####
    #probability of birth
    rate[1]<- b*X
    #probability of death
    rate[2]<- d*X
    #probability of catastrophe
    rate[3]<- c*X
    
    
    total_rate<- rate[1]+rate[2]+rate[3] #representing a0(x)
    
    #Computing tau from equation 6 based on Gillespie & Petzold 2003
    Leap_sizes[[1]]= (error*(b+d))/(abs(b^2-b*d))
    Leap_sizes[[2]]= (error*(b+d))/(abs(b*d-d^2))
    Leap_sizes[[3]]= X*(error*(b+d))^2/(b^3+ d*b^2)
    Leap_sizes[[4]]= X*(error*(b+d))^2/(b*d^2+ d^3)
    
    #estimate of the tau leap size
    tau=min(Leap_sizes[[1]],Leap_sizes[[2]],Leap_sizes[[3]],Leap_sizes[[4]])   
    
    #Switching condition
    #condition=Nc*(1/total_rate)
    
    #Running Tau-leaping
    if(total_rate*tau>Nc){#Execute tau-leaping
      
      ti <- ti + tau
      if(runif(1) < rate[3]*tau){
        X <- 0
        alive_ti<-2
      } else{
        
        X <- X + rpois(1,  rate[1]*tau) - rpois(1, rate[2]*tau)  
        
      }
      #print(paste("Pop=",X,"","time=",ti,":Results from Tauleaping"))
      
    } #end of tau-leaping
    
    
    #Running exact SSA algorithm if tau1<Nc*(1/total_rate)
    else {#Execute 1 exact SSA
      
      out=SSA_update_event_compiler(X=X,b=b,d=d,c=c)
      X=out$X
      ti <- ti +out$t_incr
      
      #print(paste("Pop=",X,"","time=",ti,":Results from SSA"))   
      
    }
    
    
    if (X < 0) break #break if there is negative population
    
    if (alive_ti == 2) {
      break
    }
    
    # Events to occur
    save_new <- which((ti >= save_ti) & !save_TF)
    
    for (i in save_new){
      pop[,i]<- pop_ti
      alive[i] <- alive_ti
      #leapsize[,i]=tau #store leap sizes
      
      
    }
    save_TF <- (ti >= save_ti)
    X<- X
    pop_ti<- sum(X)
    alive_ti <- alive_ti
    
  }
  return(list(pop=pop,alive=alive))
  
}

Hybrid_TauLeap_Gillespie2003_compiler=cmpfun(Hybrid_TauLeap_Gillespie2003)

#Print output from the Simulation Model
Hybrid_TauLeap_Gillespie2003_compiler(X0=2,b=.5,d=.3,c=.001,error=0.03,ti=0,tmax=17)

## Function of the Exact mean of the BDC process ##
Expectation_exact=function(b,d,c,t,m){
  
  roots <- sort(Re(polyroot(c(d, -(b+d+c), b))))
  v0 <- roots[1]
  v1 <- roots[2]
  #v0<-((lambda+mu+rho)-sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
  #v1<-((lambda+mu+rho)+sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
  sigma<-exp(-b*(v1-v0)*t)
  k1<-(v0*v1*(1-sigma))/(v1-(sigma*v0))
  k2<-((v1*sigma)-v0)/(v1-(sigma*v0))
  k3<-(1-sigma)/(v1-(sigma*v0))
  expectation=m*(((k1+k2)/(1-k3))^(m-1))*(k2+(k1*k3))*(1-k3)^-2
  return(expectation)
}

#t=seq(1,17) # time
#True_exact_mean=Expectation_exact(b=0.512, d=0.35, c=0.003,t,m=2)
#True_exact_mean

# Function of the Exact variance of the BDC process ##
#True variance from the pgf
Variance_analytic=function(b,d,c,t,m){
    roots <- sort(Re(polyroot(c(d, -(b+d+c), b))))
  v0 <- roots[1]
  v1 <- roots[2]
  #v0<-((lambda+mu+rho)-sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
  #v1<-((lambda+mu+rho)+sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
  sigma<-exp(-b*(v1-v0)*t)
  k1<-(v0*v1*(1-sigma))/(v1-(sigma*v0))
  k2<-((v1*sigma)-v0)/(v1-(sigma*v0))
  k3<-(1-sigma)/(v1-(sigma*v0))
    expectation=m*(((k1+k2)/(1-k3))^(m-1))*(k2+(k1*k3))*(1-k3)^-2
    
Second_derivative_pgf= (2*m*k3*(k2+k1*k3))*((k1+k2)/(1-k3))^(m-1)*(1-k3)^-3 + m*(m-1)*(k2+k1*k3)^2*
    ((k1+k2)/(1-k3))^(m-2)*(1-k3)^-4
     
Variance=(Second_derivative_pgf+ expectation)-(expectation)^2
    return(Variance)
}


#t=seq(1,17) # time
#True_exact_variance=Variance_analytic(b=0.512, d=0.35, c=0.003,t,m=2)
#True_exact_variance

# Function for computing mean pop(Exact SSA)

mean_func=function(x) mean(x,na.rm=T);sd_func=function(x) sd(x,na.rm=T);var_func=function(x) var(x,na.rm=T)
sum_func=function(x) sum(x,na.rm=T)
ExactSSA_mean_sim <- function(X0,b,d,c,ti=0,tmax=30, nrep=1e1, nc=min(16, detectCores())) {
# Birth-death process with catastrophe
# Estimated mean using simulation
#
# Estimate calculated using nrep calls to Hybrid_TauLeap_Gillespie2001_compiler(X0=2,b=b,d=d,c=c,error=error)
# mclapply is used to parallelise the calls, using nc cores (doesn't work on Windows)
#
# Returns sample mean, sample s.e. as a vector length 2

  f <- function(x) Exact_BDC_compiler(X0=X0,b=b,d=d,c=c,ti=ti,tmax=tmax)$pop
  x <- simplify2array(mclapply(1:nrep, f, mc.cores=nc))
  #print(x[,,]) #output combined as a matrix with dimension simulation output length X nrep 
  #return(x[,,])
    pop=x[, ,]
    return(list(mean_pop=apply(pop,1,mean_func),std_error=apply(pop,1,sd_func)/sqrt(nrep)))
}

ExactSSA_err <- function(X0,b,d,c,ti=1,tmax=30, nrep=1e2) {
# Birth-death process with catastrophe: Exact SSA
# A measure of the difference between the theoretical and estimated mean, and the estimation running time
#
# We calculate the squared error loss
# The results will also depend on the choice of parameters: b,d and c
#Returns mean pop, squared error loss, total_Var= sum(std err^2)
  
  time0 <- proc.time()
      sim<- ExactSSA_mean_sim(X0=X0,b=b,d=d,c=c, nrep=nrep,ti=ti,tmax=tmax)
      loss <- sum_func((Expectation_exact(b=b,d=d,c=c,t=1:tmax,m=X0) - sim$mean_pop)^2)
      sum_var <- sum_func((sim$std_error)^2)
    time1 <- proc.time()
    
  return(list(mean_pop=sim$mean_pop,loss=loss,Var=sum_var,Time=time1-time0))
}


#ExactSSA_err(X0=2,b=.5,d=.3,c=.001,ti=0,tmax=17, nrep=10)

#Hybrid tau-leaping (Gillespie 2001)
mean_func=function(x) mean(x,na.rm=T); sd_func=function(x) sd(x,na.rm=T);sum_func=function(x) sum(x,na.rm=T)

HTL2001_mean_sim <- function(X0,b,d,c,ti=0,tmax=30, error, nrep=1e1, nc=min(16, detectCores())) {
# Birth-death process with catastrophe
# Estimated mean using simulation
#
# Estimate calculated using nrep calls to Hybrid_TauLeap_Gillespie2001_compiler(X0=2,b=b,d=d,c=c,error=error)
# mclapply is used to parallelise the calls, using nc cores (doesn't work on Windows)
#
# Returns sample mean, sample s.e. and average leap size as a vector length 3
  f <- function(x) Hybrid_TauLeap_Gillespie2001_compiler(X0=2,b=b,d=d,c=c,error=error,ti=ti,tmax=tmax)$pop
  x <- simplify2array(mclapply(1:nrep, f, mc.cores=nc))
  #pop= do.call(rbind,x["pop",]) #storing population
  
   pop=x[, ,]
    return(list(mean_pop=apply(pop,1,mean_func),std_error=apply(pop,1,sd_func)/sqrt(nrep)))

}

#HTL2001_mean_sim(X0=2,b=.5,d=.3,c=.001,error=0.03,ti=0,tmax=17,nrep=5)

HTL2001sim_err <- function(X0,b,d,c,ti=0, tmax=30, error, nrep=1e1) {
# Birth-death process with catastrophe: Tau-leaping algorithm (Gillespie 2001)
# A measure of the difference between the theoretical and estimated mean, and the estimation time
#
# If tau==0 then the difference is purely due to sampling error. 
# If tau>0 then there will also be bias
# We calculate the squared error 
# The choice of ti and error bound will affect the variability of the estimates and the effect of tau-leaping
# (tau-leaping only has an effect when the population is large)
# The results will also depend on the choice of la, mu and ga
#
# We calculate the squared error loss
# The results will also depend on the choice of parameters: b,d and c
#Returns mean pop, squared error loss, total_Var= sum(std err^2)
  
 
  
  time0 <- proc.time()
    
      sim <- HTL2001_mean_sim(X0=X0,b=b,d=d,c=c,error=error,nrep=nrep,ti=ti,tmax=tmax)
    
      loss <- sum_func((Expectation_exact(b=b,d=d,c=c,t=1:tmax,m=X0) - sim$mean_pop)^2)
      sum_var <- sum_func((sim$std_error)^2)
     
  
  time1 <- proc.time()
    return(list(mean_pop=sim$mean_pop,loss=loss,Var=sum_var,Time=time1-time0))
}

#HTL2001sim_err(X0=2,b=.5,d=.3,c=.001,ti=0,tmax=17,error=0.03,nrep=5)
#HTL2003_mean_sim(X0=2,b=.5,d=.3,c=.001,error=0.03,ti=0,tmax=17,nrep=5)

#Hybrid tau-leaping (Gillespie 2001)
mean_func=function(x) mean(x,na.rm=T); sd_func=function(x) sd(x,na.rm=T);sum_func=function(x) sum(x,na.rm=T)

HTL2003_mean_sim <- function(X0,b,d,c,ti=0,tmax=30, error, nrep=1e1, nc=min(16, detectCores())) {
# Birth-death process with catastrophe
# Estimated mean using simulation
#
# Estimate calculated using nrep calls to Hybrid_TauLeap_Gillespie2001_compiler(X0=2,b=b,d=d,c=c,error=error)
# mclapply is used to parallelise the calls, using nc cores (doesn't work on Windows)
#
# Returns sample mean, sample s.e. and average leap size as a vector length 3
  f <- function(x) Hybrid_TauLeap_Gillespie2003_compiler(X0=2,b=b,d=d,c=c,error=error,ti=ti,tmax=tmax)$pop
  x <- simplify2array(mclapply(1:nrep, f, mc.cores=nc))
  #pop= do.call(rbind,x["pop",]) #storing population
  
   pop=x[, ,]
    return(list(mean_pop=apply(pop,1,mean_func),std_error=apply(pop,1,sd_func)/sqrt(nrep)))

}



HTL2003sim_err <- function(X0,b,d,c,ti=0, tmax=30, error, nrep=1e1) {
# Birth-death process with catastrophe: Tau-leaping algorithm (Gillespie 2001)
# A measure of the difference between the theoretical and estimated mean, and the estimation time
#
# If tau==0 then the difference is purely due to sampling error. 
# If tau>0 then there will also be bias
# We calculate the squared error 
# The choice of ti and error bound will affect the variability of the estimates and the effect of tau-leaping
# (tau-leaping only has an effect when the population is large)
# The results will also depend on the choice of la, mu and ga
#
# We calculate the squared error loss
# The results will also depend on the choice of parameters: b,d and c
#Returns mean pop, squared error loss, total_Var= sum(std err^2)
  
 
  
  time0 <- proc.time()
    
      sim <- HTL2003_mean_sim(X0=X0,b=b,d=d,c=c,error=error,nrep=nrep,ti=ti,tmax=tmax)
    
      loss <- sum_func((Expectation_exact(b=b,d=d,c=c,t=1:tmax,m=X0) - sim$mean_pop)^2)
      sum_var <- sum_func((sim$std_error)^2)
     
  
  time1 <- proc.time()
    return(list(mean_pop=sim$mean_pop,loss=loss,Var=sum_var,Time=time1-time0))
}

#HTL2003sim_err(X0=2,b=.5,d=.3,c=.001,ti=0,tmax=17,error=0.03,nrep=5)

error_thresholds<-c(seq(0, 0.01, 0.002), seq(0.02, 0.1, 0.01))
param_case1=c(0.5,0.3,0.001);param_case2=c(2,1,0.01);param_case3=c(3,2,0.1)

#Case 1: Hybrid Tau-leaping algorithm Gillespie 2001
out_HTL2001_case1<-mclapply(seq_along(error_thresholds),function(i)
HTL2001sim_err(X0=2,b=param_case1[1],d=param_case1[2],c=param_case1[3],ti=0,tmax=30,
                                error=error_thresholds[i],nrep=2e6),mc.cores=min(16, detectCores()))

#Case 2: Hybrid Tau-leaping algorithm Gillespie 2001
out_HTL2001_case2<-mclapply(seq_along(error_thresholds),function(i)
HTL2001sim_err(X0=2,b=param_case2[1],d=param_case2[2],c=param_case2[3],ti=0,tmax=30,
                                error=error_thresholds[i],nrep=2e6),mc.cores=min(16, detectCores()))

#Case 3: Hybrid Tau-leaping algorithm Gillespie 2001
out_HTL2001_case3<-mclapply(seq_along(error_thresholds),function(i)
HTL2001sim_err(X0=2,b=param_case3[1],d=param_case3[2],c=param_case3[3],ti=0,tmax=30,
                                error=error_thresholds[i],nrep=2e6),mc.cores=min(16, detectCores()))

#Case 1: Hybrid Tau-leaping algorithm Gillespie & Petzold 2003
out_HTL2003_case1<-mclapply(seq_along(error_thresholds),function(i)
HTL2003sim_err(X0=2,b=param_case1[1],d=param_case1[2],c=param_case1[3],ti=0,tmax=30,
error=error_thresholds[i],nrep=2e6),mc.cores=min(16, detectCores()))

#Case 2: Hybrid Tau-leaping algorithm Gillespie & Petzold 2003
out_HTL2003_case2<-mclapply(seq_along(error_thresholds),function(i)
HTL2003sim_err(X0=2,b=param_case2[1],d=param_case2[2],c=param_case2[3],ti=0,tmax=30,
                                error=error_thresholds[i],nrep=2e6),mc.cores=min(16, detectCores()))

#Case 3: Hybrid Tau-leaping algorithm Gillespie & Petzold 2003
out_HTL2003_case3<-mclapply(seq_along(error_thresholds),function(i)
HTL2003sim_err(X0=2,b=param_case3[1],d=param_case3[2],c=param_case3[3],ti=0,tmax=30,
                                error=error_thresholds[i],nrep=2e6),mc.cores=min(16, detectCores()))

Sim_means2001_case1=data.frame(matrix(NA,ncol=length(error_thresholds),nrow=30))
Sim_means2001_case2=data.frame(matrix(NA,ncol=length(error_thresholds),nrow=30))
Sim_means2001_case3=data.frame(matrix(NA,ncol=length(error_thresholds),nrow=30))
Sim_means2003_case1=data.frame(matrix(NA,ncol=length(error_thresholds),nrow=30))
Sim_means2003_case2=data.frame(matrix(NA,ncol=length(error_thresholds),nrow=30))
Sim_means2003_case3=data.frame(matrix(NA,ncol=length(error_thresholds),nrow=30))


for(i in 1:length(error_thresholds)){
  Sim_means2001_case1[,i]<- out_HTL2001_case1[[i]]$mean_pop
  Sim_means2001_case2[,i]<- out_HTL2001_case2[[i]]$mean_pop
  Sim_means2001_case3[,i]<- out_HTL2001_case3[[i]]$mean_pop
  Sim_means2003_case1[,i]<- out_HTL2003_case1[[i]]$mean_pop
  Sim_means2003_case2[,i]<- out_HTL2003_case2[[i]]$mean_pop
  Sim_means2003_case3[,i]<- out_HTL2003_case3[[i]]$mean_pop
}

names(Sim_means2001_case1)=paste("error=",error_thresholds)
names(Sim_means2001_case2)=names(Sim_means2001_case1);names(Sim_means2001_case3)=names(Sim_means2001_case1)

names(Sim_means2003_case1)=paste("error=",error_thresholds)
names(Sim_means2003_case2)=names(Sim_means2003_case1);names(Sim_means2003_case3)=names(Sim_means2003_case1)

#Saving mean population numbers for the two Tau-leaping algorithms for the 3 different cases
write.csv(Sim_means2001_case1,"Sim_means2001_case1.csv")
write.csv(Sim_means2001_case2,"Sim_means2001_case2.csv")
write.csv(Sim_means2001_case3,"Sim_means2001_case3.csv")
write.csv(Sim_means2003_case1,"Sim_means2003_case1.csv")
write.csv(Sim_means2003_case2,"Sim_means2003_case2.csv")
write.csv(Sim_means2003_case3,"Sim_means2003_case3.csv")

# Importing the results
Sim_means2001_results=Sim_means2003_results=vector("list",length=3)
Sim_means2001=Sim_means2003=vector("list",length=3)

Sim_means2001_results[[1]]<- read.csv(file="Sim_means2001_case1.csv")
Sim_means2001_results[[2]]<- read.csv(file="Sim_means2001_case2.csv")
Sim_means2001_results[[3]]<- read.csv(file="Sim_means2001_case3.csv")
Sim_means2003_results[[1]]<- read.csv(file="Sim_means2003_case1.csv")
Sim_means2003_results[[2]]<- read.csv(file="Sim_means2003_case2.csv")
Sim_means2003_results[[3]]<- read.csv(file="Sim_means2003_case3.csv")

#Tau-leaping (Gillespie 2001)
for(i in 1:3) Sim_means2001[[i]]<- Sim_means2001_results[[i]][,-1]

#Tau-leaping (Gillespie & Petzold 2003)
for(i in 1:3) Sim_means2003[[i]]<- Sim_means2003_results[[i]][,-1]

Parasite_data_summary<- function(observed_times= 1:9,Total_fish=50,
                                parameter_values, seed_number=123){

set.seed(seed_number) #for reproducibility
# sum_func<- function(x) sum(x, na.rm = T)
    
parasites_data_SSA<- data.frame(matrix(NA,nrow=Total_fish, ncol=length(observed_times)+1)) 
names(parasites_data_SSA)<- c(paste0("Parasites_t",observed_times), "SimulationMethod")
    
parasites_data_HTL2001<- data.frame(matrix(NA,nrow=Total_fish, ncol=length(observed_times)+1))
names(parasites_data_HTL2001)<- names(parasites_data_SSA)
    
parasites_data_HTL2003<- data.frame(matrix(NA,nrow=Total_fish, ncol=length(observed_times)+1))
names(parasites_data_HTL2003)<- names(parasites_data_SSA)



    for(i in 1:Total_fish){
        # Exact SSA
       parasites_SSA<- Exact_BDC_compiler(X0=2,b=parameter_values[1],
                   d=parameter_values[2],
                   c=parameter_values[3],ti=0,tmax=17)$pop
        # extracting data at time 1, 3, 5, 7, 9, 11, 13, 15 and 17
       parasites_data_SSA[i, observed_times]<- na.zero(parasites_SSA[seq(1,17, by=2)])
       parasites_data_SSA[i, length(observed_times)+1]<- "Exact SSA"

        # Hybrid TauLeap Gillespie 2001
       parasites_HTL2001<- Hybrid_TauLeap_Gillespie2001_compiler(X0=2,b=parameter_values[1],
                   d=parameter_values[2],
                   c=parameter_values[3],ti=0,tmax=17,error=0.01)$pop
        # extracting data at time 1, 3, 5, 7, 9, 11, 13, 15 and 17
       parasites_data_HTL2001[i, observed_times]<- na.zero(parasites_HTL2001[seq(1,17, by=2)])
       parasites_data_HTL2001[i, length(observed_times)+1]<- "HTL2001"

      
        # Hybrid TauLeap Gillespie 2003
       parasites_HTL2003<- Hybrid_TauLeap_Gillespie2003_compiler(X0=2,b=parameter_values[1],
                   d=parameter_values[2],
                   c=parameter_values[3],ti=0,tmax=17,error=0.01)$pop
        # extracting data at time 1, 3, 5, 7, 9, 11, 13, 15 and 17
       parasites_data_HTL2003[i, observed_times]<- na.zero(parasites_HTL2003[seq(1,17, by=2)])
       parasites_data_HTL2003[i, length(observed_times)+1]<- "HTL2003"



        }

    # Exporting the final results
    combined_data<- rbind(parasites_data_SSA, parasites_data_HTL2001,parasites_data_HTL2003)
     return(combined_data)
    }

true.parameters.case1<- c(0.5,0.3, 0.001)
true.parameters.case2<- c(2,1, 0.01)
true.parameters.case3<- c(3,2, 0.1)

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

transform_structure<- function(data){
    data.timepoints<- list(); times<- seq(1,17,by=2)

    for(t in 1:9){
        data.timepoints[[t]]<- data[, c(t,10)]
        data.timepoints[[t]][,1]<- log(data.timepoints[[t]][,1]+1) #on log scale
        names(data.timepoints[[t]])[1]<- "Parasites"
        data.timepoints[[t]]$Time<- times[t]
    }
    
    
    return(do.call("rbind",data.timepoints))
}

Simulation.Results<- NULL; transdata=NULL

# Case 1
Simulation.Results[[1]]<- Parasite_data_summary(observed_times= 1:9,Total_fish=50,
                    parameter_values=true.parameters.case1, 
                      seed_number=123)



# Case 2
Simulation.Results[[2]]<- Parasite_data_summary(observed_times= 1:9,Total_fish=50,
                    parameter_values=true.parameters.case2, 
                      seed_number=123)



# Case 3
Simulation.Results[[3]]<- Parasite_data_summary(observed_times= 1:9,Total_fish=50,
                    parameter_values=true.parameters.case3, 
                      seed_number=123)



dim(Simulation.Results[[1]])

head(Simulation.Results[[1]])

case<- 1
transdata[[case]]<- transform_structure(data=Simulation.Results[[case]])
transdata[[case]]$Day<- as.factor(transdata[[case]]$Time)
levels(transdata[[case]]$Day)<- c("Day 1","Day 3","Day 5", "Day 7", "Day 9",
                              "Day 11","Day 13","Day 15","Day 17")

transdata[[case]]$Case<- rep("Case 1", dim(transdata[[case]])[1])

head(transdata[[1]])

case<- 2
transdata[[case]]<- transform_structure(data=Simulation.Results[[case]])
transdata[[case]]$Day<- as.factor(transdata[[case]]$Time)
levels(transdata[[case]]$Day)<- c("Day 1","Day 3","Day 5", "Day 7", "Day 9",
                              "Day 11","Day 13","Day 15","Day 17")

transdata[[case]]$Case<- rep("Case 2", dim(transdata[[case]])[1])

case<- 3
transdata[[case]]<- transform_structure(data=Simulation.Results[[case]])
transdata[[case]]$Day<- as.factor(transdata[[case]]$Time)
levels(transdata[[case]]$Day)<- c("Day 1","Day 3","Day 5", "Day 7", "Day 9",
                              "Day 11","Day 13","Day 15","Day 17")

transdata[[case]]$Case<- rep("Case 3", dim(transdata[[case]])[1])



theme_set(
  theme_bw() +
    theme(legend.position = "top")
  )



p1<- do.call("rbind",transdata)%>%
  mutate(type=fct_reorder(as.factor(Time),Parasites),
        Group=fct_reorder(as.factor(SimulationMethod),Parasites)) %>%
     ggplot(fill = factor(SimulationMethod, levels=positions))+
geom_smooth(aes(x=Time, y=Parasites,color=SimulationMethod,group=SimulationMethod),
                          method = "loess",se = FALSE)+
labs(y="Parasite mean intensity") +
theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.background = element_blank())+
scale_color_manual(values=c("blue","green","orange"))+
xlim(1,17)+facet_wrap(~Case,  ncol=1)+ylim(1,17)+
  guides(color = guide_legend(title = "Comparable groups:"))


p1

# Visualise: Specify the comparisons you want
case=1
my_cols<- c("blue","green","orange")

my_comparisons <- list( c("Exact SSA", "HTL2001"),
                       c("Exact SSA", "HTL2003"),
             c("HTL2001", "HTL2003"))

#subset(transdata_reg, Time==1)
bxp1<- ggboxplot(transdata[[case]], x = "SimulationMethod", y = "Parasites",
          fill = "SimulationMethod")+scale_fill_manual(values=my_cols) #you can change "fill" to "color"

bxp1<- bxp1+stat_compare_means(comparisons = my_comparisons, hide.ns = F,
    label = c("p.signif"),palette = c('red'),label.x.npc ='center',#ref.group =".all."
  label.y.npc ='top')+ # Add pairwise comparisons p-value
  stat_compare_means(method = "kruskal.test", label.y = 10,label.x=1.5,size = 3)+     # Add global p-value
 theme(text = element_text(size = 10))+  
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
facet_wrap( ~ Day)

bxp1<- bxp1+labs(x="Simulation Methods",y="Parasite abundance (log-scale)",fill = "Simulation Methods:")+
geom_dotplot(binaxis = "y", stackdir = "center", binwidth = .05)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


bxp1



library("FactoMineR")
library("factoextra")

head(Simulation.Results[[1]])

methods<- Simulation.Results[[1]]$SimulationMethod

combined_data_parasite_numbers<- Simulation.Results[[1]][, 1:9]
combined_data_parasite_numbers<- log(combined_data_parasite_numbers+.001)

# PCA
parasite_numbers.pca <- PCA(combined_data_parasite_numbers, graph = FALSE)
print(parasite_numbers.pca)

#The sum of all the eigenvalues give a total variance
get_eigenvalue(parasite_numbers.pca )#: Extract the eigenvalues/variances of principal components

var <- get_pca_var(parasite_numbers.pca)
var

# Coordinates
head(var$coord)
# Cos2: quality on the factor map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

fviz_cos2(parasite_numbers.pca, choice = "var", axes = 1:2, sort.val="none")+
            labs(title="Cos2 of variables to Dim 1-2")



##Scree plot
#Proportion of variance explained by each principal component
#p1=fviz_eig(parasite_numbers.pca, addlabels = TRUE, ylim = c(0, 65),xlim=c(1,10),
   #         hjust = 0.3,main="A: Scree plot")+
#theme(axis.line = element_line(),
#        panel.grid.major = element_blank(),
 #       panel.grid.minor = element_blank(),

   #     panel.background = element_blank())


p1<- fviz_cos2(parasite_numbers.pca, choice = "var", axes = 1:2)+
            labs(title="A: Cos2 of variables to Dim 1-2")+
theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.background = element_blank())

#Biplot
# Color by cos2 values: quality on the factor map
#Cos2 of variables to Dim 1-2
p2<- fviz_pca_biplot(parasite_numbers.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
                   label = "var" , addEllipses=T , ellipse.level=0.95,title="B: PCA- Biplot")+
theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.background = element_blank())



gridExtra::grid.arrange(p1,p2, ncol=1,nrow=2)

#Parasite distribution/numbers  - (Individual- PCA)
p1<- fviz_pca_ind(parasite_numbers.pca,
             # label = "none", # hide individual labels
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = methods, # color by groups
             palette = my_cols,
             addEllipses = TRUE, # Concentration ellipses
            #ellipse.type = "confidence",
             label = "var",
             legend.title = "Simulation Methods:",
             axes=c(1,2),
            repel = TRUE,
            title=""
             )+
theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        panel.background = element_blank())


p1+theme(legend.position="top")

library(vegan)

#The function computes dissimilarity indices that are useful for or popular with community ecologists
methods<- Simulation.Results[[1]]$SimulationMethod

combined_data_parasite_numbers<- Simulation.Results[[1]][, 1:9]
combined_data_parasite_numbers<- log(combined_data_parasite_numbers+.00001)


diss_indices <- vegdist(combined_data_parasite_numbers,na.rm = TRUE,method="euclidean")



methods<- as.factor(methods)

levels(methods)


## Calculate multivariate dispersions
multivariate_dispersion_mod <- betadisper(d=diss_indices, group=methods,  
                                          bias.adjust = T)
multivariate_dispersion_mod 

## Perform test
anova(multivariate_dispersion_mod)[1, ]


permutest(multivariate_dispersion_mod, permutations = 99)

TukeyHSD(multivariate_dispersion_mod)

plot(multivariate_dispersion_mod, ellipse = F, hull = T,label = TRUE, label.cex = .54,main="",col=my_cols, 
xlab="PCoA 1",ylab="PCoA 2",lwd=2, pch=c(6,8, 19),seg.lwd=.5,ylim=c(-15,15)) # 1 sd data ellipse

text(0, 20, "Multivariate homogeneity test of variances in the simulated data: p-value=0.972", col="black", font=2)
# mtext(LETTERS[1], adj=0, line=1,font=2)

distance_df<- data.frame(distance_centroid=multivariate_dispersion_mod$distances, 
                         SimulationMethod=multivariate_dispersion_mod$ group )

dim(distance_df)
head(distance_df)

my_cols<- c("blue","green","orange")

my_comparisons <- list( c("Exact SSA", "HTL2001"),
                       c("Exact SSA", "HTL2003"),
             c("HTL2001", "HTL2003"))

#subset(transdata_reg, Time==1)
bxp1<- ggboxplot(distance_df, x = "SimulationMethod", y = "distance_centroid",
          fill = "SimulationMethod")+scale_fill_manual(values=my_cols) #you can change "fill" to "color"

bxp1<- bxp1+stat_compare_means(comparisons = my_comparisons, hide.ns = F,
    label = c("p.signif"),palette = c('red'),label.x.npc ='center',#ref.group =".all."
  label.y.npc ='top')+ # Add pairwise comparisons p-value
  stat_compare_means(method = "kruskal.test", label.y = 50,label.x=.8,size = 4)+     # Add global p-value
 theme(text = element_text(size = 11))+  
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))


bxp1<- bxp1+labs(x="Simulation Methods",y="Distance to Centroid",fill = "Simulation Methods:")+
geom_dotplot(binaxis = "y", stackdir = "center", binwidth = .5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


bxp1


