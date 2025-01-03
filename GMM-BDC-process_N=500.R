# Method-of-moment is considered as an ad-hoc estimation to the MLE


library(parallel)
library(gmm)
library(compiler)


setwd("/home/clement/Clement_R_Results_folder/GMM_MLE_results_time")
#options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

ncores=16
ncores
## Function of the Exact mean of the BDC process ##
First_moment=function(b,d,c,t,m){
  
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

#t=seq(1,30) # time
#First_moment_values=First_moment(b=0.512, d=0.35, c=0.003,t,m=2)
#First_moment_values




#Second moment
Second_moment=function(b,d,c,t,m){
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
  
  Second_moment_results<- Variance + expectation^2
  return(Second_moment_results)
}


#t=seq(1,30) # time
#Second_moment_values=Second_moment(b=0.512, d=0.35, c=0.003,t,m=2)
#Second_moment_values

#Third moment

Third_moment=function(b,d,c,t,m){
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
  
  Third_derivative_pgf<- 6*m*(k2+k1*k3)*(k3^2)*(((k1+k2)/(1-k3))^(m-1))*(1-k3)^(-4)+
    6*m*(m-1)*((k2+k1*k3)^2)*k3*(((k1+k2)/(1-k3))^(m-2))*(1-k3)^(-5)+
    m*(m-1)*(m-2)*((k2+k1*k3)^3)*(((k1+k2)/(1-k3))^(m-3))*(1-k3)^(-6)   
  
  Variance=(Second_derivative_pgf+ expectation)-(expectation)^2
  
  Second_moment_results<- Variance + expectation^2
  
  Third_moment_results<- Third_derivative_pgf+(3*Second_moment_results)-(2*expectation)
  return(Third_moment_results)
}


#t=seq(1,30) # time
#Third_moment_values=Second_moment(b=0.512, d=0.35, c=0.003,t,m=2)
#Third_moment_values



#Function for constants and PGF
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
Prob_catastrophe=function(lambda,mu,rho,t,z=1,m=2){
  constant=1-PGF_z_compiler(lambda=lambda,mu=mu,rho=rho,t=t,z=z,m=m)
  return(constant)
}




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
#parameters<- c(0.5,0.3,0.001)
#Exact_BDC_compiler(X0=2,b=parameters[1],d=parameters[2],c=parameters[3],ti=0,tmax=17)


#Function to simulate data for group of fish (nrep)
#For GMM
#Simulation that return both the parasite numbers and survival status
ExactSSA_sim <- function(X0,b,d,c,ti=0,tmax=30, nrep=1e1) {


  f <- function(x) Exact_BDC_compiler(X0=X0,b=b,d=d,c=c,ti=ti,tmax=tmax)
  x <- lapply(1:nrep, f)
  Out<-t(do.call("rbind",x)) #rbinding simulated data on nrep fish (each column for a fish)
  pop<-t(do.call("rbind",Out["pop", ]))
  alive<-t(do.call("rbind",Out["alive", ])) 
    return(list(pop=pop,alive=alive))
}

data<-ExactSSA_sim(X0=2,b=.5,d=.3,c=.001,ti=0,tmax=17, nrep=5)
#data


#Set the catastrophe state -1 to 0
zero.catastrophe <- function (x) {
  x[x<0] <- 0
  return(x)
}

#function for sample moments
sample_mean_1st<- function(x) sum(x)/length(x)
sample_mean_2nd<- function(x) sum(x^2)/length(x)
sample_mean_3rd<- function(x) sum(x^3)/length(x)


# Function to return estimates from the Method of moments

Simul_data_nruns_GMM_own<- function(X0,b,d,c,ti=0,tmax=30, nrep=1e1,nruns=2, nc=ncores) {
  
  Sim_data<-array(dim = c(nruns,9,nrep))
  survival_data<-array(dim = c(nruns,9,nrep))
  Prob_catastrophe_sample<- matrix(0,nrow=9,ncol=nruns)
  
  x1_out <- mclapply(1:nruns,
                     function(x) ExactSSA_sim(X0=X0,b=b,d=d,c=c,ti=ti,tmax=tmax,nrep=nrep), mc.cores=nc)
  
  
  
  for (k in 1:nruns){
    for(i in 1:nrep){
      # saving results for each 100 realizations and 
      Sim_data[k, ,i]<- x1_out[[k]]$pop[, i][seq(1,17,by=2)]
      survival_data[k, ,i]<- x1_out[[k]]$alive[, i][seq(1,17,by=2)]
      #print( survival_data[k, ,i])
      
    }
    
  }
  
  #Computing catastrophic probability analytically and based on the sample data                   
  time<-seq(1,17,by=2)
  
  #computing sample prob of catastrophe                   
  
  time_index<- seq_along(time)
  for (k in 1:nruns){ 
    for(i in time_index){
      if(any(survival_data[k,i ,]==2)==TRUE){
        #print(paste("time=",time[i]))
        fish_dead_sim<-length(which(survival_data[k,i ,]==2))
        #print( fish_dead_sim)
        Prob_catastrophe_sample[i,k]<-fish_dead_sim/dim(survival_data[k, ,])[2]
      }
    }
    
  }
  
  
  g_objectivefunc_firstStep <- function(x,prob_sample,fixed=c(FALSE,FALSE,FALSE)) {
    Prob_catastrophe_analytical<- rep(NA,length=length(time))
    params<-fixed
    function(p){
      params[!fixed]<-p
      #The three parameters to be optimized
      b1<-params[1]
      d1<-params[2]
      c1<-params[3]
      
      
      #Computing theoritical prob of catastrophe
      for(i in seq_along(time)){
        Prob_catastrophe_analytical[i]<- Prob_catastrophe(lambda=b1,mu=d1,rho=c1,t=time[i])                 
      } 
      
      m1 <- First_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_1st)
      m2 <- Second_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_2nd)
      m3 <- Third_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_3rd)
      
      Catastrophe_Prob<- Prob_catastrophe_analytical- prob_sample
      
      gbar_theta<-c(mean(m1),mean(m2),mean(m3),mean(Catastrophe_Prob))
      
      Objective_func<- t(gbar_theta)%*%gbar_theta
      
    } 
    
  }
  
  
  #First step of GMM
  GMM_firstStep<-function(prob_sample,x){
    objec_func<- g_objectivefunc_firstStep(x=x,prob_sample=prob_sample)
    initial<-c(2, 1, 0.001)# initial values to optimize over
    estimates=constrOptim(initial, objec_func, NULL, 
                          ui=rbind(c(1,0,0),  # lambda>0 
                                   c(0,1,0),    # mu >0
                                   c(0,0,1) # rho > 0
                          ),       
                          ci=c(0,0,0),method='Nelder-Mead')$par # the thresholds of the constraints (RHS)$par
    
    return(estimates)
    
  }
  
  GMM_resultsStep1<- mclapply(1:nruns,function(k) 
    GMM_firstStep(prob_sample=Prob_catastrophe_sample[,k],x=as.data.frame(Sim_data[k, ,])),mc.cores=nc) 
  
  
  #Second-step of the optimization
  
  Weight<-function(x,prob_sample,estimate1){ 
    est_step1<- c(estimate1)
    Prob_catastrophe_analytical1<- rep(NA,length=length(time))
    #Computing theoritical prob of catastrophe
    for(i in seq_along(time)){
      Prob_catastrophe_analytical1[i]<- Prob_catastrophe(lambda=est_step1[1],mu=est_step1[2],
                                                         rho=est_step1[3],t=time[i])                 
    } 
    
    m1 <- First_moment(b=est_step1[1],d=est_step1[2],c=est_step1[3],t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_1st)
    m2 <- Second_moment(b=est_step1[1],d=est_step1[2],c=est_step1[3],t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_2nd)
    m3 <- Third_moment(b=est_step1[1],d=est_step1[2],c=est_step1[3],t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_3rd)
    Catastrophe_Prob<- Prob_catastrophe_analytical1- prob_sample
    
    g<-cbind(m1,m2,m3,Catastrophe_Prob)
    
    covariance_matrix<- cov(g)
    #Setting off-diagonals to 0 to obtain an invertible weighting (diagonal) matrix 
    #by assuming that the moment conditions are uncorrelated
    covariance_matrix[lower.tri(covariance_matrix)] <- 0
    covariance_matrix[upper.tri(covariance_matrix)] <- 0
    weightmatrix<- solve(covariance_matrix)
  }
  
  
  weighting_matrix_cov<- lapply(1:nruns,function(k) Weight(x=as.data.frame(Sim_data[k, ,]),
                                                           prob_sample=Prob_catastrophe_sample[,k],estimate1=GMM_resultsStep1[[k]]))
  
  #Second optimization step                               
  g_objectivefunc_2ndStep <- function(x,prob_sample,weighting_matrix,fixed=c(FALSE,FALSE,FALSE)) {
    Prob_catastrophe_analytical<- rep(NA,length=length(time))
    params<-fixed
    function(p){
      params[!fixed]<-p
      #The three parameters to be optimized
      b1<-params[1]
      d1<-params[2]
      c1<-params[3]
      
      
      #Computing theoritical prob of catastrophe
      for(i in seq_along(time)){
        Prob_catastrophe_analytical[i]<- Prob_catastrophe(lambda=b1,mu=d1,rho=c1,t=time[i])                 
      } 
      
      m1 <- First_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_1st)
      m2 <- Second_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_2nd)
      m3 <- Third_moment(b=b1, d=d1,c=c1,t=seq(1,17,by=2),m=2)-apply(zero.catastrophe(x),1,sample_mean_3rd)
      
      Catastrophe_Prob<- Prob_catastrophe_analytical- prob_sample
      
      gbar_theta<-c(mean(m1),mean(m2),mean(m3),mean(Catastrophe_Prob))
      
      Objective_func<- t(gbar_theta)%*%weighting_matrix%*%gbar_theta
      
    } 
    
  }
  
  #second step of GMM
  GMM_2ndStep<-function(prob_sample,x,weighting_matrix){
    objec_func<- g_objectivefunc_2ndStep(x=x,prob_sample=prob_sample,weighting_matrix=weighting_matrix)
    initial<-c(2, 1, 0.001)# initial values to optimize over
    estimates=constrOptim(initial, objec_func, NULL, 
                          ui=rbind(c(1,0,0),  # lambda>0 
                                   c(0,1,0),    # mu >0
                                   c(0,0,1) # rho > 0
                          ),       
                          ci=c(0,0,0),method='Nelder-Mead')$par # the thresholds of the constraints (RHS)$par
    
    return(estimates)
    
  }
  
  time0=proc.time()
  GMM_resultsStep2<- mclapply(1:nruns,function(k) 
    GMM_2ndStep(prob_sample=Prob_catastrophe_sample[,k],x=as.data.frame(Sim_data[k, ,]),
                weighting_matrix=weighting_matrix_cov[[k]]),mc.cores=nc) 
  time1=proc.time()-time0
  
  GMM_estimates<-do.call("rbind",GMM_resultsStep2)
  time<-sum(as.vector(time1)[-3])
  return(list(GMM_estimates=GMM_estimates,Computational_time=time))                              
}



#A function to compute the mean of MM estimates, bias, variance,MSE & Confidence intervals

Summary_peformance<-function(MME_data,true_parameter_values){
    #Computing the mean estimates
    means_estimates<-apply(MME_data,2,mean)
    MM_estimates<- as.vector(means_estimates)
    
    #computing the variance
    var_estimates<-apply(MME_data,2,var)
    var_MME<- as.vector(var_estimates)
    
    #Computing bias for each parameter
    bias_birth<-  MM_estimates[1]-true_parameter_values[1]
    bias_death<-  MM_estimates[2]-true_parameter_values[2] 
    bias_catastrophe<-  MM_estimates[3]-true_parameter_values[3]
    bias_MME<-c(bias_birth, bias_death,bias_catastrophe)
    #Total_bias<- bias_birth+ bias_death+bias_catastrophe
 
      
    #Computing MSE
    
    MSE_birth<- var_MME[1]+bias_birth^2
    MSE_death<- var_MME[2]+bias_death^2
    MSE_catastrophe<- var_MME[3]+bias_catastrophe^2
    MSE_MME<-c(MSE_birth,MSE_death, MSE_catastrophe)
    
    #Computing confidence interval
    n<- dim(MME_data)[1]
  number_of_parameters<- dim(MME_data)[2]
  MM_mean_estimates=stderr=Lower_lim95=Upper_lim95=rep(NA,number_of_parameters) 
 for (i in 1:number_of_parameters){
 #Expected ML estimates       
 MM_mean_estimates[i]<-   MM_estimates[i]
    
  #std error= sqrt(Var/n)  
 stderr[i]<- sqrt(var_MME[i]/n)
  
  #Lower limits  
  Lower_lim95[i]<- MM_mean_estimates[i]- (1.96* stderr[i])
    
  #Upper limits  
  Upper_lim95[i]<- MM_mean_estimates[i]+ (1.96* stderr[i])
    }
       
    
output_accuracy<- data.frame(GMM_estimates=MM_estimates,true_parameter_values=true_parameter_values,
    bias_GMM=bias_MME,var_GMM=var_MME, MSE_GMM=MSE_MME)
    rownames(output_accuracy)<- c("birth rate","death_rate","catastrophe_rate")
    
     output_accuracy$Lower_lim95<-Lower_lim95
     output_accuracy$Upper_lim95<-Upper_lim95

    
    return(output_accuracy)
}


#Simulation for one realization (50 fish) for 100 replicates**


#Checking for consistency of Generalized  Method of method estimator
#Increasing the number of simulation for one realization to 100 for 100 replicates**

Time_results2<- data.frame(matrix(NA,nrow=3,ncol=1))
names(Time_results2)<-c("GMM_time") 
rownames(Time_results2)<-c("case1","case2","case3")

#case 1
options(warn=-1)

GMM_results_case1_100<-Simul_data_nruns_GMM_own(X0=2,b=0.5,d=0.3,c=.001,ti=0,tmax=17, nrep=100,nruns=100)

write.csv(GMM_results_case1_100$GMM_estimates,"GMM_results_case1_100.csv")

#Saving accuracy measures
GMM_case1_accuracy_100<-Summary_peformance(MME_data=GMM_results_case1_100$GMM_estimates,
                                           true_parameter_values=c(0.5,0.3,0.001))
write.csv(GMM_case1_accuracy_100,"GMM_case1_accuracy_100.csv")


#case 2
options(warn=-1)

GMM_results_case2_100<-Simul_data_nruns_GMM_own(X0=2,b=2,d=1,c=0.01,ti=0,tmax=17, nrep=100,nruns=100)

write.csv(GMM_results_case2_100$GMM_estimates,"GMM_results_case2_100.csv")

#Saving accuracy measures
GMM_case2_accuracy_100<-Summary_peformance(MME_data=GMM_results_case2_100$GMM_estimates,
                                           true_parameter_values=c(2,1,0.01))
write.csv(GMM_case2_accuracy_100,"GMM_case2_accuracy_100.csv")


#case 3
options(warn=-1)

GMM_results_case3_100<-Simul_data_nruns_GMM_own(X0=2,b=3,d=2,c=0.1,ti=0,tmax=17, nrep=100,nruns=100)
write.csv(GMM_results_case3_100$GMM_estimates,"GMM_results_case3_100.csv")

Time_results2[,1]<- c(GMM_results_case1_100$Computational_time,GMM_results_case2_100$Computational_time,
                      GMM_results_case3_100$Computational_time)
write.csv(Time_results2,"Time_experiment2_GMM.csv")


#Saving accuracy measures
GMM_case3_accuracy_100<-Summary_peformance(MME_data=GMM_results_case3_100$GMM_estimates,
                                           true_parameter_values=c(3,2,0.1))
write.csv(GMM_case3_accuracy_100,"GMM_case3_accuracy_100.csv")







#Checking for consistency of Generalized Method of method estimator
#Increasing the number of simulation for one realization to 500 for 100 replicates**
Time_results3<- data.frame(matrix(NA,nrow=3,ncol=1))
names(Time_results3)<-c("GMM_time") 
rownames(Time_results3)<-c("case1","case2","case3")
#case 1
options(warn=-1)
GMM_results_case1_500<-Simul_data_nruns_GMM_own(X0=2,b=0.5,d=0.3,c=.001,ti=0,tmax=17, nrep=500,nruns=100)
write.csv(GMM_results_case1_500$GMM_estimates,"GMM_results_case1_500.csv")

#Saving accuracy measures
GMM_case1_accuracy_500<-Summary_peformance(MME_data=GMM_results_case1_500$GMM_estimates,
                                           true_parameter_values=c(0.5,0.3,0.001))
write.csv(GMM_case1_accuracy_500,"GMM_case1_accuracy_500.csv")


#case 2
options(warn=-1)

GMM_results_case2_500<-Simul_data_nruns_GMM_own(X0=2,b=2,d=1,c=0.01,ti=0,tmax=17, nrep=500,nruns=100)
write.csv(GMM_results_case2_500$GMM_estimates,"GMM_results_case2_500.csv")

#Saving accuracy measures
GMM_case2_accuracy_500<-Summary_peformance(MME_data=GMM_results_case2_500$GMM_estimates,
                                           true_parameter_values=c(2,1,0.01))
write.csv(GMM_case2_accuracy_500,"GMM_case2_accuracy_500.csv")


#case 3
options(warn=-1)
GMM_results_case3_500<-Simul_data_nruns_GMM_own(X0=2,b=3,d=2,c=0.1,ti=0,tmax=17, nrep=500,nruns=100)
write.csv(GMM_results_case3_500$GMM_estimates,"GMM_results_case3_500.csv")

Time_results3[,1]<- c(GMM_results_case1_500$Computational_time,GMM_results_case2_500$Computational_time,
                      GMM_results_case3_500$Computational_time)
write.csv(Time_results3,"Time_experiment3_GMM.csv")
#Saving accuracy measures
GMM_case3_accuracy_500<-Summary_peformance(MME_data=GMM_results_case3_500$GMM_estimates,
                                           true_parameter_values=c(3,2,0.1))
write.csv(GMM_case3_accuracy_500,"GMM_case3_accuracy_500.csv")
