library(gmp) #Dealing with larger factorial limits with arbitrary precision
library(compiler)# byte code compilation
library(data.table)
library(parallel)
library(JuliaCall) #for julia commands & scripts

library("maxLik")#for maximum likelihood estimation/optimization
#library('latex2exp') #For adding LaTeX symbols to R plots

setwd("/home/clement/Clement_R_Results_folder/GMM_MLE_results_time")

#options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

# Calculate the number of cores for parLapply
#no_cores <- min(16, detectCores())

# Initiate cluster
#cluster <- makeCluster(no_cores)
#cluster
ncores=16
ncores


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



ExactSSA_sim <- function(X0,b,d,c,ti=0,tmax=30, nrep=1e1) {
  
  
  f <- function(x) Exact_BDC_compiler(X0=X0,b=b,d=d,c=c,ti=ti,tmax=tmax)$pop
  x <- lapply(1:nrep, f)
  pop<-t(do.call("rbind",x)) #rbinding simulated data on nrep fish (each column for a fish)
  return(pop)
}



#Simulated 100 realization of data of 50 fish with MLE
#Using Julia log-likelihood function

Simul_data_nruns_MLE<- function(X0,b,d,c,ti=0,tmax=30, nrep=1e1,nruns=2, nc=ncores) {
  
  Sim_data<-array(dim = c(nruns,9,nrep))
  
  x1_out <- mclapply(1:nruns,
                     function(x) ExactSSA_sim(X0=X0,b=b,d=d,c=c,ti=ti,tmax=tmax,nrep=nrep), mc.cores=nc)
  
  
  for (k in 1:nruns){
    for(i in 1:nrep){
      # saving results for each 100 realizations and 
      Sim_data[k, ,i]<- x1_out[[k]][, i][seq(1,17,by=2)] 
    }
    
  }
  
  
  #Computing MLE based on Julia function 
  time0=proc.time()
  MLE_estimates_Julia <- lapply(1:nruns, 
             function(k) MLE_BDC_Julia_compiler(Parasite_numbers=as.data.frame(Sim_data[k, ,])))  #as lists
  time1=proc.time()-time0
  
  LogLike_MLE_Julia <- as.data.frame(do.call("rbind",MLE_estimates_Julia)) 
  colnames(LogLike_MLE_Julia)<-c("b_MLE","d_MLE","c_MLE")
  
  
  #returns MLE estimates for each realization  
   
  time<-sum(as.vector(time1)[-3])
  return(list(LogLike_MLE_Julia=LogLike_MLE_Julia,Computational_time=time))
  
}

Summary_peformance<-function(MLE_data,true_parameter_values){
  #Computing the mean estimates
  means_estimates<-apply(MLE_data,2,mean)
  ML_estimates<- as.vector(means_estimates)
  
  #computing the variance
  var_estimates<-apply(MLE_data,2,var)
  var_MLE<- as.vector(var_estimates)
  
  #Computing bias for each parameter
  bias_birth<-  ML_estimates[1]-true_parameter_values[1]
  bias_death<-  ML_estimates[2]-true_parameter_values[2] 
  bias_catastrophe<-  ML_estimates[3]-true_parameter_values[3]
  bias_MLE<-c(bias_birth, bias_death,bias_catastrophe)
  #Total_bias<- bias_birth+ bias_death+bias_catastrophe
  
  
  #Computing MSE
  #for R
  MSE_birth<- var_MLE[1]+bias_birth^2
  MSE_death<- var_MLE[2]+bias_death^2
  MSE_catastrophe<- var_MLE[3]+bias_catastrophe^2
  MSE_MLE<-c(MSE_birth,MSE_death, MSE_catastrophe)
  
  #Total MSE
  #for R
  #MSE_total_R<- MSE_birth_R+MSE_death_R+MSE_catastrophe_R
  #Computing confidence interval
  n<- dim(MLE_data)[1]
  number_of_parameters<- dim(MLE_data)[2]
  MLE_mean_estimates=stderr=Lower_lim95=Upper_lim95=rep(NA,number_of_parameters) 
  for (i in 1:number_of_parameters){
    #Expected ML estimates       
    MLE_mean_estimates[i]<-   ML_estimates[i]
    
    #std error= sqrt(Var/n)  
    stderr[i]<- sqrt(var_MLE[i]/n)
    
    #Lower limits  
    Lower_lim95[i]<- MLE_mean_estimates[i]- (1.96*stderr[i])
    
    #Upper limits  
    Upper_lim95[i]<- MLE_mean_estimates[i]+ (1.96*stderr[i])
  }
  
  
  output_accuracy<- data.frame(ML_estimates,true_parameter_values,bias_MLE,var_MLE, MSE_MLE)
  rownames(output_accuracy)<- c("birth rate","death_rate","catastrophe_rate")
  
  output_accuracy$Lower_lim95<-Lower_lim95
  output_accuracy$Upper_lim95<-Upper_lim95
  
  return(output_accuracy)
}

Time_results1_MLE<- data.frame(matrix(NA,nrow=3,ncol=1))
names(Time_results1_MLE)<-c("MLE_time") 
rownames(Time_results1_MLE)<-c("case1","case2","case3")

Output1=Simul_data_nruns_MLE(X0=2,b=0.5,d=0.3,c=0.001,ti=0,tmax=17,nrep=50,nruns=100)

write.csv(Output1$LogLike_MLE_Julia,"MLE_Case1_Results.csv")



Accuracy_MLE_Case1<- Summary_peformance(MLE_data=Output1$LogLike_MLE_Julia,
                                        true_parameter_values=c(0.5,0.3,0.001))
write.csv(Accuracy_MLE_Case1,"Accuracy_MLE_Case1.csv")



Output2=Simul_data_nruns_MLE(X0=2,b=2,d=1,c=0.01,ti=0,tmax=17,nrep=50,nruns=100)

write.csv(Output2$LogLike_MLE_Julia,"MLE_Case2_Results.csv")




Accuracy_MLE_Case2<-Summary_peformance(MLE_data=Output2$LogLike_MLE_Julia,
                                       true_parameter_values=c(2,1,0.01))
write.csv(Accuracy_MLE_Case2,"Accuracy_MLE_Case2.csv")



Output3=Simul_data_nruns_MLE(X0=2,b=3,d=2,c=0.1,ti=0,tmax=17,nrep=50,nruns=100)

write.csv(Output3$LogLike_MLE_Julia,"MLE_Case3_Results.csv")


Time_results1_MLE[,1]<- c(Output1$Computational_time,Output2$Computational_time,
                          Output3$Computational_time)
write.csv(Time_results1_MLE,"Time_experiment1_MLE.csv")

Accuracy_MLE_Case3<-Summary_peformance(MLE_data=Output3$LogLike_MLE_Julia,
                                       true_parameter_values=c(3,2,0.1))
write.csv(Accuracy_MLE_Case3,"Accuracy_MLE_Case3.csv")
