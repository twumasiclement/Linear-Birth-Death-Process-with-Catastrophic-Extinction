####Setting the plot size and resolution (300 dpi)
options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

# Loading relevant R packages
library('latex2exp') #For adding LaTeX symbols to R plots
library(gmp) #Dealing with larger factorial limits
library(compiler)# byte code compilation
library(parallel)
library(data.table)
library(RColorBrewer)
library(JuliaCall) #for julia commands & scripts

# Number of cores to use for parallelization (adjust as needed)
num_cores <- detectCores() - 1  # Using one less core to avoid overloading
num_cores

# Set up the Julia environment
julia_setup()

# julia_command('readdir("/Users/clementtwumasi/Desktop/")')
# julia_command('readdir("/Users/clementtwumasi/Desktop/BDC_Paper_Journal_of_Applied_Probability /Main_R_codes")')

# Define the location of the BDCfit.jl script
julia_command('BDCdir = "/Users/clementtwumasi/Desktop/BDC_Paper_Journal_of_Applied_Probability /Main_R_codes/BDC_Loglik.jl"')

# Push the directory containing the script to LOAD_PATH
julia_command("push!(LOAD_PATH, dirname(BDCdir))")

# Include the script to make the module available
julia_command('include(BDCdir)')

# Now, use the module after it is included
julia_command("using .BDCfit")


# Estimating $\beta_{j}^{n}=\frac{\gamma_{j}^{n}}{n!}$ from loops 
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
PGF_z_compiler(lambda=0.513,mu=0.35,rho=0.003,t=1,z=1,m=2)

#Estimating C(t)=P(catastrophe resulting in 0 population)
Constant_added=function(lambda,mu,rho,t,z=1,m){
constant=1-PGF_z_compiler(lambda=lambda,mu=mu,rho=rho,t=t,z=z,m=m)
return(constant)
    }

PGF_z_correct=function(lambda,mu,rho,t,z,m){
   #v0<-((lambda+mu+rho)-sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
   #v1<-((lambda+mu+rho)+sqrt( ((lambda+mu+rho)^2)-4*mu*lambda))/(2*lambda)
    constants=BDCconsts(lambda,mu,rho,t)
    v0<- constants$v0
    v1<- constants$v1
    sigma<- constants$sigma
    num<-(v1*(1-sigma))+(z*(v1*sigma-v0))
    den<- v1-(sigma*v0)-(z*v0*(1-sigma))
    return( (num/den)^m)
}
PGF_z_correct_compiler=cmpfun(PGF_z_correct)
PGF_z_correct_compiler(lambda=0.513,mu=0.35,rho=0.003,t=1,z=1,m=2)


# The transition function for the BDC (with the Julia ProbBDC code)
Pn_m_t<-function(n,m,lambda, mu, rho,t){ 
 constants=BDCconsts(lambda=lambda, mu=mu, rho=rho,t=t)
    C_t=Constant_added(lambda=lambda, mu=mu, rho=rho,t=t,m=m)
    k1=constants$k1;k2=constants$k2;k3=constants$k3
    
    if(n==0 & m>=1){return((k1^m)+C_t)} # P(X(t) = 0 | X(0) = m) (due natural & catastrophe )
    if(n==0 & m==0){return(1)} #given 0 parasites the probability of having 0 parasites at any time t is 1
    if(n>0 & m==0){return(0)}  #given 0 parasites the probability of having  parasites n>0 at any time t is 0
  
    #If n>=1 and m>=1 compute the transition probability based on the Julia code
    #lambda, mu, rho, t, m=mmax, n=nmax
    if(n>=1 & m>=1){
     Prob<- tail(c(julia_call("ProbBDC", lambda, mu, rho, t,m, n)), n=1)
        return(Prob)
    }
    
    
}

Pn_m_t_compiler=cmpfun(Pn_m_t)

Pn_m_t_compiler(n=20,m=10,lambda=0.513, mu=0.35, rho=0.003,t=2)

o<-par(mar=c(0,4,2,2))#Run this before the code below

nf<-layout(matrix(1:4, nrow=2,ncol=2))


#First half of the plots (Plot1)
t=seq(0,100) # time
result_new=list()
params=c(0.512,0.35,0.003)
for (i in 1:7){
result_new[[i]]=numeric(length(t))
    }
#Given m=1
m=1

for (i in seq_along(t)){
result_new[[1]][i]= Pn_m_t_compiler(n=0,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[2]][i]=Pn_m_t_compiler(n=1,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[3]][i]= Pn_m_t_compiler(n=2,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[4]][i]=Pn_m_t_compiler(n=3,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))  
result_new[[5]][i]=Pn_m_t_compiler(n=4,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[6]][i]= Pn_m_t_compiler(n=5,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))
result_new[[7]][i]=Pn_m_t_compiler(n=10,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))
}

par(o)#
o<-par(mar=c(0,4,2,2)) 

colours=rainbow(20)[c(1,4,7,11,14,17,19)]    #brewer.pal(n=7, name="Dark2")
plot(t,result_new[[1]],type="l",col=colours[1],lwd=2.5,ylim=c(0,1.1),xaxt = "n",ylab="",xlab="",main="")
  
 for(i in 2:7){   
    lines(t,result_new[[i]],type="l",col=colours[i],lwd=2.5)
}

text(50,1.1, expression(paste("Parameter values:",lambda,"=0.512",";","", mu,"=0.35",";",rho,"=0.003") ),
    col="black",lwd=3)

text <- c("P{X(t)=0|X(0)=1}","P{X(t)=1|X(0)=1}","P{X(t)=2|X(0)=1}","P{X(t)=3|X(0)=1}",
         "P{X(t)=4|X(0)=1}","P{X(t)=5|X(0)=1}","P{X(t)=10|X(0)=1}")

legend(x=25,y=.75,legend = text,col=colours,cex=1, horiz = FALSE,pt.cex = 1,
       box.lwd = 2,fill=colours,ncol=1,title="Transition probability")


mtext(text = "Probability",
      side = 2,#side 2 = left
      line = 3,cex=.9)

#First half of the plots (Plot2)

t=seq(0,100) # time
result_new=list()
params=c(0.512,0.35,0.003)
for (i in 1:7){
result_new[[i]]=numeric(length(t))
    }
#Given m=2
m=2

for (i in seq_along(t)){
result_new[[1]][i]= Pn_m_t_compiler(n=0,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[2]][i]=Pn_m_t_compiler(n=1,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[3]][i]= Pn_m_t_compiler(n=2,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[4]][i]=Pn_m_t_compiler(n=3,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))  
result_new[[5]][i]=Pn_m_t_compiler(n=4,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[6]][i]= Pn_m_t_compiler(n=5,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))
result_new[[7]][i]=Pn_m_t_compiler(n=10,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))
}

par(mar=c(4,4,0,2))
colours=rainbow(20)[c(1,4,7,11,14,17,19)]    #brewer.pal(n=7, name="Dark2")
plot(t,result_new[[1]],type="l",col=colours[1],lwd=2.5,ylim=c(0,1.1),ylab="",xlab="Time (in days)",main="")
  
 for(i in 2:7){   
    lines(t,result_new[[i]],type="l",col=colours[i],lwd=2.5)
}

text(50,1.1, expression(paste("Parameter values:",lambda,"=0.512",";","", mu,"=0.35",";",rho,"=0.003") ),
    col="black",lwd=3)

text <- c("P{X(t)=0|X(0)=2}","P{X(t)=1|X(0)=2}","P{X(t)=2|X(0)=2}","P{X(t)=3|X(0)=2}",
         "P{X(t)=4|X(0)=2}","P{X(t)=5|X(0)=2}","P{X(t)=10|X(0)=2}")
#legend(par("usr")[1],par("usr")[3],text,col=colours,xjust=0, yjust=2.4,bty = "n",lwd=3,ncol=3)
legend(x=25,y=.75,legend = text,col=colours,cex=1, horiz = FALSE,pt.cex = 1,
       box.lwd = 2,fill=colours,ncol=1,title="Transition probability")


mtext(text = "Probability",
      side = 2,#side 2 = left
      line = 3,cex=.9)

#Second half of the plots (Plot1)
#par(mfrow=c(3,3),mar=c(4,4,1,1),xpd=NA,oma=c(2,0,0,0))

t=seq(0,100) # time
result_new=list()
params=c(0.512,0.35,0.003)
for (i in 1:7){
result_new[[i]]=numeric(length(t))
    }
#Given m=5
m=5

for (i in seq_along(t)){
result_new[[1]][i]= Pn_m_t_compiler(n=0,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[2]][i]=Pn_m_t_compiler(n=1,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[3]][i]= Pn_m_t_compiler(n=2,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[4]][i]=Pn_m_t_compiler(n=3,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))  
result_new[[5]][i]=Pn_m_t_compiler(n=4,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[6]][i]= Pn_m_t_compiler(n=5,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))
result_new[[7]][i]=Pn_m_t_compiler(n=10,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))
}

o<-par(mar=c(0,4,2,2)) 
colours=rainbow(20)[c(1,4,7,11,14,17,19)]    #brewer.pal(n=7, name="Dark2")
plot(t,result_new[[1]],type="l",col=colours[1],lwd=2.5,ylim=c(0,1.1),xaxt = "n",ylab="",xlab="",main="")
  
 for(i in 2:7){   
    lines(t,result_new[[i]],type="l",col=colours[i],lwd=2.5)
}

text(50,1.1, expression(paste("Parameter values:",lambda,"=0.512",";","", mu,"=0.35",";",rho,"=0.003") ),
    col="black",lwd=3)

text <- c("P{X(t)=0|X(0)=5}","P{X(t)=1|X(0)=5}","P{X(t)=2|X(0)=5}","P{X(t)=3|X(0)=5}",
         "P{X(t)=4|X(0)=5}","P{X(t)=5|X(0)=5}","P{X(t)=10|X(0)=5}")
#legend(par("usr")[1],par("usr")[3],text,col=colours,xjust=0, yjust=2.4,bty = "n",lwd=3,ncol=3)

legend(x=25,y=.75,legend = text,col=colours,cex=1, horiz = FALSE,pt.cex = 1,
       box.lwd = 2,fill=colours,ncol=1,title="Transition probability")



#Second half of the plots (Plot2)
#par(mfrow=c(3,3),mar=c(4,4,1,1),xpd=NA,oma=c(2,0,0,0))

t=seq(0,100) # time
result_new=list()
params=c(0.512,0.35,0.003)
for (i in 1:7){
result_new[[i]]=numeric(length(t))
    }
#Given m=10
m=10

for (i in seq_along(t)){
result_new[[1]][i]= Pn_m_t_compiler(n=0,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[2]][i]=Pn_m_t_compiler(n=1,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[3]][i]= Pn_m_t_compiler(n=2,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[4]][i]=Pn_m_t_compiler(n=3,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))  
result_new[[5]][i]=Pn_m_t_compiler(n=4,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i])) 
result_new[[6]][i]= Pn_m_t_compiler(n=5,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))
result_new[[7]][i]=Pn_m_t_compiler(n=10,m=m,lambda=as.numeric(params[1]), mu=as.numeric(params[2]), 
                                    rho=as.numeric(params[3]),t=as.numeric(t[i]))
}

par(mar=c(4,4,0,2))
colours=rainbow(20)[c(1,4,7,11,14,17,19)]    #brewer.pal(n=7, name="Dark2")
plot(t,result_new[[1]],type="l",col=colours[1],lwd=2.5,ylim=c(0,1.1),ylab="",xlab="Time (in days)",main="")
  
 for(i in 2:7){   
    lines(t,result_new[[i]],type="l",col=colours[i],lwd=2.5)
}

text(50,1.1, expression(paste("Parameter values:",lambda,"=0.512",";","", mu,"=0.35",";",rho,"=0.003") ),
    col="black",lwd=3)

text <- c("P{X(t)=0|X(0)=10}","P{X(t)=1|X(0)=10}","P{X(t)=2|X(0)=10}","P{X(t)=3|X(0)=10}",
         "P{X(t)=4|X(0)=10}","P{X(t)=5|X(0)=10}","P{X(t)=10|X(0)=10}")

legend(x=25,y=.75,legend = text,col=colours,cex=1, horiz = FALSE,pt.cex = 1,
       box.lwd = 2,fill=colours,ncol=1,title="Transition probability")

par(o)

#True mean from the pgf
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

#Computing the transition probabilities at m=2
n_values=0:1000 #possible values of n=0,1,2,3, .... (finite sum)
t=seq(1,30,length.out=30) # time
 Prob_distribution=list()
time0=proc.time()
for(i in seq_along(n_values)){Prob_distribution[[i]]=numeric(length(t))} #Creating list of vectors of 0's 
   for (i in t){
       for(k in seq_along(n_values)){
       Prob_distribution[[k]][i]=Pn_m_t_compiler(n=as.numeric(n_values[k]),
                                            m=2,lambda=0.512, mu=0.35, rho=0.003,t=i) 
                  }  
   }

proc.time()-time0

# setting working directory
setwd("/Users/clementtwumasi/Desktop/BDC_Paper_Journal_of_Applied_Probability /Transition_Proof_saved_results")

#A function to convert NAs to 0
na.zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
}



#Combining transition probabilities from n=1 to n=1000
Analytical_probabilities<-data.frame(matrix(data=0,ncol =length(n_values),nrow = 30 )) #nsims x timepoints

for(i in seq_along(n_values)) names(Analytical_probabilities)[i]<-paste("n=",n_values[i])

Analytical_means<-numeric(length(t))
for (k in seq_along(n_values)){
Analytical_probabilities[,k]<-as.data.frame(na.zero(Prob_distribution[[k]]))
}

for (i in 1:30){
Analytical_means[i]=sum((n_values[seq_along(n_values)])*Analytical_probabilities[i,])
}

# exporting the analytical probabililities
write.csv(Analytical_probabilities, "Analytical_probabilities.csv")

write.csv(Analytical_means,"Analytical_means.csv")

Pn_m_t_compiler(n=0,m=1,lambda=0.512, mu=0.35, rho=0.003,t=0)

t_values=0:100

for(k in seq_along(t_values))
   print(Pn_m_t_compiler(n=0,m=3,lambda=0.512, mu=0.35, rho=0.003,t=as.numeric(t_values[k])))


#Importing the transition probabilities based on the derived transition function
n_values=0:1000
Analytical_probs_data<- read.csv(file="Analytical_probabilities.csv")
Analytical_probs<-Analytical_probs_data[,-1]

for(i in seq_along(n_values)) names(Analytical_probs)[i]<-paste("n=",n_values[i])
dim(Analytical_probs)
Analytical_probs[, 1:9] #viewing first 9 columns (rows represent time=1 to 30)

#Importing the expectation E(X) based on the transition probabilities
Analytical_mean<- read.csv(file="Analytical_means.csv")
Analytical_means<-Analytical_mean[[2]]
Analytical_means

Analytical_variance<-numeric(length(t))
for (i in 1:30){
Analytical_variance[i]= sum(((n_values[seq_along(n_values)])^2)*Analytical_probabilities[i,])-
    sum((n_values[seq_along(n_values)])*Analytical_probabilities[i,])^2
}

write.csv(Analytical_variance,"Analytical_variance.csv")

#Importing the expectation E(X) based on the transition probabilities
Analytical_var<- read.csv(file="Analytical_variance.csv")
Analytical_variance<-Analytical_var[[2]]
Analytical_variance

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

t=seq(1,30) # time
True_exact_mean=Expectation_exact(b=0.512, d=0.35, c=0.003,t,m=2)
True_exact_mean

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


t=seq(1,30) # time
True_exact_variance=Variance_analytic(b=0.512, d=0.35, c=0.003,t,m=2)
True_exact_variance

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
Exact_BDC_compiler(X0=2,b=.512,d=.35,c=.003,ti=0,tmax=17)

ExactSSA_sim <- function(X0,b,d,c,ti=0,tmax=30, nrep=1e1) {


  f <- function(x) Exact_BDC_compiler(X0=X0,b=b,d=d,c=c,ti=ti,tmax=tmax)$pop
  x <- lapply(1:nrep, f)
  pop<-t(do.call("rbind",x)) #rbinding simulated data on time by nsim
      return(pop)
}

time0=proc.time()
Sim_combine<-ExactSSA_sim(X0=2,b=0.512,d=0.35,c=0.003,ti=0,tmax=30,nrep=1e6) 
proc.time()-time0

#A function to convert NAs to 0
na.zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
}

Mean_func=function(x){return(mean(x,na.rm=T))}
Variance_func=function(x){return(var(x,na.rm=T))}

Expectation_Simulation_Main<-apply(na.zero(Sim_combine),1,Mean_func)
round(Expectation_Simulation_Main,4)

write.csv(Expectation_Simulation_Main,"Expectation_Simulation_Main.csv")

# importing Expectation_Simulation_Main.csv
Expectation_Simulation_Main<-read.csv(file="Expectation_Simulation_Main.csv")
Expectation_Simulation_Main<-Expectation_Simulation_Main[[2]]
Expectation_Simulation_Main

Variance_Simulation_Main<-apply(na.zero(Sim_combine),1,Variance_func)
round(Variance_Simulation_Main,4)

write.csv(Variance_Simulation_Main,"Variance_Simulation_Main.csv")

Variance_Simulation_Main<-read.csv(file="Variance_Simulation_Main.csv")
Variance_Simulation_Main<-Variance_Simulation_Main[[2]]
Variance_Simulation_Main

Conf_int95<- function(mean_data,var_data,nsim=1e6){
    lower<- mean_data-(1.96* sqrt(var_data/nsim))
   upper<- mean_data+(1.96*sqrt(var_data/nsim))
    lower_var<- var_data-(1.96* var_data*sqrt(2/(nsim-1)))
    upper_var<- var_data+(1.96* var_data*sqrt(2/(nsim-1)))
    return(list(lower=lower,upper=upper,lower_var=lower_var,upper_var=upper_var))
}

Conf_interval<- Conf_int95(mean_data=Expectation_Simulation_Main,var_data=Variance_Simulation_Main)

Conf_lower95<-Conf_interval$lower
Conf_lower95

Conf_upper95<-Conf_interval$upper
Conf_upper95


Conf_lower95_var<-Conf_interval$lower_var
Conf_lower95_var

Conf_upper95_var<-Conf_interval$upper_var
Conf_upper95_var

options(warn=-1)
#par(xpd=NA,oma=c(2,0,0,0))
o<-par(mar=c(0,4,2,2))#Run this before the code below

nf<-layout(matrix(1:2, nrow=2,ncol=1))

#first plot
o<-par(mar=c(0,4,2,2))

plot(t,Expectation_Simulation_Main,type="b",col="green",lwd=3,ylim=c(1.5,9),ylab="Mean",
     xlab="",xaxt="n")
lines(t,True_exact_mean,type="l",col="black",lwd=3)
lines(t,Analytical_means,type="l",col="red",lwd=3)#from transition function

# Add error bars
arrows(x0=t, y0=Conf_lower95, x1=t, y1=Conf_upper95, code=3, 
     angle=90, length=0.1,col ="blue",lwd=1)

text(8,9,
expression(paste("Parameters values:",lambda,"=0.512",",",mu,"=0.35",",",rho,"=0.003") ),col="blue",lwd=3)

text(5,8.5,TeX('Initial population: $X_0=2$') ,col="blue",lwd=3)

text <- c("Based on 1million exact simulations","Based on true mean of the B-D-C model","Based on the derived transition function")
#legend(par("usr")[1],par("usr")[3],text,col=c("green","black","red"),xjust=0, yjust=1.4,bty = "n",lwd=3,ncol=1)

legend(x=8,y=5,legend = text,col=c("green","black","red")
      ,text.width = strwidth(text)[1]*.78, cex=0.8,box.lwd = 2,title="Mean comparison",
      fill=c("green","black","red"))


#second plot
par(mar=c(4,4,0,2))

#par(xpd=NA,oma=c(2,0,0,0))

plot(t,Variance_Simulation_Main,type="b",col="green",lwd=3,ylab="Variance",xlab="Time (in days)")
lines(t,True_exact_variance ,type="l",col="black",lwd=3)
lines(t,Analytical_variance,type="l",col="red",lwd=3) #from transition function

# Add error bars
arrows(x0=t, y0=Conf_lower95_var, x1=t, y1=Conf_upper95_var, code=3, 
     angle=90, length=0.1,col ="blue",lwd=1)

text <- c("Based on 1million exact simulations","Based on true variance of the B-D-C model","Based on the derived transition function")

legend(par("usr")[1],par("usr")[3],text,col=c("green","black","red"),xjust=0, yjust=1.4,bty = "n",lwd=3,ncol=1)
#text(8,520,
#expression(paste("Parameters values:",lambda,"=0.512",",",mu,"=0.35",",",rho,"=0.003") ),col="blue",lwd=3)

#text(5,480,TeX('Initial population: $X_0=2$') ,col="blue",lwd=3)

legend(x=13,y=200,legend = text,col=c("green","black","red")
      ,text.width = strwidth(text)[1]*.83, cex=0.8,box.lwd = 2,title="Variance comparison",
      fill=c("green","black","red"))

# Graphical comparison of means
#par(xpd=NA,oma=c(2,0,0,0))
o<-par(mar=c(0,4,2,2))#Run this before the code below

nf<-layout(matrix(1:2, nrow=2,ncol=1))

#first plot
o<-par(mar=c(0,4,2,2))

plot(t,Expectation_Simulation_Main,type="l",col="green",lwd=3,ylim=c(1.5,10),ylab="Mean",
     xlab="",xaxt="n")
lines(t,True_exact_mean,type="l",col="red",lwd=3)
#lines(t,Analytical_means,type="l",col="red",lwd=3)#from transition function

# text(8,10,
# expression(paste("Parameters values:",lambda,"=0.512",",",mu,"=0.35",",",rho,"=0.003") ),
     # col="black",lwd=3,font=4)

# text(5.5,9.5,TeX('Initial population size: $X_0=2$') ,col="black",lwd=3, font=3)

text <- c("Based on 1 million exact simulations","Based on true mean of the B-D-C process")

#legend(par("usr")[1],par("usr")[3],text,col=c("green","black","red"),xjust=0, yjust=1.4,bty = "n",lwd=3,ncol=1)

legend(x=8,y=5,legend = text,col=c("green","red")
      ,text.width = strwidth(text)[1]*.8, cex=0.8,box.lwd = 1.2,title="Mean comparison",
      fill=c("green","red"))


#second plot
par(mar=c(4,4,0,2))

#par(xpd=NA,oma=c(2,0,0,0))

plot(t,Variance_Simulation_Main,type="l",col="green",lwd=3,ylab="Variance",xlab="Time (in days)")
lines(t,True_exact_variance ,type="l",col="red",lwd=3)
#lines(t,Analytical_variance,type="l",col="red",lwd=3) #from transition function


text <- c("Based on 1million exact simulations","Based on true variance of the B-D-C process")
#legend(par("usr")[1],par("usr")[3],text,col=c("red","black"),xjust=0, yjust=1.4,bty = "n",lwd=3,ncol=1)
#text(8,520,
#expression(paste("Parameters values:",lambda,"=0.512",",",mu,"=0.35",",",rho,"=0.003") ),col="blue",lwd=3)

#text(5,480,TeX('Initial population: $X_0=2$') ,col="blue",lwd=3)

legend(x=13,y=200,legend = text,col=c("green","red")
      ,text.width = strwidth(text)[1]*.85, cex=0.8,box.lwd = 1.2,title="Variance comparison",
      fill=c("green","red"))

# Graphical comparison of variances
par(xpd=NA,oma=c(2,0,0,0))

plot(t,Variance_Simulation_Main,type="l",col="green",lwd=3,ylab="Variance",xlab="Time (in days)")
lines(t,True_exact_variance ,type="l",col="black",lwd=3)
lines(t,Analytical_variance,type="l",col="red",lwd=3) #from transition function


text <- c("Based on 1million exact simulations","Based on true variance of the B-D-C model","Based on the derived transition function")

#legend(par("usr")[1],par("usr")[3],text,col=c("green","black","red"),xjust=0, yjust=1.4,bty = "n",lwd=3,ncol=1)
text(8,520,
expression(paste("Parameters values:",lambda,"=0.512",",",mu,"=0.35",",",rho,"=0.003") ),col="blue",lwd=3)

text(5,480,TeX('Initial population: $X_0=2$') ,col="blue",lwd=3)

legend(x=13,y=200,legend = text,col=c("green","black","red")
      ,text.width = strwidth(text)[1]*.8, cex=0.8,box.lwd = 2,title="Variance comparison",
      fill=c("green","black","red"))

# Comparing probability distribution between the Analytical transition function and Monte Carlo estimates from the exact simulation
t<- 1:30
nsims=1e6
breaks=vector("list",length=length(t))
Actual_probabilities=vector("list",length=length(t))
counts=vector("list",length=length(t))
Actual_probabilities=vector("list",length=length(t))
Data_probabilities=vector("list",length=length(t))
MonteCarlo_estimates=vector("list",length=length(t))
span_cut=vector("list",length=length(t))

for(time in t){

breaks[[time]]=seq(min(na.zero(Sim_combine[time,])),max(na.zero(Sim_combine[time,]))+1,by=1)
Actual_probabilities[[time]]=numeric(length(breaks[[time]][-length(breaks[[time]])]))
counts[[time]]=breaks[[time]][-length(breaks[[time]])]

for (i in 1:length(counts[[time]])){
Actual_probabilities[[time]][i]=Pn_m_t_compiler(n=as.numeric(counts[[time]][i]),
                            m=2,lambda=0.512, mu=0.35, rho=0.003,t=as.numeric(time))
}
    
    
Data_probabilities[[time]]=as.data.table(matrix(0,nrow=(length(breaks[[time]])-1), ncol=3,byrow=T))
names(Data_probabilities[[time]])=c("Population","Actual","Estimate")
span_cut[[time]]=cut(as.numeric(na.zero(Sim_combine[time,])), 
                     breaks[[time]], len = length(breaks[[time]])+1,right=F)
MonteCarlo_estimates[[time]]=table(span_cut[[time]])/nsims
Data_probabilities[[time]][,1]=breaks[[time]][-length(breaks[[time]])] #unique values of population count
Data_probabilities[[time]][,2]=Actual_probabilities[[time]]
Data_probabilities[[time]][,3]=MonteCarlo_estimates[[time]]

    }

time=1
head(Data_probabilities[[time]])

Conf_interval95=function(p_hat){
    lower=p_hat-(1.96*sqrt((p_hat*(1-p_hat))/nsims))
    upper=p_hat+(1.96*sqrt((p_hat*(1-p_hat))/nsims))
    return(list(lower_limit=lower,upper_limit=upper))
}

#Adding lower and upper bounds of  Monte Carlo estimates
x=y_est=y_act=y.sd_lower=y.sd_upper=vector("list",length=length(t))
for (i in t){
Data_probabilities[[i]]$lower_bound=Conf_interval95(p_hat=Data_probabilities[[i]]$Estimate)$lower_limit
Data_probabilities[[i]]$upper_bound=Conf_interval95(p_hat=Data_probabilities[[i]]$Estimate)$upper_limit
x[[i]]=as.vector(Data_probabilities[[i]][[1]])
y_est[[i]]=as.vector(Data_probabilities[[i]][[3]])
y_act[[i]]=as.vector(Data_probabilities[[i]][[2]])
y.sd_lower[[i]]=as.vector(Data_probabilities[[i]][[4]])
y.sd_upper[[i]]=as.vector(Data_probabilities[[i]][[5]])
}

#Adding an indicator to indicate whether the actual values falls in the confidence interval
for(time in t){
  for (i in 1:dim(Data_probabilities[[time]])[1]){
    
 Data_probabilities[[time]]$Within_CI<-Data_probabilities[[time]]$Actual>Data_probabilities[[time]]$lower_bound&
    Data_probabilities[[time]]$Actual<Data_probabilities[[time]]$upper_bound
    
}
    }

time=1
head( Data_probabilities[[time]])

Within_CI<-rep(NA,length=length(t))
Within_CI[1]=1
for(time in t[-1]){
Within_CI[time]<-as.numeric(table(Data_probabilities[[time]]$Within_CI)/
                              sum(table(Data_probabilities[[time]]$Within_CI)))[2]
       }

unlist(Within_CI)


cat("The proportion of times the theoritical estimate was found within the 
    estimated confidence interval is:","",
      round(mean(Within_CI[1:6]),3),"+-","",round(sd(Within_CI[1:6]),3))


c(0.916-0.052,0.916+0.052)

Combined_data<-rbind(Data_probabilities[[1]],Data_probabilities[[2]],Data_probabilities[[3]],
                     Data_probabilities[[4]],Data_probabilities[[5]],Data_probabilities[[6]])





print(
paste("The proportion of times the theoritical estimate was found within the estimated C.I is","",
      round(as.numeric(table(Combined_data$Within_CI)/
                              sum(table(Combined_data$Within_CI)))[2],3)))


write.csv(Combined_data, "Combined_data.probabilities.csv")

par(mfrow=c(3,2),mar=c(4,4,2,1))
x_position<- c(18,30,36,52,62,80)
for(time in 1:6){
# Plot
plot(x[[time]],y_est[[time]], xlab="Population size", ylab="Probability estimate",type="b", pch=3, 
     cex=.7,col="blue",lwd=1,main="",ylim=c(0,0.4),las=1)

lines(x[[time]],y_act[[time]], type="p",pch=17, cex=.7,col="red",lwd=3)

text <- c("Monte Carlo estimates","Actual values")
#legend(par("usr")[1],par("usr")[3],text,col=c("green","black","red"),xjust=0, yjust=1.4,bty = "n",lwd=3,ncol=1)
#mtext("Population",side=3,col="black",line=1.5,cex=.8)  

    if(time==3){
legend(x=10,y=0.3,legend = text,col=c("blue","red"),box.lwd = 1.2
      ,text.width = strwidth(text)[1]*1, cex=1.1,title="Transition function",pch=c(3,17))
    }

    if(time==4){
legend(x=20,y=0.3,legend = text,col=c("blue","red"),box.lwd = 1.2
      ,text.width = strwidth(text)[1]*1, cex=1.1,title="Transition function",pch=c(3,17))
    }

    
text(x_position[time],0.39,paste("t =",time),font=2,cex=1.2)

    
    }
