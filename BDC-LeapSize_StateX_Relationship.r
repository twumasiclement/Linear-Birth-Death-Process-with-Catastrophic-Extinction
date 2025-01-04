####Setting the plot size and resolution (300 dpi)
options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

#Based on Gillespie 2001
tau_TotalRate_2001<-function(X,b,d,c,error){
    rate1<- b*X
    #probability of death
    rate2<- d*X
    #probability of catastrophe
    rate3<- c*X
    
    total_rate<- rate1+rate2+rate3 #representing a0(x)
    
     tau<-(error*(b+d))/(abs(b-d)*max(b,d)) 
    leap_condition<- 2/total_rate
    return(list(tau=tau,leap_conditions=leap_condition))
}

#Based on Gillespie & Petzold 2003
tau_TotalRate_2003<-function(X,b,d,c,error){  
   rate1<- b*X
    #probability of death
    rate2<- d*X
    #probability of catastrophe
    rate3<- c*X
    
    total_rate<- rate1+rate2+rate3 #representing a0(x)
     #Computing tau from equation 6 based on Gillespie & Petzold 2003
    Leap_sizes1= (error*(b+d))/(abs(b^2-b*d))
    Leap_sizes2= (error*(b+d))/(abs(b*d-d^2))
    Leap_sizes3= X*(error*(b+d))^2/(b^3+ d*b^2)
    Leap_sizes4= X*(error*(b+d))^2/(b*d^2+ d^3)
    
    #estimate of the tau leap size
    tau=min(Leap_sizes1,Leap_sizes2,Leap_sizes3,Leap_sizes4)  
     leap_condition<- 1/(10*total_rate)
    return(list(tau=tau,leap_condition= leap_condition))
}

#Based on Gillespie & Petzold 2003
tau_TotalRate_2003_2<-function(X,b,d,c,error){  
   rate1<- b*X
    #probability of death
    rate2<- d*X
    #probability of catastrophe
    rate3<- c*X
    
    total_rate<- rate1+rate2+rate3 #representing a0(x)
     #Computing tau from equation 6 based on Gillespie & Petzold 2003
    Leap_sizes1= (error*(b+d))/ (abs(b-d)*max(b,d))
    Leap_sizes2= X*(error*(b+d))^2/((b+d)*max(b^2,d^2))
    
    #estimate of the tau leap size
    tau=min(Leap_sizes1,Leap_sizes2)  
     leap_condition<- 1/(10*total_rate)
    return(list(tau=tau,leap_condition= leap_condition))
}

state<-seq(100,300,by=10)

 for (i in seq_along(state)) print(tau_TotalRate_2003(X=state[i],b=0.5,d=0.3,c=0.001,error=0.01)$tau)

state<- seq(100,300,by=10)

for (i in seq_along(state)) print(tau_TotalRate_2003_2(X=state[i],b=0.5,d=0.3,c=0.001,error=0.01)$tau)


tau2001=tau2003=results2001=results2003=leapCondition2001=leapCondition2003=NULL

state<- seq(1,1000,by=10)
tau2001[["case1"]]=leapCondition2001[["case1"]]=numeric(length(state))
tau2001[["case2"]]=leapCondition2001[["case2"]]=numeric(length(state))
tau2001[["case3"]]=leapCondition2001[["case3"]]=numeric(length(state))

for(i in seq_along(state)) {
results2001[["case1"]]<- tau_TotalRate_2001(X=state[i],b=0.5,d=0.3,c=0.001,error=0.01)
results2001[["case2"]]<- tau_TotalRate_2001(X=state[i],b=2,d=1,c=0.01,error=0.01)  
results2001[["case3"]]<- tau_TotalRate_2001(X=state[i],b=3,d=2,c=0.1,error=0.01)
    
tau2001[["case1"]][i]<-results2001[["case1"]]$tau
tau2001[["case2"]][i]<-results2001[["case2"]]$tau
tau2001[["case3"]][i]<-results2001[["case3"]]$tau
    
leapCondition2001[["case1"]][i]<-results2001[["case1"]]$leap_condition
leapCondition2001[["case2"]][i]<-results2001[["case2"]]$leap_condition
leapCondition2001[["case3"]][i]<-results2001[["case3"]]$leap_condition
}

state<- seq(1,1000,by=10)
tau2003[["case1"]]=leapCondition2003[["case1"]]=numeric(length(state))
tau2003[["case2"]]=leapCondition2003[["case2"]]=numeric(length(state))
tau2003[["case3"]]=leapCondition2003[["case3"]]=numeric(length(state))

for(i in seq_along(state)) {
results2003[["case1"]]<- tau_TotalRate_2003(X=state[i],b=0.5,d=0.3,c=0.001,error=0.01)
results2003[["case2"]]<- tau_TotalRate_2003(X=state[i],b=2,d=1,c=0.01,error=0.01)  
results2003[["case3"]]<- tau_TotalRate_2003(X=state[i],b=3,d=2,c=0.1,error=0.01)
    
tau2003[["case1"]][i]<-results2003[["case1"]]$tau
tau2003[["case2"]][i]<-results2003[["case2"]]$tau
tau2003[["case3"]][i]<-results2003[["case3"]]$tau
    
leapCondition2003[["case1"]][i]<-results2003[["case1"]]$leap_condition
leapCondition2003[["case2"]][i]<-results2003[["case2"]]$leap_condition
leapCondition2003[["case3"]][i]<-results2003[["case3"]]$leap_condition
}

# setting working directory
setwd("/Users/clementtwumasi/Desktop/BDC_Paper_Journal_of_Applied_Probability /LeapSize_conditions_results")

data=NULL
data[["case1"]]= data.frame(state=state,tau2001=tau2001[["case1"]],
                  leapCondition2001=leapCondition2001[["case1"]],tau2003=tau2003[["case1"]],
                 leapCondition2003=leapCondition2003[["case1"]])

data[["case2"]]= data.frame(state=state,tau2001=tau2001[["case2"]],
                  leapCondition2001=leapCondition2001[["case2"]],tau2003=tau2003[["case2"]],
                 leapCondition2003=leapCondition2003[["case2"]])

data[["case3"]]= data.frame(state=state,tau2001=tau2001[["case3"]],
                  leapCondition2001=leapCondition2001[["case3"]],tau2003=tau2003[["case3"]],
                 leapCondition2003=leapCondition2003[["case3"]])

write.csv(data[["case1"]],"LeapSize_State_Condition_Results1.csv")
write.csv(data[["case2"]],"LeapSize_State_Condition_Results2.csv")
write.csv(data[["case3"]],"LeapSize_State_Condition_Results3.csv")

state[which(tau2001$case1==tau2003$case1)]

state[which(tau2001$case2==tau2003$case2)]

#For tau-leaping (For Gillespie) X>=41 and no leaping for X<=31
#state[which(tau2001>leapCondition2001)]
paste("For Gillespie 2001, leaping starts when state x>=",min(state[which(tau2001[["case1"]]>leapCondition2001[["case1"]])]),
      
      "","where the fixed leap size=","",tau2001[["case1"]][which(tau2001[["case1"]]>leapCondition2001[["case1"]])[1]],"",
     "and leap condition >=",leapCondition2001[["case1"]][which(tau2001[["case1"]]>leapCondition2001[["case1"]])[1]])

#0.00672	0.0059449498
#For tau-leaping (For Gillespie) X>=21 and no leaping for X<=11
#state[which(tau2003>leapCondition2003)]

paste("For Gillespie & Petzold 2003, leaping starts when state x>=",
      min(state[which(tau2003[["case1"]]>leapCondition2003[["case1"]])]),
      
      "","where the leap size >=","",tau2003[["case1"]][[which(tau2003[["case1"]]>leapCondition2003[["case1"]])[1]]],"",
     "and leap condition >=",leapCondition2003[["case1"]][[which(tau2003[["case1"]]>leapCondition2003[["case1"]])[1]]])

#For tau-leaping (For Gillespie) X>=51 and no leaping for X<=41
#state[which(tau2001>leapCondition2001)]
paste("For Gillespie 2001, leaping starts when state x>=",min(state[which(tau2001[["case2"]]>leapCondition2001[["case2"]])]),
      
      "","where the fixed leap size=","",tau2001[["case2"]][which(tau2001[["case2"]]>leapCondition2001[["case2"]])[1]],"",
     "and leap condition >=",leapCondition2001[["case2"]][which(tau2001[["case2"]]>leapCondition2001[["case2"]])[1]])


#For tau-leaping (For Gillespie) X>=31 and no leaping for X<=21
#state[which(tau2003>leapCondition2003)]

paste("For Gillespie & Petzold 2003, leaping starts when state x>=",
      min(state[which(tau2003[["case2"]]>leapCondition2003[["case2"]])]),
      
      "","where the leap size >=","",tau2003[["case2"]][[which(tau2003[["case2"]]>leapCondition2003[["case2"]])[1]]],"",
     "and leap condition >=",leapCondition2003[["case2"]][[which(tau2003[["case2"]]>leapCondition2003[["case2"]])[1]]])

#For tau-leaping (For Gillespie) X>=31 and no leaping for X<=21
#state[which(tau2001>leapCondition2001)]
paste("For Gillespie 2001, leaping starts when state x>=",min(state[which(tau2001[["case3"]]>leapCondition2001[["case3"]])]),
      
      "","where the fixed leap size=","",tau2001[["case3"]][which(tau2001[["case3"]]>leapCondition2001[["case3"]])[1]],"",
     "and leap condition >=",leapCondition2001[["case3"]][which(tau2001[["case3"]]>leapCondition2001[["case3"]])[1]])


#For tau-leaping (For Gillespie) X>=21 and no leaping for X<=11
#state[which(tau2003>leapCondition2003)]

paste("For Gillespie & Petzold 2003, leaping starts when state x>=",
      min(state[which(tau2003[["case3"]]>leapCondition2003[["case3"]])]),
      
      "","where the leap size >=","",tau2003[["case3"]][[which(tau2003[["case3"]]>leapCondition2003[["case3"]])[1]]],"",
     "and leap condition >=",leapCondition2003[["case3"]][[which(tau2003[["case3"]]>leapCondition2003[["case3"]])[1]]])



## add extra space to right margin of plot within frame
#par(mar=c(5, 4, 4, 6) + 0.1)

o<-par(mar=c(0.5, 4, 0, 6) + 0.1)
#o<-par(mar=c(0,4,2,2))#Run this before the code below

nf<-layout(matrix(1:3, nrow=3,ncol=1))


#### FIRST PLOT ####
## Plot first set of data and draw its axis

plot(state,tau2001[["case1"]],type="l",col="blue",lwd=2,ylab="",xlab="",
    ylim=c(0,.1),axes=FALSE)

lines(state,tau2003[["case1"]],col="red",lwd=2,ylab="",xlab="")
abline(v=41,col="black",lty=3,lwd=1.5)
abline(v=21,col="grey",lty=3,lwd=1.5)


text1<- c("For HTL2001, leaping started at state x>40",
         "For HTL2003, leaping started at state x>20")
legend(x = 200, y = 0.05, legend = text1,
       col = c("black", "grey"), lty = 3, lwd = 1.5, bty = "n", cex = 1)

# text(400,0.04, paste("For Gillespie 2001, leaping started at state x>40"),col="black",font=2)
# text(440,0.03, paste("For Gillespie & Petzold 2003, leaping started at state x>20"),col="grey",font=2)







axis(2, ylim=c(0,0.1),col="black",las=1,xaxt = "n")  ## las=1 makes horizontal labels

mtext(expression(paste("leap size","(",tau,")") ),side=2,line=2.5,cex=.8)
box()


text(240,.1, "Case 1-",col="black",lwd=3,font=2)

text(430,.0995, expression(paste("Parameter values:"
                ,lambda,"=0.5",";","", mu,"=0.3",";",rho,"=0.001") ),
    col="black",lwd=3,font=2)

text(285,.094, expression(paste("Error bound:"
                ,epsilon,"=0.01") ),col="black",lwd=3,font=2)



## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(state,leapCondition2003[["case1"]],type="l",col="orange",lwd=2,ylab="",
    xlab="",axes=FALSE,ylim=c(0,2.5))
lines(state,leapCondition2001[["case1"]],type="l",col="green",lwd=2)
mtext("leap condition",side=4,col="black",line=3,cex=.8) 
axis(4, ylim=c(0,2.5), col="black",col.axis="black",las=1)


## Draw the time axis
axis(1,pretty(range(state),10),xaxt = "n")
#mtext("state x",side=1,col="black",line=2.3)  



########SECOND PLOT #######
#par(mar=c(0, 4, 4, 6) + 0.1)
#par(mar=c(0,4,2,2)
#o<-par(mar=c(0,4,0,2))

#o<-par(mar=c(0, 4, 2, 6) + 0.1)
par(mar=c(0.5, 4, 0, 6) + 0.1)
## Plot first set of data and draw its axis

plot(state,tau2001[["case2"]],type="l",col="blue",lwd=2,ylab="",xlab="",
    ylim=c(0,.1),axes=FALSE)

lines(state,tau2003[["case2"]],col="red",lwd=2,ylab="",xlab="")

abline(v=51,col="black",lty=3,lwd=1.5)


text2<- c("For HTL2001, leaping started at state x>50",
         "For HTL2003, leaping started at state x>30")
legend(x = 200, y = 0.075, legend = text2,
       col = c("black", "grey"), lty = 3, lwd = 1.5, bty = "n", cex = 1)



# text(400,0.07, paste("For Gillespie 2001, leaping started at state x>50"),col="black",font=2)
# text(440,0.06, paste("For Gillespie & Petzold 2003, leaping started at state x>30"),col="grey",font=2)


abline(v=31,col="grey",lty=3,lwd=1.5)


axis(2, ylim=c(0,0.1),col="black",las=1,xaxt = "n")  ## las=1 makes horizontal labels

mtext(expression(paste("leap size","(",tau,")") ),side=2,line=2.5,cex=.8)
box()


text(240,.1, "Case 2-",col="black",lwd=3,font=2)

text(410,.0995, expression(paste("Parameter values:"
                ,lambda,"=2",";","", mu,"=1",";",rho,"=0.01") ),
    col="black",lwd=3,font=2)

text(285,.094, expression(paste("Error bound:"
                ,epsilon,"=0.01") ),col="black",lwd=3,font=2)


text<-c("HTL2001: leap size vs state x","HTL2003: leap size vs state x",
       "HTL2001: leap condition vs state x","HTL2003: leap condition vs state x")

 legend(x=280,y=.055,legend = text,col=c("blue","red","green","orange"),bty = "n",
       cex=1.05, horiz = FALSE,pt.cex = 1,fill=c("blue","red","green","orange"),ncol=1)





## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(state,leapCondition2003[["case2"]],type="l",col="orange",lwd=2,ylab="",
    xlab="",axes=FALSE,ylim=c(0,2.5))
lines(state,leapCondition2001[["case2"]],type="l",col="green",lwd=2)
mtext("leap condition",side=4,col="black",line=3,cex=.8) 
axis(4, ylim=c(0,2.5), col="black",col.axis="black",las=1)


## Draw the time axis
axis(1,pretty(range(state),10),xaxt = "n")
#mtext("state x",side=1,col="black",line=2.3)  




########THIRD PLOT #######
#par(mar=c(4,4,0,2))
#par(mar=c(5, 4, 0, 6) + 0.1)
par(mar=c(3, 4, 0, 6) + 0.1)
## Plot first set of data and draw its axis

plot(state,tau2001[["case3"]],type="l",col="blue",lwd=2,ylab="",xlab="",
    ylim=c(0,.1),axes=FALSE)

lines(state,tau2003[["case3"]],col="red",lwd=2,ylab="",xlab="")

abline(v=31,col="black",lty=3,lwd=1.5)


text3<- c("For HTL2001, leaping started at state x>30",
         "For HTL2003, leaping started at state x>20")
legend(x = 200, y = 0.06, legend = text3,
       col = c("black", "grey"), lty = 3, lwd = 1.5, bty = "n", cex = 1)

# text(400,0.06, paste("For Gillespie 2001, leaping started at state x>30"),col="black",font=2)
# text(440,0.05, paste("For Gillespie & Petzold 2003, leaping started at state x>20"),col="grey",font=2)


abline(v=21,col="grey",lty=3,lwd=1.5)


axis(2, ylim=c(0,0.1),col="black",las=1,xaxt = "n")  ## las=1 makes horizontal labels

mtext(expression(paste("leap size","(",tau,")") ),side=2,line=2.5,cex=.8)
box()


text(240,.1, "Case 3-",col="black",lwd=3,font=2)

text(405,.0995, expression(paste("Parameter values:"
                ,lambda,"=3",";","", mu,"=2",";",rho,"=0.1") ),
    col="black",lwd=3,font=2)

text(285,.094, expression(paste("Error bound:"
                ,epsilon,"=0.01") ),col="black",lwd=3,font=2)



## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(state,leapCondition2003[["case3"]],type="l",col="orange",lwd=2,ylab="",
    xlab="",axes=FALSE,ylim=c(0,2.5))
lines(state,leapCondition2001[["case3"]],type="l",col="green",lwd=2)
mtext("leap condition",side=4,col="black",line=3, cex=.8) 
axis(4, ylim=c(0,2.5), col="black",col.axis="black",las=1)


## Draw the time axis
axis(1,pretty(range(state),10))
mtext("state x",side=1,col="black",line=2.3,cex=.8)  



par(o)
