library(rtdists)
library(tidyverse)
library(gtools)
library(lattice)
library(fitdistrplus)
library(logspline)
library(ggplot2)
library(gridExtra)
library(gridGraphics)
library(grid)
library(RColorBrewer)

T = 40 #number of decision timesteps (weeks)
# L's refer to value at which first lockdown is switched to (number of deaths)
# E'S refer to value of easing of lockdown is switched to (proportion of cases)
L1=100
L2=300
L3=500

E0=0
E1=0.12
E2=0.3
E3=0.5

sits=array(c(L1,E0,0,
             L1,E1,0,
             L1,E2,0,
             L1,E3,0,
             L2,E0,0,
             L2,E1,0,
             L2,E2,0,
             L2,E3,0,
             L3,E0,0,
             L3,E1,0,
             L3,E2,0,
             L3,E3,0,
             L1,0,1,
             L2,0,1,
             L3,0,1),
           dim=c(3,15))
strat_names=c("L1_E0_0",
              "L1_E1_0",
              "L1_E2_0",
              "L1_E3_0",
              "L2_E0_0",
              "L2_E1_0",
              "L2_E2_0",
              "L2_E3_0",
              "L3_E0_0",
              "L3_E1_0",
              "L3_E2_0",
              "L3_E3_0",
              "L1_E0_1",
              "L2_E0_1",
              "L3_E0_1")

setwd("~/Documents/covidDSS")
pop=read.csv("population.csv") #input population sizes
lyl=read.csv("lifeexpectancy.csv")  #input life expectancy sizes
rates<-read.csv("death+recovery.csv") #input transmission rates between age groups
intialcases<-read.csv("cases.csv") #input inisiation of cases
delay<-read.csv("delay.csv")  #input lyl due to cancer 
DEP<-read.csv("dep.csv") #input lyl due to deprevation 

death_array=c(lyl$Males,lyl$Females)
cases=intialcases$X06.Mar

t = 7 #how many discrete timesteps between each decisions days
Runs=T*t-1
agecon=rep(rates$contacts,3)
con=cbind(rep(3,14),cbind(0.5*agecon,agecon))
beta=0.023
agedeath=rates$death
death=rep(agedeath,2)
agerecovery=rates$recovery
recovery=rep(agerecovery,2)

## 3 is no lockdown
## 2 is intermediate lockdown
## 1 is complete lockdown
lag=7 #how many days behind are the tests from the truth 
D_lag=0
L_val=1.05


pops=as.matrix(pop[,3:16])
S=array(rep(0,12*14*t*T),dim=c(12,14,t*T)) #3d array location/classification/time
I=array(rep(0,12*14*t*T),dim=c(12,14,t*T)) #3d array location/classification/time
R=array(rep(0,12*14*t*T),dim=c(12,14,t*T)) #3d array location/classification/time
D=array(rep(0,12*14*t*T),dim=c(12,14,t*T)) #3d array location/classification/time
I[,,1]=pops/rowSums(pops)*cases/rowSums(pops)
S[,,1]=pops/rowSums(pops)-I[,,1]

coviddeaths=array(0,dim=c(14,15))
strat=array(rep(3,T*15),dim=c(T,15)) #which strategy is active
par(mfrow=c(3,4))
for(j in 1:15){
  L_time=0
  
  for (i in 1:Runs){
    if(i>D_lag){
      if((sum(D[,,(i-D_lag)]*rowSums(pops))>sits[1,j])&&(L_time==0)){
        L_time=i
        if(sits[3,j]==0){
          strat[i%/%7+1,j]=1}
        else{
          strat[i%/%7+1,j]=2}
      }
    }
    
    if((i%%7==0)&&(L_time!=0)){
      if((i>L_time)&&(sits[3,j]==0)){
        a=0
        if (strat[i/7,j]!=1){
          if(sum(I[,,(i-lag)]*rowSums(pops))>L_val*sum(I[,,(i-7-lag)]*rowSums(pops))){
            L_time<<-i
            strat[i/7+1,j]=1
            a=1
          }}
        if (strat[i/7,j]==1){
          if(max(colSums(colSums(aperm(I[,,(L_time-lag):(i-lag)]*rowSums(pops), c(2,1,3)))))*sits[2,j]>sum(I[,,(i-lag)]*rowSums(pops))){
            strat[i/7+1,j]=2
            a=1
          }}
        if(a==0){
          strat[i/7+1,j]=strat[i/7,j]
        }}
      
      else{
        strat[i/7+1,j]=strat[i/7,j]
      }
    }
    S[,,i+1]=S[,,i]-t(con[,strat[i%/%7+1,j]]*beta*t(S[,,i]))*rowSums(I[,,i])
    I[,,i+1]=I[,,i]+t(con[,strat[i%/%7+1,j]]*beta*t(S[,,i]))*rowSums(I[,,i])-t((death+recovery)*t(I[,,i]))
    R[,,i+1]=R[,,i]+t(recovery*t(I[,,i]))
    D[,,i+1]=D[,,i]+t(death*t(I[,,i]))
  }
  S_all=colSums(colSums(aperm(S, c(2,1,3)))*rowSums(pops))
  I_all=colSums(colSums(aperm(I, c(2,1,3)))*rowSums(pops))
  R_all=colSums(colSums(aperm(R, c(2,1,3)))*rowSums(pops))
  D_all=colSums(colSums(aperm(D, c(2,1,3)))*rowSums(pops))
  
  #plot(S_all,col='green',ylim=c(0,R_all[Runs]+D_all[Runs]),xlab = "Day",ylab="Number of people",main=strat_names[j])
  #points(I_all,col='red')
  #points(R_all,col='blue')
  #points(D_all,col='black')
  
  #plot(I_all,col='red',ylim=c(0,1.2*max(max(I_all),max(D_all))),xlab = "Day",ylab="Number of people",main=strat_names[j])
  #points(D_all,col='black')
  #legend("topright" , legend=c("Infected","Dead"), col=c("red","black"),lty=1:1,cex=1)
  
  
  #plot(I_all,col='red',xlab = "Day",ylab="Number of people",main=strat_names[j])
  #points(D_all,col='black')
  #legend("topright" , legend=c("Infected","Dead"), col=c("red","black"),lty=1:1,cex=1)
  
  #legend(150, 50000000, legend=c("Suceptible","Infected","Recovered","Dead"), col=c("green", "red","blue","black"),lty=1:1,cex=1)
  coviddeaths[,j]=colSums(D[,,Runs]*pops)
  
}
strat
delaydeaths=array(0,dim=c(14,15))
stratcount=array(0,dim=c(3,15))
DEPLYL=array(0,dim=c(14,15))

for (j in 1:15){
  stratcount[1,j]=sum(strat[,j]==1)
  stratcount[2,j]=sum(strat[,j]==2)
  stratcount[3,j]=sum(strat[,j]==3)
  delaydeaths[,j]=as.numeric(delay[1,]*stratcount[1,j]+delay[2,]*stratcount[2,j]+delay[3,]*stratcount[3,j])
  DEPLYL[,j]=as.numeric(DEP[1,]*stratcount[1,j]+DEP[2,]*stratcount[2,j]+DEP[3,]*stratcount[3,j])
}
delay
coviddeaths
delaydeaths


par(mfrow=c(1,1))

grouped_coviddeaths=coviddeaths[1:7,]+coviddeaths[8:14,]
grouped_delaydeaths=delaydeaths[1:7,]+delaydeaths[8:14,]
groupedcovidLYL=coviddeaths[1:7,]*death_array[1:7]+coviddeaths[8:14,]*death_array[8:14]
groupeddelayLYL=delaydeaths[1:7,]*death_array[1:7]+delaydeaths[8:14,]*death_array[8:14]
groupedDEPLYL=DEPLYL[1:7,]+DEPLYL[8:14,]


labels=colnames(pops)[1:7]
labels=sub('.', '', labels)
rownames(grouped_coviddeaths)=labels

#colnames(grouped_coviddeaths)=strat_names
#rownames(grouped_coviddeaths)=labels

#colnames(grouped_delaydeaths)=strat_names
#rownames(grouped_delaydeaths)=labels

colnames(groupedcovidLYL)=strat_names
rownames(groupedcovidLYL)=labels

colnames(groupeddelayLYL)=strat_names
rownames(groupeddelayLYL)=labels

colnames(groupedDEPLYL)=strat_names
rownames(groupedDEPLYL)=labels


y=grouped_coviddeaths
mybarplot<-function(y,main_name,ylab_name){
  barplot(y,
          xlim=c(0, ncol(y) + 10),
          col=brewer.pal(nrow(y), "Paired"),
          ylab=ylab_name,
          legend.text=TRUE,
          args.legend=list(
            x=ncol(y) + 10,
            y=max(colSums(y)),
            bty = "n"
          ),las=3,
          main=main_name
  )
}
par(mfrow=c(1,3))

#mybarplot(grouped_coviddeaths,"Deaths due to Covid-19 under different strategies","Total deaths")
#mybarplot(grouped_delaydeaths,"Deaths due to delays in cancer treatments under different strategies","Total deaths")
#mybarplot(grouped_coviddeaths+grouped_delaydeaths,"Total deaths under different strategies","Total deaths")

mybarplot(groupedcovidLYL,"Life years lost due to covid under different strategies","Total life years lost")
mybarplot(groupeddelayLYL,"life years lost due to delays in cancer treatments under different strategies","Total life years lost")
mybarplot(groupedDEPLYL,"life years lost due to deprivation under different strategies","Total life years lost")
mybarplot(groupedcovidLYL+groupeddelayLYL+groupedDEPLYL,"Total life years lost under different strategies","Total life years lost")

par(mfrow=c(1,1))
plot(colSums(delaydeaths),colSums(coviddeaths),col="blue",xlab="Cancer deaths",ylab = "Covid deaths", main="Comparision of different strategies on causes of death")
pos_vector <- rep(3, length(strat_names))
pos_vector[strat_names %in% c("L2_E1_0","L1_E2_0","L3_E3_0","L1_E0_1")] <- 4
pos_vector[strat_names %in% c("L3_E2_0","L2_E3_0","L2_E0_1")] <- 1
pos_vector[strat_names %in% c("L2_E0_0")] <- 2

text(colSums(delaydeaths),colSums(coviddeaths), labels=strat_names,cex=0.8, font=2,pos=pos_vector)

#plot(colSums(delaydeaths)+colSums(coviddeaths))

