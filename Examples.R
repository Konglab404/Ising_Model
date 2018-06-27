#This Script present some examples for ising model;
#20180627 
#Wanghaotian
library(magrittr)
#Construct A 1000 x 1000 "Magn Matrix" as the example matrix
#MAg<-matrix(sample(x = c(1,-1),replace = T,size = 100*100),
#            nrow = 100,ncol = 100)

#Define a function to change the matrix status randomly
ChangeS<-function(m,loc){
  ma<-m
  ma[loc]<-ma[loc] * -1
#  print(x)
#  print(y)
  return(ma)
}

#Example-A Fix model
##Define a function to calculate the ensemble energy for fixed couping constrant J and magFeild H
##H is 1 or -1
CalcuE<-function(m,J,H){
  Sm<-sum(m[1:nrow(m)-1,]*m[2:nrow(m),]) + #The sum of m<i,j>
    sum(m[,1:ncol(m)-1]*m[,2:ncol(m)])
  E<- -J *Sm - H *sum(m)
  return(E)
}

## Define a function to perform Monte Carlo evolution
## mNew is the new matrix, mOld is the old matrix
## T is the temperature, k is the fixed Boltzmann Constant 1.3806488 x 10^(23)
## The return is a bool result, if output is True,  matrix can be evolved, nor the matirx wont be evolved
accp<-function(mNew,mOld,k = 1,T = 2){
  mu<-min(exp((CalcuE(mOld,J = 0.5,H=0) - CalcuE(mNew,J=0.5,H=0))/(k*T)),1)
  sz<-runif(1,0,1)
  return(sz < mu)
}

MC_step<-function(m){
  dt<-ncol(m)*nrow(m)
  neworder<-1:dt
  mTemp<-m
  for(i in 1:dt){
    mNew<-ChangeS(mTemp,neworder[i])
    p<-accp(mNew,mTemp,T = 2)
    if(p == T){mTemp<-mNew}
  }
  return(mTemp)
}

##The simulation of above expamle
MAg<-matrix(sample(c(-1,1),size = 2500,replace = T),
            nrow = 50,ncol = 50)
orin<-MAg
pheatmap(orin,cluster_rows = F,cluster_cols = F)
#Energy<-data.frame(time = 1:3000,energy = 0,mean_spin = 0,accp=T,
#                   x = 0,y=0)
for(i in 1:1000){
  MAg<-MC_step(MAg)
}
pheatmap(MAg,cluster_rows = F,cluster_cols = F)
#plot(Energy$time,Energy$energy,type = "l")






