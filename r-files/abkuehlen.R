source("readFiles.R")
source("linearfit.R")
source("lorentzfit.R")
source("round.R")
source("druck.R")

library("plotrix")

BIfactor=3.363*10^(-4)
omegaBfactor=1.4838*9.274009994*10^(-24)/(1.054571800*10^(-34))

N0 = 18
N45= 18
N90= 18

offset=273.15
T0=c(18,17,16,15,14,13,12,11,10,9,8,7,6,4,3,2,1,-1)+offset
T45=c(18,17,16,15,14,13,12,10,9,8,7,6,5,4,3,2,1,-1)+offset
T90=c(18,16,15,14,13,12,11,10,9,8,7,6,4,3,2,1,0,-3)+offset

P0=druck(T0)
P45=druck(T45)
P90=druck(T90)

tau0=c()
tau45=c()
tau90=c()

intense0=c()
intense45=c()
intense90=c()

stau0=c()
stau45=c()
stau90=c()

messung0grad=array(dim=c(2000,3,N0))
messung45grad=array(dim=c(2000,3,N45))
messung90grad=array(dim=c(2000,3,N90))

prefix="abkuehlung/"
prefix2="messung"

for(i in 1:54){
  n=i
  if(i<10)
    n=paste("0",i,sep="")
  mess=readCSV(paste(prefix,prefix2,n,sep=""))
  messung=array(dim=c(2000,3))
  messung[,1]=mess$x
  messung[,2]=mess$y1
  messung[,3]=mess$y2
  
  if(i%%3==1)
    messung0grad[,,i%/%3+1]=messung
  if(i%%3==2)
    messung45grad[,,i%/%3+1]=messung
  if(i%%3==0)
    messung90grad[,,i%/%3]=messung
  
}
plotprefix="../figures/cold"

table0=starttable()
for(i in 1:N0){
  t=messung0grad[,1,i]
  U=messung0grad[,2,i]
  Int=messung0grad[,3,i]
  #plot(t,messung0grad[,2,i],cex=0.4,pch=4,bty="l",xlab="t / s", ylab="U / V")
  fit=linearfit(data.frame(x=t,y=U),bereich=c(20,30))
  b=fit[1]
  a=fit[2]
  sb=fit[3]
  sa=fit[4]
  I=a*t+b
  sI=sqrt((sa*t)^2+sb^2)
  B=BIfactor*I
  sB=BIfactor*sI
  omegaL=omegaBfactor*B
  somegaL=omegaBfactor*sB
  #plotlinear(fit,c(0,50))
  #grid()
  png(paste(plotprefix,"0-",i,".png",sep=""))
  plot(omegaL,Int,xlab="omega_L / Hz", ylab="Intensität / V",cex=0.6,pch=4,bty="l")
  grid()
  fit=lorentzfit(data.frame(x=omegaL,y=Int))
  plotlorentz(fit,c(-6,6)*10^7)
  title(paste("Abkühlvorgang - Nr. ",i," - T=",T0[i],"K; Pol=0°",sep=""))
  dev.off()
  tau0[i]=abs(fit["tau","Estimate"])
  stau0[i]=fit["tau","Std. Error"]
  intense0[i]=abs(getlorentzvalue(fit,fit["omega","Estimate"])-fit["D","Estimate"])
  table0=paste(table0,makeline(fit,T0[i]),sep="\n")
}
table0=paste(table0,endtable("Abkühlvorgang, Polarisation: 0°",label="cold0"),sep="\n")
cat(table0,file="../tables/cold0.tex")

table90=starttable()
for(i in 1:N90){
  #print(i)
  t=messung90grad[,1,i]
  U=messung90grad[,2,i]
  Int=messung90grad[,3,i]
  #plot(t,U,cex=0.4,pch=4,bty="l",xlab="t / s", ylab="U / V")
  fit=linearfit(data.frame(x=t,y=U),bereich=c(20,30))
  b=fit[1]
  a=fit[2]
  sb=fit[3]
  sa=fit[4]
  I=a*t+b
  sI=sqrt((sa*t)^2+sb^2)
  B=BIfactor*I
  sB=BIfactor*sI
  omegaL=omegaBfactor*B
  somegaL=omegaBfactor*sB
  #plotlinear(fit,c(0,50))
  #grid()
  png(paste(plotprefix,"90-",i,".png",sep=""))
  plot(omegaL,Int,xlab="omega_L / Hz", ylab="Intensität / V",cex=0.6,pch=4,bty="l")
  grid()
  fit=lorentzfit(data.frame(x=omegaL,y=Int),neg=TRUE)
  plotlorentz(fit,c(-6,6)*10^7)
  title(paste("Abkühlvorgang - Nr. ",i," - T=",T90[i],"K; Pol=90°",sep=""))
  dev.off()
  
  tau90[i]=abs(fit["tau","Estimate"])
  stau90[i]=fit["tau","Std. Error"]
  intense90[i]=abs(getlorentzvalue(fit,fit["omega","Estimate"])-fit["D","Estimate"])
  table90=paste(table90,makeline(fit,T0[i]),sep="\n")
}
table90=paste(table90,endtable("Abkühlvorgang, Polarisation: 90°",label="cold90"),sep="\n")
cat(table90,file="../tables/cold90.tex")

table45=starttable()
for(i in 1:N45){
  t=messung45grad[,1,i]
  U=messung45grad[,2,i]
  Int=messung45grad[,3,i]
  #plot(t,U,cex=0.4,pch=4,bty="l",xlab="t / s", ylab="U / V")
  fit=linearfit(data.frame(x=t,y=U),bereich=c(20,30))
  b=fit[1]
  a=fit[2]
  sb=fit[3]
  sa=fit[4]
  I=a*t+b
  sI=sqrt((sa*t)^2+sb^2)
  B=BIfactor*I
  sB=BIfactor*sI
  omegaL=omegaBfactor*B
  somegaL=omegaBfactor*sB
  #plotlinear(fit,c(0,50))
  #grid()
  png(paste(plotprefix,"45-",i,".png",sep=""))
  plot(omegaL,Int,xlab="omega_L / Hz", ylab="Intensität / V",cex=0.6,pch=4,bty="l")
  grid()
  fit=dispersionsfit(data.frame(x=omegaL,y=Int))
  plotdisp(fit,c(-6,6)*10^7)
  title(paste("Abkühlvorgang - Nr. ",i," - T=",T45[i],"K; Pol=45°",sep=""))
  dev.off()
  
  tau45[i]=abs(fit["tau","Estimate"])
  stau45[i]=fit["tau","Std. Error"]
  intense45[i]=1/2*abs(optimize(function(x){getdispvalue(fit,x)},c(-2,2)*10^7,maximum=TRUE)$objective-optimize(function(x){getdispvalue(fit,x)},c(-2,2)*10^7,maximum=FALSE)$objective)
  table45=paste(table45,makeline(fit,T0[i]),sep="\n")
}
table45=paste(table45,endtable("Abkühlvorgang, Polarisation: 45°",label="cold45"),sep="\n")
cat(table45,file="../tables/cold45.tex")


fit0=linearfit(data.frame(x=P0,y=tau0,sy=stau0),weighted=TRUE)
fit45=linearfit(data.frame(x=P45,y=tau45,sy=stau45),weighted=TRUE)
fit90=linearfit(data.frame(x=P90,y=tau90,sy=stau90),weighted=TRUE)

plotCI(P0,tau0,uiw=stau0,cex=0.6,pch=4,bty="l",xlab="p / Pa", ylab="Tau / s")
plotlinear(fit0,c(0,300))
grid()
title("Polarisation 0°")
plotCI(P90,tau90,uiw=stau90,cex=0.6,pch=4,bty="l",xlab="p / Pa", ylab="Tau / s")
plotlinear(fit90,c(0,300))
grid()
title("Polarisation 90°")
plotCI(P45,tau45,uiw=stau45,cex=0.6,pch=4,bty="l",xlab="p / Pa", ylab="Tau / s")
plotlinear(fit45,c(0,300))
grid()
title("Polarisation 45°")

fitInt0=linearfit(data.frame(x=T0,y=intense0),weighted=FALSE)
fitInt45=linearfit(data.frame(x=T45,y=intense45),weighted=FALSE)
fitInt90=linearfit(data.frame(x=T90,y=intense90),weighted=FALSE)

plot(T0,intense0,cex=0.6,pch=4,bty="l",xlab="T / K", ylab="Amplitude / V")
plotlinear(fitInt0,c(0,300))
grid()
title("Polarisation 0°")
plot(T90,intense90,cex=0.6,pch=4,bty="l",xlab="T / K", ylab="Amplitude / V")
plotlinear(fitInt90,c(0,300))
grid()
title("Polarisation 90°")
plot(T45,intense45,cex=0.6,pch=4,bty="l",xlab="T / K", ylab="Amplitude / V")
plotlinear(fitInt45,c(0,300))
grid()
title("Polarisation 45°")


tau_ext=c()
stau_ext=c()

tau_ext[1]=fit0[1]
stau_ext[1]=fit0[3]
tau_ext45=fit45[1]
stau_ext45=fit45[3]
tau_ext[2]=fit90[1]
stau_ext[2]=fit90[3]

tau_end=weighted.mean(tau_ext,1/stau_ext^2)*10^9
stau_end=sqrt(1/sum(1/stau_ext^2))*10^9

results=roundfunc(c(tau_end,stau_end))

cat(paste("Tau2 = (",results[1]," +- ",results[2],") ns\n",sep=""))