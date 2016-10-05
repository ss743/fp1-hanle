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

chi0=c()
chi45=c()
chi90=c()


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
  png(paste(plotprefix,"0-",i,".png",sep=""),width=768,height=512)
  plot(omegaL,Int,xlab="omega_L / Hz", ylab="Intensität / V",cex=0.6,pch=4,bty="l")
  grid()
  fit=lorentzfit(data.frame(x=omegaL,y=Int))
  plotlorentz(fit,c(-6,6)*10^7)
  title(paste("Abkühlvorgang - Nr. ",i," - T=",T0[i],"K; Pol=0°",sep=""))
  dev.off()
  chi0[i]=fit[5,1]
  tau0[i]=abs(fit["tau","Estimate"])
  stau0[i]=fit["tau","Std. Error"]
  intense0[i]=fit["D","Estimate"]
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
  png(paste(plotprefix,"90-",i,".png",sep=""),width=768,height=512)
  plot(omegaL,Int,xlab="omega_L / Hz", ylab="Intensität / V",cex=0.6,pch=4,bty="l")
  grid()
  fit=lorentzfit(data.frame(x=omegaL,y=Int),neg=TRUE)
  plotlorentz(fit,c(-6,6)*10^7)
  title(paste("Abkühlvorgang - Nr. ",i," - T=",T90[i],"K; Pol=90°",sep=""))
  dev.off()
  chi90[i]=fit[5,1]
  
  tau90[i]=abs(fit["tau","Estimate"])
  stau90[i]=fit["tau","Std. Error"]
  intense90[i]=fit["D","Estimate"]
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
  png(paste(plotprefix,"45-",i,".png",sep=""),width=768,height=512)
  plot(omegaL,Int,xlab="omega_L / Hz", ylab="Intensität / V",cex=0.6,pch=4,bty="l")
  grid()
  fit=dispersionsfit(data.frame(x=omegaL,y=Int))
  plotdisp(fit,c(-6,6)*10^7)
  title(paste("Abkühlvorgang - Nr. ",i," - T=",T45[i],"K; Pol=45°",sep=""))
  dev.off()
  
  chi45[i]=fit[5,1]
  tau45[i]=abs(fit["tau","Estimate"])
  stau45[i]=fit["tau","Std. Error"]
  intense45[i]=fit["D","Estimate"]
  table45=paste(table45,makeline(fit,T0[i]),sep="\n")
}
table45=paste(table45,endtable("Abkühlvorgang, Polarisation: 45°",label="cold45"),sep="\n")
cat(table45,file="../tables/cold45.tex")


fit0=linearfit(data.frame(x=P0,y=tau0,sy=stau0),weighted=TRUE)
fit45=linearfit(data.frame(x=P45,y=tau45,sy=stau45),weighted=TRUE)
fit90=linearfit(data.frame(x=P90,y=tau90,sy=stau90),weighted=TRUE)

cat(paste("\nK_tau, 0°: Chi_quadrat=",fit0[5]*10^9,sep=""))
cat(paste("\nK_tau, 45°: Chi_quadrat=",fit45[5]*10^9,sep=""))
cat(paste("\nK_tau, 90°: Chi_quadrat=",fit90[5]*10^9,sep=""))

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

plotlindata(fit0,"K, 0°")
plotlindata(fit45,"K, 45°")
plotlindata(fit90,"K, 90°")

si0 =abs(intense0 *0.03)
si45=abs(intense45*0.03)
si90=abs(intense90*0.03)

sT0 =0.5
sT45=0.5
sT90=0.5

fitInt0=linearfit(data.frame(x=T0,y=intense0,sy=si0),weighted=TRUE)
fitInt45=linearfit(data.frame(x=T45,y=intense45,sy=si45),weighted=TRUE)
fitInt90=linearfit(data.frame(x=T90,y=intense90,sy=si90),c(271,300),weighted=TRUE)

cat(paste("\nK_i, 0°: Chi_quadrat=",fitInt0[5],sep=""))
cat(paste("\nK_i, 45°: Chi_quadrat=",fitInt45[5],sep=""))
cat(paste("\nK_i, 90°: Chi_quadrat=",fitInt90[5],sep=""))

plotCI(T0,intense0,uiw=si0,cex=0.6,pch=4,bty="l",xlab="T / K", ylab="Grundintensität / V")
plotCI(T0,intense0,uiw=sT0,err="x",cex=0.6,pch=4,add=TRUE)
plotlinear(fitInt0,c(0,300))
grid()
title("Polarisation 0°")
plotCI(T90,intense90,uiw=si90,cex=0.6,pch=4,bty="l",xlab="T / K", ylab="Grundintensität / V")
plotCI(T90,intense90,uiw=sT90,err="x",cex=0.6,pch=4,add=TRUE)
plotlinear(fitInt90,c(0,300))
grid()
title("Polarisation 90°")
plotCI(T45,intense45,uiw=si45,cex=0.6,pch=4,bty="l",xlab="T / K", ylab="Grundintensität / V")
plotCI(T45,intense45,uiw=sT45,err="x",cex=0.6,pch=4,add=TRUE)
plotlinear(fitInt45,c(0,300))
grid()
title("Polarisation 45°")

cat("\n")
plotlindata(fitInt0,"Kalt - Polarisation 0°")
plotlindata(fitInt45,"Kalt - Polarisation 45°")
plotlindata(fitInt90,"Kalt - Polarisation 90°")


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
results45=roundfunc(c(tau_ext45,stau_ext45))*10^9
results90=roundfunc(c(tau_ext[2],stau_ext[2]))*10^9
results0=roundfunc(c(tau_ext[1],stau_ext[1]))*10^9

cat(paste("Tau2_0  = (",results0[1]," +- ",results0[2],") ns\n",sep=""))
cat(paste("Tau2_90 = (",results90[1]," +- ",results90[2],") ns\n",sep=""))
cat(paste("Tau2 = (",results[1]," +- ",results[2],") ns\n",sep=""))
cat(paste("Tau2_45 = (",results45[1]," +- ",results45[2],") ns\n",sep=""))

tau2_45=results45[1]
stau2_45=results45[2]