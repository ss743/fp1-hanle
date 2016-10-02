lorentzfit<-function(input,neg=FALSE){
  
  lorentz <- y~D+C/2*(2-(2*(x-omega)*tau)^2/(1+(2*(x-omega)*tau)^2))
  
  D0=max(input$y)
  C0=min(input$y)-max(input$y)
  omega0=input$x[which.max(input$y)]
  tau0=10^(-7)
  
  if(neg){
    C0=max(input$y)-min(input$y)
    D0=-min(input$y)+max(input$y)
    omega0=input$x[which.min(input$y)]
    
  }
  
  #print(c(D0,C0,omega0,tau0))
  
  #plot(function(x){D0+C0*((2*(x-omega0)*tau0)^2/(1+(2*(x-omega0)*tau0)^2))},-6*10^(7),6*10^(7),add=TRUE,col="green")
  
  try({
    fit=nls(lorentz,input,start=list(D=D0,C=C0,omega=omega0,tau=tau0))
    return(summary(fit)$parameters)
  })
  
  return(NULL)
}

dispersionsfit<-function(input){
  
  disp <- y~D+C*((2*(x-omega)*tau)/(1+(2*(x-omega)*tau)^2)^2)
  
  D0=(min(input$y)+max(input$y))/2
  C0=max(input$y)-min(input$y)
  omega0=(input$x[which.max(input$y)]+input$x[which.min(input$y)])/2
  tau0=10^(-7)
  

  #print(c(D0,C0,omega0,tau0))
  
  #plot(function(x){D0+C0*((2*(x-omega0)*tau0)/(1+(2*(x-omega0)*tau0)^2)^2)},-6*10^(7),6*10^(7),add=TRUE,col="green")
  
  try({
    fit=nls(disp,input,start=list(D=D0,C=C0,omega=omega0,tau=tau0))
    return(summary(fit)$parameters)
  })
  
  return(NULL)
}


plotlorentz<-function(fitdata,bereich){
  D<-fitdata["D","Estimate"]
  C<-fitdata["C","Estimate"]
  omega<-fitdata["omega","Estimate"]
  tau<-fitdata["tau","Estimate"]
  
  try({plot(function(x){D+C/2*(2-(2*(x-omega)*tau)^2/(1+(2*(x-omega)*tau)^2))},-6*10^(7),6*10^(7),add=TRUE,col="red")})
}

plotdisp<-function(fitdata,bereich){
  D<-fitdata["D","Estimate"]
  C<-fitdata["C","Estimate"]
  omega<-fitdata["omega","Estimate"]
  tau<-fitdata["tau","Estimate"]
  
  try({plot(function(x){D+C*((2*(x-omega)*tau)/(1+(2*(x-omega)*tau)^2)^2)},-6*10^(7),6*10^(7),add=TRUE,col="red")})
}

getlorentzvalue<-function(fitdata,x){
  D<-fitdata["D","Estimate"]
  C<-fitdata["C","Estimate"]
  omega<-fitdata["omega","Estimate"]
  tau<-fitdata["tau","Estimate"]

  return(D+C/2*(2-(2*(x-omega)*tau)^2/(1+(2*(x-omega)*tau)^2)))    
}

getdispvalue<-function(fitdata,x){
  D<-fitdata["D","Estimate"]
  C<-fitdata["C","Estimate"]
  omega<-fitdata["omega","Estimate"]
  tau<-fitdata["tau","Estimate"]
  
  return(D+C*((2*(x-omega)*tau)/(1+(2*(x-omega)*tau)^2)^2))    
}

getlorentzformula<-function(fitdata){
  D0<-fitdata["D","Estimate"]
  C0<-fitdata["C","Estimate"]
  omega0<-fitdata["omega","Estimate"]
  tau0<-fitdata["tau","Estimate"]
  sD<-fitdata["D","Std. Error"]
  sC<-fitdata["C","Std. Error"]
  somega<-fitdata["omega","Std. Error"]
  stau<-fitdata["tau","Std. Error"]
  
  D=roundfunc(D0,sD)
  C=roundfunc(C0,sC)
  omega=roundfunc(omega0,somega)
  tau=roundfunc(tau0,stau)
  
  txt=paste("(",D[1],"+-",D[2],") V + (",C[1]/2,"+-",C[2]/2,") V * (2-(2(omega_L-(",omega[1],"+-",omega[2],")Hz)*(",tau[1],"+-",tau[2],")s)^2/(1+(2(omega_L-(",omega[1],"+-",omega[2],")Hz)*(",tau[1],"+-",tau[2],")s)^2))",sep="")
  
  return(txt)
}

getdispformula<-function(fitdata){
  D0<-fitdata["D","Estimate"]
  C0<-fitdata["C","Estimate"]
  omega0<-fitdata["omega","Estimate"]
  tau0<-fitdata["tau","Estimate"]
  sD<-fitdata["D","Std. Error"]
  sC<-fitdata["C","Std. Error"]
  somega<-fitdata["omega","Std. Error"]
  stau<-fitdata["tau","Std. Error"]
  
  D=roundfunc(c(D0,sD))
  C=roundfunc(c(C0,sC))
  omega=roundfunc(c(omega0,somega))
  tau=roundfunc(c(tau0,stau))
  
  txt=paste("(",D[1],"+-",D[2],") V + (",C[1],"+-",C[2],") V * (2-(2(omega_L-(",omega[1],"+-",omega[2],")Hz)*(",tau[1],"+-",tau[2],")s)/(1+(2(omega_L-(",omega[1],"+-",omega[2],")Hz)*(",tau[1],"+-",tau[2],")s)^2)^2)",sep="")
  
  return(txt)
}