starttable <- function(){
  return("\\begin{table}[h!]\n\\centering\n\\begin{tabular}{|c|c|c|c|c|}\n\\hline\n$T/K$&$D/V$&$C/V$&$\\omega/Hz$&$\\tau/s$\\\\\\hline\\hline")
}

makeline<-function(fitdata,T){
  D0<-fitdata["D","Estimate"]
  C0<-fitdata["C","Estimate"]
  omega0<-fitdata["omega","Estimate"]
  tau0<-fitdata["tau","Estimate"]
  sD<-fitdata["D","Std. Error"]
  sC<-fitdata["C","Std. Error"]
  somega<-fitdata["omega","Std. Error"]
  stau<-fitdata["tau","Std. Error"]
  
  D=pm(roundfunc(c(D0,sD)))
  C=pm(roundfunc(c(C0,sC)))
  omega=pm(roundfunc(c(omega0,somega)))
  tau=pm(roundfunc(c(tau0,stau)))
  
  return(paste("$",T,"$&$",D,"$&$",C,"$&$",omega,"$&$",tau,"$\\\\\\hline",sep=""))
}

endtable <- function(caption=""){
  return(paste("\\end{tabular}\n\\caption{",caption,"}\n\\end{table}\n",sep=""))
}

pm <- function(vals){
  return(paste("(",vals[1],"\\pm",vals[2],")",sep=""))
}