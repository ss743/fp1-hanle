bildsize=c(384,256)


source("aufwaermen.R")
tau1=tau_end
stau1=stau_end
cat("\n")
source("abkuehlen.R")
tau2=tau_end
stau2=stau_end


tau=weighted.mean(c(tau1,tau2),1/c(stau1,stau2)^2)
stau=sqrt(1/sum(1/c(stau1,stau2)^2))

results=roundfunc(c(tau,stau))

cat(paste("\nTau = (",results[1]," +- ",results[2],") ns\n",sep=""))

tau45=weighted.mean(c(tau1_45,tau2_45),1/c(stau1_45,stau2_45)^2)
stau45=sqrt(1/sum(1/c(stau1_45,stau2_45)^2))

results=roundfunc(c(tau45,stau45))

cat(paste("\nTau45 = (",results[1]," +- ",results[2],") ns\n",sep=""))