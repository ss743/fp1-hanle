druck <- function(T){
  Tc=1764
  pc=167000000
  a1=-4.57618368
  a2=-1.40726277
  a3=2.36263541
  a4=-31.0889985
  a5=58.0183959
  a6=-27.6304546
  
  Tr=(1-T/Tc)
  
  return(pc*exp(Tc/T*(a1*Tr+a2*Tr^(1.89)+a3*Tr^2+a4*Tr^8+a5*Tr^(8.5)+a6*Tr^9)))
}

sdruck <- function(p,T,sT){
  sig = matrix(c( 0.002223462, -0.039761786,  0.038596482, -0.057150539,  0.104904983, -0.048841224,
                      -0.039761786,  0.713627443, -0.693021642,  1.050624836, -1.931983155,  0.901047084,
                      0.038596482, -0.693021642,  0.673047757, -1.023299883,  1.882148471, -0.877990861,
                      -0.057150539,  1.050624836, -1.023299883,  1.805980893, -3.35831368,   1.583127006,
                      0.104904983, -1.931983155,  1.882148471, -3.35831368,   6.249787823, -2.948350916,
                      -0.048841224,  0.901047084, -0.877990861,  1.583127006, -2.948350916,  1.391862893),
                   nrow=6, ncol=6, 
                   byrow=TRUE
    )
  
  Tc=1764
  pc=167000000
  a=c()
  a[1]=-4.57618368
  a[2]=-1.40726277
  a[3]=2.36263541
  a[4]=-31.0889985
  a[5]=58.0183959
  a[6]=-27.6304546
  
  x=c(1,1.89,2,8,8.5,9)
  
  Tr=(1-T/Tc)
  
  dpdT=-1/Tc*p/(1-Tr)*sum((a*Tr^x)*(-1/(1-Tr)+x/Tr))
  dpda=p/(1-Tr)*Tr^x
  dpda1=p/(1-Tr[1])*Tr[1]^x[1]#dpda[1]*dpda*sig[1,]
  dpda2=p/(1-Tr[2])*Tr[2]^x[2]#dpda[2]*dpda*sig[2,]
  dpda3=p/(1-Tr[3])*Tr[3]^x[3]#dpda[3]*dpda*sig[3,]
  dpda4=p/(1-Tr[4])*Tr[4]^x[4]#dpda[4]*dpda*sig[4,]
  dpda5=p/(1-Tr[5])*Tr[5]^x[5]#dpda[5]*dpda*sig[5,]
  dpda6=p/(1-Tr[6])*Tr[6]^x[6]#dpda[6]*dpda*sig[6,]
  
  dpda11=dpda1*dpda1*sig[1,1]
  dpda12=dpda1*dpda2*sig[1,2]
  dpda13=dpda1*dpda3*sig[1,3]
  dpda14=dpda1*dpda4*sig[1,4]
  dpda15=dpda1*dpda5*sig[1,5]
  dpda16=dpda1*dpda6*sig[1,6]

  dpda21=dpda2*dpda1*sig[2,1]
  dpda22=dpda2*dpda2*sig[2,2]
  dpda23=dpda2*dpda3*sig[2,3]
  dpda24=dpda2*dpda4*sig[2,4]
  dpda25=dpda2*dpda5*sig[2,5]
  dpda26=dpda2*dpda6*sig[2,6]

  dpda31=dpda3*dpda1*sig[3,1]
  dpda32=dpda3*dpda2*sig[3,2]
  dpda33=dpda3*dpda3*sig[3,3]
  dpda34=dpda3*dpda4*sig[3,4]
  dpda35=dpda3*dpda5*sig[3,5]
  dpda36=dpda3*dpda6*sig[3,6]

  dpda41=dpda4*dpda1*sig[4,1]
  dpda42=dpda4*dpda2*sig[4,2]
  dpda43=dpda4*dpda3*sig[4,3]
  dpda44=dpda4*dpda4*sig[4,4]
  dpda45=dpda4*dpda5*sig[4,5]
  dpda46=dpda4*dpda6*sig[4,6]

  dpda51=dpda5*dpda1*sig[5,1]
  dpda52=dpda5*dpda2*sig[5,2]
  dpda53=dpda5*dpda3*sig[5,3]
  dpda54=dpda5*dpda4*sig[5,4]
  dpda55=dpda5*dpda5*sig[5,5]
  dpda56=dpda5*dpda6*sig[5,6]

  dpda61=dpda6*dpda1*sig[6,1]
  dpda62=dpda6*dpda2*sig[6,2]
  dpda63=dpda6*dpda3*sig[6,3]
  dpda64=dpda6*dpda4*sig[6,4]
  dpda65=dpda6*dpda5*sig[6,5]
  dpda66=dpda6*dpda6*sig[6,6]

  dpda1sum=dpda11+dpda12+dpda13+dpda14+dpda15+dpda16
  dpda2sum=dpda21+dpda22+dpda23+dpda24+dpda25+dpda26
  dpda3sum=dpda31+dpda32+dpda33+dpda34+dpda35+dpda36
  dpda4sum=dpda41+dpda42+dpda43+dpda44+dpda45+dpda46
  dpda5sum=dpda51+dpda52+dpda53+dpda54+dpda55+dpda56
  dpda6sum=dpda61+dpda62+dpda63+dpda64+dpda65+dpda66
  
  sp=sqrt(dpdT^2*sT^2+dpda1sum+dpda2sum+dpda3sum+dpda4sum+dpda5sum+dpda6sum)
  
  return(sp)
}
