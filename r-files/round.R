roundfunc <- function(vals){
    x=vals[1]
    xerr=vals[2]
    n=0
    for(i in -20:20){
      a=round(xerr,i)*10^i
      if(a==1){
        n=i+1
        return(c(round(x,n),round(xerr,n)))
      }
      if(a==2){
        if(xerr*10^i<1.95){
          n=i+1
        } else {
          n=i
        }
        return(c(round(x,n),round(xerr,n)))
      }
      if(a>2){
        n=i
        return(c(round(x,n),round(xerr,n)))
      }
  }
return(vals)

}