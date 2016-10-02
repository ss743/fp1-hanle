readCSV <- function(name){
  
  daten=read.csv(paste("../data/",name,"_HM1508.csv",sep=""))
  
  return(data.frame(x=daten[[1]],y1=daten[[2]],y2=daten[[3]]))
}