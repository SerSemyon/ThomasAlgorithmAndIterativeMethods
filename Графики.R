lf<-list.files(pattern = '\\.csv')


for (i in 1:length(lf)){
  x<-read.csv2(lf[i])
  if(ncol(x)==2){
    png(filename=paste0(lf[i],".png"))
    plot(x[,1],x[,2], xlab = "Êîëè÷åñòâî èòåğàöèé", ylab = "Ïîãğåøíîñòü")
    lines(x[,1],x[,2])
    dev.off()
  }
}
for (i in 1:length(lf)){
  x<-read.csv2(lf[i])
  if(ncol(x)==4){
    png(filename=paste0(lf[i],".png"))
    plot(x[,1],x[,3], xlab = "x", ylab = "y")
    lines(x[,1],x[,2])
    dev.off()
  }
}

x<-read.csv2("Îìåãà.csv")
png(filename=paste0("Îìåãà.png"))
plot(x[,1],x[,2], xlab = "Îìåãà", ylab = "Êîëè÷åñòâî èòåğàöèé")
lines(x[,1],x[,2])
dev.off()
