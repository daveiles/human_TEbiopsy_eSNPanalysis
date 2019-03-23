# PZSautosomes function edited from 'Plot_Zygosity_Single.R' to plot each SNP, 
# without any summarization. Function limited to human autosomes.

PZSautosomes <- function(tab){
  require(zoo)      # ensure necessary packages loaded
  require(gplots)
  
  tbl=tab
  
  {
    centromere_pos=c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,
                     17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7)
    centromere_pos=centromere_pos*1000000
    chr_size=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
               146364022,141213431,135534747,135006516,133851895,115169878,107349540,
               102531392,90354753,81195210,78077248, 59128983,63025520,48129895,51304566)
    lb=c(seq(from = 150, to = 10, by = -20),seq(from = 10, to = 150, by = 20))
    plot(c(1:22),rep(0,22),ylim=c(-150000000,150000000),cex.axis=0.7,xlab = "Chromosome",
         ylab="Position (Mbp)",xaxt = "n",yaxt="n")
    axis(2, seq(from = 150000000, to = -150000000, by = -20000000), labels=lb, cex.axis=0.7, las=1)
    mx=22
  }
  
  for (i in 1:mx) {
    chr <- i
    tbl2 <- tbl[tbl$chr==chr,]
    trsnp <- tbl2[tbl2$snpType == 1,]
    trsnp0 <- tbl2[tbl2$snpType == 0,]
    fl <- tbl2[tbl2$snpType == -1,]
    
    if (dim(fl)[1] > 0){
      for (j in 1:dim(fl)[1]){
        x = centromere_pos[i] - fl$start[j]
        lines(c(i,i - 0.4), c(x,x), col="DodgerBlue", lwd=0.1)
      }
    } 
    
    if(dim(trsnp)[1] > 0){
      for (j in 1:dim(trsnp)[1]){
        x = centromere_pos[i] - trsnp$start[j]
        lines(c(i,i + 0.4), c(x,x), col="Red", lwd=0.1)
      }
    } 
    
    if (dim(trsnp0)[1] > 0){
      for (j in 1:dim(trsnp0)[1]){
        x = centromere_pos[i] - trsnp0$start[j]
        lines(c(i,i + 0.4), c(x,x), col="Magenta", lwd=0.1)
      }
    }
    
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=10,col="black")  
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=8,col="gray50")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=7,col="gray53")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=6,col="gray59")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=4,col="gray75")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=2,col="gray85")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=1,col="gray90")
    
    points(i,0,pch=16,col="grey13")
    if (i<23) {text(i,centromere_pos[i]+18000000,i,cex=1)}
  }
}
