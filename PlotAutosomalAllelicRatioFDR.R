# Edited from 'PlotGenome.R'.
# function to plot allelic ratios along human autosomes.

PlotAutosomalAllelicRatioFDR <- function(tab, window, Ylim) {
  
  require(zoo)      # ensure necessary packages loaded
  require(gplots)
  
  # set up graphics
  x11(width=17, height=5)
  par(oma=c(1, 1, 0, 10))
  
  # define chromosomes
  centromere_pos=c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,
                   17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7)
  centromere_pos=centromere_pos*1000000
  chr_size=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
             146364022,141213431,135534747,135006516,133851895,115169878,107349540,
             102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
  
  chr_total=0
  for (i in 1:(length(chr_size)-1)){
    chr_total=c(chr_total,sum(chr_size[1:i]))
  }
  
  genome_size = sum(chr_size)
  
  # ensure sex chromosomes omitted from vcf data
  tab <- tab[tab$chr %in% 1:22,]
  tab$start = tab$start + chr_total[tab$chr]
  
  # get colours for AR plot
  col = rollmedian(tab$chr, window)%%2 + 1
  col[col==1] = "dodgerblue4"
  col[col==2] = "dodgerblue1"
  
  # start plot
  plot(rollmedian(tab$start, 151), rollmedian(tab$MajorMinor, 151),
       col="dimgrey", pch=15, cex=0.4, xlim=c(1, genome_size), ylim=c(1, Ylim), type="l", xaxt = "n",
       xlab="chromosome", ylab="allelic ratio", cex.lab=1.5, bty='n')
  points(rollmedian(tab$start, 151),rollmedian(tab$MajorMinor, 151),
         col=col, pch=15, cex=0.4, ylim=c(1,Ylim))
  
  # calculate p.values of allelic ratios
  g <- tab$MajorMinor
  
  ttest <- function(x){
    ttt <- t.test(x, g, alternative = "greater")$p.value
    return(ttt)
  }
  
  # calculate p.values of AR skew
  tt <- rollapply(g, width=window, FUN=ttest)
  
  tt <- p.adjust(tt, "fdr")
  tt <- -1*log(tt,10)
  cl <- colorpanel(50000, "white","grey40")
  cl2 <- colorpanel(50000, "white","grey40", "red2")
  pmax <- ceiling(max(tt))
  if(pmax < 2) {pmax = 2}
  l <- length(rollmedian(tab$MajorMinor, window))
  xx <- rollmedian(tab$start, window)
  
  # plot out p.values of allelic ratios and annotate
  # regions with FDR < 0.01 highlighted in red
  for(i in 1:l){
    if (pmax == 2) {
      lines(c(xx[i], xx[i]), c(2.5,2.65), col=cl[round(tt[i]*(50000/pmax))])
    }
    if (pmax > 2) {
      lines(c(xx[i],xx[i]), c(2.5, 2.65), col=cl2[round(tt[i]*(50000/max(tt)))])
      if (tt[i] > 2) {lines(c(xx[i],xx[i]),c(2.7,2.85),col="red2", lwd=2)}
    }
  }
  
  xleft <- genome_size + 115000000
  xright <- xleft + 55000000
  ybottom <- 1.5
  ytop <- 2.5
  
  if (pmax > 2) {
    for(i in 1:100){
      rect(xleft, ybottom, xright, ytop, xpd=T, col=cl2[i*500], border = NA )
      ybottom = ybottom + 0.01
    } 
    ybottom <- 1.5
    rect(xleft, ybottom, xright, ytop, xpd=T, border = "grey28" )
    text(xright, 2.5, pos=4, "-Log10(P-Value)", font=2, col='grey10', xpd=T, srt=-90)
    text(xleft + 27500000, 2.5, pos=3, pmax, font=2, col='grey10', xpd=T)
    text(xleft + 27500000, 1.5, pos=1, "0", font=2, col='grey10', xpd=T)
  }
  
  if (pmax <= 2) {
    for(i in 1:100){
      rect(xleft, ybottom, xright, ytop, xpd=T, col=cl[i*500], border = NA )
      ybottom = ybottom + 0.01
    }
    ybottom <- 1.5
    rect(xleft, ybottom, xright, ytop, xpd=T, border = "grey28" )
    text(xright, 2.5, pos=4, "-Log10(P-Value)", font=2, col='grey10', xpd=T, srt=-90)
    text(xleft + 27500000, 2.5, pos=3, "2", font=2, col='grey10', xpd=T)
    text(xleft + 27500000, 1.5, pos=1, "0", font=2, col='grey10', xpd=T)
  }
  
  # plot out chromosomes
  for(i in 1:22) {
    if (i>1){lines(c(chr_total[i], chr_total[i]), c(1,Ylim), col="gray48", lwd=0.8)}
    lines(c(chr_total[i] + centromere_pos[i], 
            chr_total[i] + centromere_pos[i]), c(1, Ylim), col="gray55", lty=4, lwd=0.8)
    text(chr_total[i]+chr_size[i]/2, 0.8, i, xpd = TRUE, cex=0.85, font=2)
  }
  
  sum=0
  for (i in 1:(length(chr_size))){    
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=10,col="black")  
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=8,col="gray50")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=7,col="gray53")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=6,col="gray59")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=4,col="gray75")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=2,col="gray85")
    lines(c(sum+20000000,sum+chr_size[i]-20000000),c(1,1),lwd=1,col="gray90")
    lines(c(sum+centromere_pos[i],sum+centromere_pos[i]),c(1.01,0.99),col="grey13",lwd=2)
    
    sum=sum+chr_size[i]
  }
  # clear up memory
  rm(tab, i)
}
