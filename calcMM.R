# Function edited from MajorMinorCalc.R
# minDP = minimum depth of read coverage at a variant
# maxDP = maximum depth of read coverage at a variant
# minAF = minimum allele frequency

calcMM <- function(x, minDP = 20, maxDP = 1e+05, minAF = 0.2) {
  
  xf <- x[x$DP >= minDP,]
  xf <- xf[xf$DP <= maxDP,]
  AF1 <- xf$AD1/xf$DP
  AF2 <- xf$AD2/xf$DP
  xf <- data.frame("chr" = xf$chr, 
                   "start" = xf$position, 
                   "AF1" = AF1, 
                   "AF2" = AF2,
                   "DP" = xf$DP)
  ftx <- xf[xf$AF1 >= minAF,]
  ftx <- ftx[ftx$AF2 >= minAF,]
  sortedx <- Sort_major_minor(data=ftx, col1=3, col2=4)
  MajorMinor <- sortedx$AF1/sortedx$AF2
  sortedx["MajorMinor"] <- MajorMinor
  return(sortedx)
}
