# short function to sort MajorMinor calls

Sort_major_minor <- function(data, col1, col2){
  for (i in 1:dim(data)[1]) {
    if (data[i,col1] < data[i,col2]) {
      save <- data[i,col2]
      data[i,col2] <- data[i,col1]
      data[i,col1] <- save
    }
  }
  return(data)
}
