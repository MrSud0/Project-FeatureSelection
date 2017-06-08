


#this function will take as input a transposed mat and a factor with genoData
#it will return the selected genes on a matrix based on rfe feature selection using function crossval

FSelection <- function (esetMatT,genoData) {
  
  #finds the most important features/genes
  res.rfe <- crossval(esetMatT,genoData,DEBUG=TRUE, theta.fit=fit.rfe)
  
  #we save the most important features/genes on a matrix 
  x <- as.matrix(extractFeatures(res.rfe,toFile=FALSE))
  
  
  temp <- colnames(esetMatT)  #columns' names for if statements below
  
  selectedMat <- 0 #this will be the selected genes matrix
  
  #checks which genes are the important ones and creates a matrix with only those
  for( y in 1  : nrow(x))
  {
    for (i in 1 :ncol(esetMatT))
    {
      if(temp[i]==x[y] & y == 1)
      {selectedMat<-as.matrix(esetMatT[,i])
      break} 
      if(temp[i]==x[y] & y != 1)
      {selectedMat<-cbind(selectedMat,esetMatT[,i]) 
      break}
    }
  }
  
  #we put the names of the genes to the columns
  colnames(selectedMat)<-t(x[,1]);
  
  return(selectedMat)
  
}