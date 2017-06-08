################Fisher Score Feature Selection on normalized data #############


FisherFS<-function(esetMatT,genoData)
{
  genoDataMat <- as.matrix(genoData)
 
  
  score <- 0
  selectedMat <- 0
  finalMat <- 0
  scoreSorted <- 0
  
  for(y in 1 : ncol(esetMatT))   
  {
    matPlus <- 0
    matMinus <- 0
    
    for(i in 1 : nrow(genoDataMat))
    { 
      if(genoDataMat[i] == "ER (IHC): -")
      {matMinus <- rbind(matMinus,esetMatT[i,y]) }
      else
      {matPlus <- rbind(matPlus,esetMatT[i,y]) }
      
      
    }
    
    score[y] <- ( mean(matPlus[2:length(matPlus)]) - mean(matMinus[2:length(matMinus)]))  ^ 2
    names(score[y])
  }
  
score <- as.data.frame(score)
tempNames <- as.matrix(names(esetMatT[1,]))

score <- cbind(score,tempNames)
scoreSorted <- score[mixedorder(score$score,decreasing = TRUE),]

selectedMat <- scoreSorted[1:98,]
flag = FALSE

for ( y in 1 : nrow(tempNames))
{
  for (i in 1: nrow(selectedMat))
  {
    if( tempNames[y] == selectedMat[i,2] & flag == FALSE)
     {finalMat <- esetMatT[,y]
      flag=TRUE
      break}
    if( tempNames[y] == selectedMat[i,2] ) 
    { finalMat <- cbind(finalMat,esetMatT[,y])
    break }
      
      
  }
  
}
colnames(finalMat) <- selectedMat[,2]  
 
  

return(finalMat)


}