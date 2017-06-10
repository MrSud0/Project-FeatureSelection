################Fisher Score Feature Selection on normalized data #############

#####this function requires gtools package ######
FisherFS<-function(esetMatT,genoData)
{
  genoDataMat <- as.matrix(genoData)
 
  
  score <- NULL
  selectedMat <- NULL
  finalMat <- NULL
  scoreSorted <- NULL
  
  
  for(y in 1 : ncol(esetMatT))   
  {
    matPlus <- NULL
    matMinus <- NULL
    
    for(i in 1 : nrow(genoDataMat))
    { 
      if(genoDataMat[i] == "ER (IHC): -")
      {matMinus <- rbind(matMinus,esetMatT[i,y]) }
      else
      {matPlus <- rbind(matPlus,esetMatT[i,y]) }
      
      
    }
    
    score[y] <- ( mean(matPlus[1:length(matPlus)]) - mean(matMinus[1:length(matMinus)]))  ^ 2
    
  }
  
score <- as.data.frame(score)
tempNames <- as.matrix(names(esetMatT[1,]))

score <- cbind(score,tempNames)
scoreSorted <- score[mixedorder(score$score,decreasing = TRUE),]

selectedMat <- scoreSorted[1:100,]


for ( y in 1 :  nrow(tempNames))
{
  for (i in 1: nrow(selectedMat))
  {
    if( tempNames[y] == selectedMat[i,2] ) 
    { finalMat <- cbind(finalMat,esetMatT[,y])
      break }
      
      
  }
  
}
colnames(finalMat) <- selectedMat[,2]  
 
  

return(finalMat)


}