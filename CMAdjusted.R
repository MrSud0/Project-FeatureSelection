###Confusion Matrix + accuracy for svm model adjusted ##
###predicted is the data we received from the predict method
###data is the original data (genoData)

CMAdjusted <- function (predicted,data){

#we transform the data to matrix for easier use 
predicted <- as.matrix(predicted)
data <- as.matrix(data)


cm1<-0


#replace +/- classes to match our model class with an if statement to reduced workload

  

if(data[1] == "ER -" || data[1] == "ER +" ){
  for(i in 1:length(data))
   {
  
  if(data[i] == "ER -")
   {data[i] <- "ER (IHC): -"}
  else
   {data[i] <- "ER (IHC): +"}

  }
}

if(data[1] == "ER- breast cancer" || data[1] == "ER+ breast cancer" ){
  for(i in 1:length(data))
    {
  
  if(data[i] == "ER- breast cancer")
    {data[i] <- "ER (IHC): -"}
  else
    {data[i] <- "ER (IHC): +"}
  
  }
}

##we gather values for our matrix 
cm1<- matrix(0,2,2 ,dimnames = list(  c("ER +", "ER -"),c("ER +" ,"ER -" ) ) )
##cm[1] TP
##cm[2] FN
##cm[3] FP
##cm[4] TN

for(i in 1:length(predicted))
{
  if(predicted[i] == "ER (IHC): +" & predicted[i] == data[i]){
    cm1[1] <- cm1[1]+1
  }
  if(predicted[i] == "ER (IHC): +" & predicted[i] != data[i]){
    cm1[3] <- cm1[3]+1
  }
  if(predicted[i] == "ER (IHC): -" & predicted[i] == data[i]){
    cm1[4] <- cm1[4]+1
  }
  if(predicted[i] == "ER (IHC): -" & predicted[i] != data[i]){
    cm1[2] <- cm1[2]+1
  }
  
   
} 

  accuracy <- (cm1[1]+cm1[4])/(cm1[1]+cm1[4]+cm1[3]+cm1[2])

  cm1 <- cbind(cm1,rbind("accuracy",accuracy))
  
  return(cm1)
    
  
}
    
  
  
  
