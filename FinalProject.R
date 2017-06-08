##################Project###########################


#Importing Libraries







#Downloads "GDSXXXX" and saves it in a global variable GEOdata
GDSdownload("GDS4057") 
#converting GDS to express set
eset<-GDS2eSet(GEOdata,do.log2=TRUE) 
#we take only the samples needed
genoData<-eset@phenoData@data$`genotype/variation` 
esetReduced<-eset[,genoData %in% c("ER (IHC): -","ER (IHC): +")]

esetMat<-(exprs(esetReduced))

#we create the factor genoData for crossval
genoData<-factor(genoData[genoData %in% c("ER (IHC): -","ER (IHC): +")]) 

#we transpose the matrix for the crossval input
esetMatT<-t(esetMat)


#we replace the NaN values with zeroes to avoid nan errors
esetMatT[!is.finite(esetMatT)] <- 0 




#this function will take as input a transposed mat and a factor with genoData
#it will return the selected genes on a matrix based on rfe feature selection using function crossval

crossValMat <- FSelection(esetMatT,genoData)
#printing a heatmap of the selected features
heatmap(CrossValMat)

#this function will take as input a transposed mat nad a factor with genoData
#it will return the selected genes on a matrix based on the Fisher Score
fisherMat <- FisherFS(esetMatT,genoData)

#printing a heatmap of the selected features
heatmap(FisherMat)

#control variable used by svm model
control<-trainControl(method = "repeatedcv",number = 10, repeats = 3)

#we create and train out a SVM model based on our crossValMat and the genoData

modelSvm1<-train(x=CrossValMat,y=genoData ,method="svmLinearWeights",trControl = control)

#we create and train out a SVM model based on our fisherMat and the genoData
modelSvm2<-train(x=FisherMat,y=genoData ,method="svmLinearWeights",trControl = control)


#predict data for GDS4057 using the model trained on the crossval selected data (basically testing the model with the data used to train it)
pr <- predict(modelSvm1, newdata = esetMatT)
#create confusion matrix and calculate accuracy 
CMAdjusted(pr,genoData)
View(cm1)
View(accuracy)

#predict data for GDS4057 using the model trained on the Fisher Score selected data (basically testing the model with the data used to train it)
pr2 <- predict(modelSvm2, newdata = esetMatT)
#create confusion matrix and calculate accuracy 
CMAdjusted(pr2,genoData)
View(cm1)
View(accuracy)


##testing the second dataset 4056 ##
GDSdownload("GDS4056")
eset2<-GDS2eSet(GEOdata,do.log2=TRUE)
genoData2<-eset2@phenoData@data$`genotype/variation`
#we dont need to reduce the samples on this one as they dont contain unwanted values
esetMat2<-exprs(eset2)
esetMat2T<-t(esetMat2)
esetMat2T[!is.finite(esetMat2T)] <- 0

#predict data for GDS4056 using the model trained on the crossval selected data
pr3 <-predict(modelSvm1, newdata = esetMat2T)
#create confusion matrix and calculate accuracy  
CMAdjusted(pr3,genoData)
View(cm1)
View(accuracy)
####predict data for GDS4056 using the model trained on the Fisher Score selected data
pr4 <-predict(modelSvm2, newdata = esetMat2T)
CMAdjusted(pr4,genoData)
View(cm1)
View(accuracy)

#we print the names of the selected features ,crossval, so we can use them for the venn diagram
write.table(names(crossValMat[1,]), file="mymatrix.txt", row.names=FALSE, col.names=FALSE)
#we print the names of the selected features ,crossval, so we can use them for the venn diagram
write.table(names(fisherMat[1,]), file="mymatrix2.txt", row.names=FALSE, col.names=FALSE)
