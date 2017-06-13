
#this function checks if following packages are loaded , if not , it will download them and load them

DLP<-function()
{

if (!require("GEOquery")) {
  source("http://www.bioconductor.org/biocLite.R")
  biocLite("GEOquery")
  library(GEOquery)
}

if (!require("Biobase")) {
  source("http://www.bioconductor.org/biocLite.R")
  biocLite("Biobase")
  library(Biobase)
}


if (!require("gtools")) {
  install.packages("gtools", dependencies = TRUE)
  library(gtools)
}

if (!require("pathClass")) {
  install.packages("pathClass", dependencies = TRUE)
  library(pathClass)
}
if (!require("caret")) {
  install.packages("caret", dependencies = TRUE)
  library(caret)
}
if (!require("e1071")) {
  install.packages("e1071", dependencies = TRUE)
  library(e1071)
}


}