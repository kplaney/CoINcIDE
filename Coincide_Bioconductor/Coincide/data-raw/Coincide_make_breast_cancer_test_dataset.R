

#make test dataset
library("curatedBreastData")
#load up the full datasets with real gene expression data
data(curatedBreastDataExprSetList)

#know dataset #5 (GSE9983) has a few NAs. cherry-pick these out to test NA
#impute function.
indices <- c()
for(g in 1:nrow(exprs(curatedBreastDataExprSetList[[5]]))){
  if(any(is.na(exprs(curatedBreastDataExprSetList[[5]])[g,]))){
    indices <- append(indices,g)
  }
}

#so it's patient 127 with 3 NAs:
which(is.na(exprs(curatedBreastDataExprSetList[[5]])[indices[1],]))
#147
which(is.na(exprs(curatedBreastDataExprSetList[[5]])[indices[2],]))
#147
which(is.na(exprs(curatedBreastDataExprSetList[[5]])[indices[3],]))
#147


#OK now that we know this is a good test dataset for a processing function,
#let's take a small subset of it
#add on some extra non-NA samples and indices
#actually use the Coincide function to create an S4 expression set
library("Coincide")
exprData <- exprs(curatedBreastDataExprSetList[[5]][append(c(1:20), indices),
                  c(147, 148, 149, 150)])
#in phenoData: patients/samples are in the rows
phenData <- pData(curatedBreastDataExprSetList[[5]])[c(147, 148, 149, 150), ]
featData <- rownames(exprData)
featData <- data.frame(featData)
colnames(featData) <- "gene"
testData1 = createS4exprSet(expr = exprData, phenoData = phenData, 
                            featureData = featData, featureDataFieldName = 
                              "gene")

#now take another dataset - the first one - and subset it
#just take first 25 genes, first 5 patients
exprData <- exprs(curatedBreastDataExprSetList[[1]][c(1:25), c(1:5)])
phenData <- pData(curatedBreastDataExprSetList[[1]])[c(1:5), ]
featData <- rownames(exprData)
featData <- data.frame(featData)
colnames(featData) <- "gene"
testData2 = createS4exprSet(expr = exprData, phenoData = phenData, 
                            featureData = featData, featureDataFieldName = 
                              "gene")

breastExData <- list(data1 = testData1, data2 = testData2)
#then save as RData object
savePath <- "~/Dropbox/github/Coincide/Coincide_Bioconductor/Coincide/data/"
save(breastExData, file = paste0(savePath, "/breastExData.RData"))

