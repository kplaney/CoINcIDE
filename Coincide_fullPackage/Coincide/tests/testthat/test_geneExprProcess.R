
context("Gene Expression Processing")

test_that("tesing the filterAndImputeSamples function", {
  
  numSamples <- 30
  numGenes <- 100
  tmpExpr <- matrix(data=rnorm(numSamples*numGenes,mean=2,sd=0.01),ncol=numSamples,nrow=numGenes,
                    dimnames=list(paste0("gene_",c(1:numGenes)),paste0("sample_",c(1:numSamples))))
  
  tmpData <- list()
  tmpData$expr <- tmpExpr
  tmpData$keys <- rownames(tmpExpr)
  #make one "patient" all NAs
  tmpData$expr[,1] <- NA
  #make one gene all NAs
  tmpData$expr[50,1:30] <- NA
  #make some random NAs
  tmpData$expr[10,10] <- NA
  
  output <- filterAndImputeSamples(study=tmpData,studyName = "study",
                                   outputFile = "createTestTrainSetsOutput.txt",
                                   impute=TRUE, knnFractionSize=.01,
                                   fractionSampleNAcutoff=(10/numSamples),
                                   fractionGeneNAcutoff = (10/numGenes),
                                   exprIndex="expr",sampleCol=TRUE,
                                   returnErrorRate=TRUE)
  
  expect_equal(ncol(output$exprFilterImpute), numSamples-1)
  expect_equal(nrow(output$exprFilterImpute), numGenes-1)
  expect_true(is.numeric(output$exprFilterImpute))

})
