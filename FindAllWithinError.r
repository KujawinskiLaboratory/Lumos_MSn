#write up a version of FindAllWithinError for R
#based on SimStitch code, modified in 2007 to find all mz values within the error set
#KLongnecker, 23 August 2017
FindAllWithinError <- function(oneList,maxErr,target)

##to run this...do this in the console:
#source("FindAllWithinError.r")
#r <- FindAllWithinError(oneList,maxErr,target)

# oneList <- peakTable$mzmed
# maxErr <- 1
# #target <- 196.076276
# target <- 172.9564
{
#sort input list
inputList <- sort(oneList)
idxSortInput <- sort(oneList,index.return=TRUE)$ix


if (length(target)==1) {
  ti <- 1
  } else {
  stop('This version only handles one mz value')}

e = abs((inputList - target[ti])/target[ti])
ii <- which(e <=maxErr/1e6)

FindAllWithinError <- idxSortInput[ii]
}