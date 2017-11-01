#using a script to run through the peak picking steps (XCMS), isotopes (CAMERA), and compMS2Miner to match MS2
#Running one sample at a time because I do not want compMS2Miner to work on a full set
#KL 8/24/2017

fullOutput_v2 <- function(oneFile, ionMode,msDataDir)
  
{
#most basic - only peak picking
library(xcms)
library(snow)
library(CAMERA)
library(compMS2Miner)
library(MetMSLine)
library(XLConnect)

#this seems to help issues with continually have to rinstall Rtools (8/2017)
#Sys.setenv(PATH="%PATH%;C:/RBuildTools/3.4/bin;c:/RBuildTools/3.4/ming3_32/bin")
library(devtools)
devtools::find_rtools()

##the first bit is just the usual XCMS steps for peak picking, grouping, and so on
# ## How many CPU cores has your machine (or cluster) ?
snowparam <- SnowParam(workers = 4, type = "SOCK")

#For negative ion mode: ppm = 2,3 seems best
#For positive ion mode: ppm = 3,4 seems best
xs<-xcmsSet(oneFile, method = "centWave",ppm = 4 ,snthresh = 10,
  prefilter = c(3,500), mzCenterFun = "wMean",integrate = 2,
  verbose.columns = FALSE, peakwidth=c(15,80), fitgauss= TRUE, noise = 500,
  mzdiff=-0.005,BPPARAM=snowparam) #,sleep = 0.00001)

#setting up the parameters
out <- as.data.frame(xs@peaks)
sOut <- out[order(out$mz),]

nSamples <- length(oneFile)
xsa<-xsAnnotate(xs)

#group the features initially just by retention time
xsaF <-groupFWHM(xsa)

#figure out which features also have a matching 13C feature. Have to enter both
#the relative error (ppm) and the absolute error (mzabs)
xsaFI <-findIsotopes(xsaF,ppm=1.5,mzabs = 0.0001,minfrac = nSamples)

#now group by the correlations based on (1) intensity, (2) EIC, (3) isotopes...
xsaC <-groupCorr(xsaFI,cor_eic_th=0.75,pval=0.05, graphMethod="hcs",
 calcIso = TRUE, calcCiS = TRUE, calcCaS = FALSE)

if (ionMode =='pos') {
  #setup the file to also look for adducts, only go with the primary adducts for the moment
  file <-system.file("rules/primary_adducts_pos.csv",package = "CAMERA")
  rules <-read.csv(file)
  an <-findAdducts(xsaC,polarity = "positive",rules=rules,ppm=1.5)
} else if (ionMode == 'neg') {
  #setup the file to also look for adducts, only go with the primary adducts for the moment
  file <-system.file("rules/primary_adducts_neg.csv",package = "CAMERA")
  rules <-read.csv(file)
  an <-findAdducts(xsaC,polarity = "negative",rules=rules,ppm=1.5)
  
} 

#do some housecleaning
rm(xsa,xsaF,xsaFI,xsaC)

#now set up the results in a format that compMS2miner can use
#Generate result
peakTable <- getPeaklist(an)
# rearrange columns
names(peakTable)[names(peakTable) == 'mz'] <- 'mzmed'
names(peakTable)[names(peakTable) == 'rt'] <- 'rtmed'

peakTable <- cbind(1:nrow(peakTable), peakTable)
colnames(peakTable)[1] <- "EICno"

colIndxTmp <- grep("EICno|mzmed|rtmed|adduct|into", colnames(peakTable))
peakTable <- peakTable[, c(colIndxTmp, setdiff(1:ncol(peakTable), colIndxTmp))]

#correct to (8/21/2017, KLongnecker)
peakTable[, c(1:3)] <- apply(peakTable[, 1:3], 2, as.numeric)
peakTable$EICno <- as.integer(peakTable$EICno)

#use this to change the file name into the one file I am working on (into is the peak intensity information, could also use intb)
changeN <- grep("into", colnames(peakTable))
colnames(peakTable)[changeN]<- basename(oneFile)
rm(changeN)

##get rid of unneeded columns bc causes issues later...only need first four columns
peakTable <- peakTable[,-(6:dim(peakTable)[2])]

# observation names
sampStr <- 'mtab_Method'

obsNames <- colnames(peakTable)[grep(sampStr, colnames(peakTable))]
#blanks <- colnames(peakTable)[grep(blankStr, colnames(peakTable))]
#qcs <- colnames(peakTable)[grep(qcStr, colnames(peakTable))]

##may turn these on later

# # zero fill
# peakTable <- MetMSLine::zeroFill(peakTable, c(obsNames, blanks))
# #peakTable <- MetMSLine::zeroFill(peakTable, c(obsNames))
# 
# # # blank fold change greater than 2
# # peakTable <- MetMSLine::blankSub(peakTable, blankNames=blanks, sampNames=obsNames,
# #                                  thresh=2)
# 
# # cv calculation using QC samples (retain features >= 30 %)
# peakTable <- MetMSLine::cvCalc(peakTable, qcs, thresh=30)
# # log transform
# peakTable <- MetMSLine::logTrans(peakTable, c(obsNames, blanks))



# compMS2Miner
# 1. Pre-processing: compMS2 object construction, noise filtration, composite spectra generation and substructure annotation
#directories/parameters

## *i.* compMS2 object construction from peak table and MS/MS files.


#note: even if the incoming peakTable only has one sample, the msDataDir portion of the command will find all possible mzML files
nCores <- parallel::detectCores()
compMS2_object <- compMS2Construct(MS1features=peakTable,
                                   msDataDir = usePath,
                                   MS2files = oneFile, #only match the file of interest
                                   nCores=nCores,
                                   mode = ionMode, precursorPpm = 10,
                                   ret = 20, TICfilter = 10000, verbose=TRUE) ##TRUE?



## *ii.* dynamic noise filter
compMS2_object <- deconvNoise(compMS2_object, "DNF")

## *iii.* intra-spectrum ion grouping and signal summing

compMS2_object <- combineMS2(compMS2_object, "Ions")

## *iv.* inter-spectrum ion grouping and signal summing
# LEAVE THIS ON...even if there is only one spectra, otherwise later steps will barf
compMS2_object <- combineMS2(compMS2_object, "Spectra", specSimFilter=0.8)


## *vi.* annotate possible substructures with internal database (?SubStructure_masses)

compMS2_object <- subStructure(compMS2_object, 'Annotate', Frag_mzabs=0.05)
# identify most probable substructure annotation based on total relative intensity
# explained
compMS2_object <- subStructure(compMS2_object, "prob")
# summary of most probable substructure annotation
mostProbSubStr <- subStructure(compMS2_object, "probSummary")

# 2. Combinatorial Metabolite Identification
subStrMassShift <- c(42.010565, 119.004101, 176.03209, 255.988909,
                     305.068159, 57.021464, 161.014666, 79.956817)
names(subStrMassShift) <- c("acetyl", "cysteine", "glucuronide",
                            "glucuronide sulfate", "glutathione", "glycine",
                            "mercapturate", "sulfate")

#  common mandatory esi adducts i.e. added to those detected in the MS1 data by CAMERA

if(ionMode == 'neg'){
  #mandEsiAdducts <- c('[M-H]-', '[M+Cl]-', '[M-H2O-H]-', '[M+Na-2H]-', '[M+K-2H]-',	'[M-H+HCOOH]-',	
   #                   '[M-H+CH3COOH+Na]-', '[M-H+CH3COOH]-')
  mandEsiAdducts <- c('[M-H]-', '[M+Cl]-')
} else if (ionMode == 'pos') {
  #mandEsiAdducts <- c('[M+H]+', '[M+NH4]+', '[M+Na]+','[M+CH3OH+H]+', '[M+CH3COONa]+',
  #                     '[M+K]+')
  mandEsiAdducts <- c("[M+H]+","[M+NH4]+")
  
  
}

#browser()
# what if any substructure detected?
subStrMassShift <- subStrMassShift[names(subStrMassShift) %in% names(mostProbSubStr)]
# internal database annotation, can also set metDB=LMSD
# compMS2_object <- metID(compMS2_object, "dbAnnotate", SubStrs=subStrMassShift,
#                            metDB=HMDB, includeElements=c('C', 'H', 'N', 'O', 'P',
#                                                          'S', 'F', 'Cl', 'Br',
#                                                          'I'),
#                            esiAdducts=mandEsiAdducts, MS1adducts=TRUE)

#having issues with MS1adducts=TRUE
compMS2_object <- metID(compMS2_object, "dbAnnotate", SubStrs=subStrMassShift,
                           metDB=HMDB, includeElements=c('C', 'H', 'N', 'O', 'P',
                                                         'S', 'F', 'Cl', 'Br',
                                                         'I'),
                           esiAdducts=mandEsiAdducts, MS1adducts=TRUE)

lbMspFiles <- paste0('https://raw.githubusercontent.com/WMBEdmands/mspFiles/master/MoNA-export-Libraries_-_LipidBlast_SMILES_', 1:5, '.msp')
for(i in 1:5){
  compMS2_object <- metID(compMS2_object, 'matchSpectralDB', mspFile=lbMspFiles[i]) 
}

## select most probable annotations based on substructures detected
compMS2_object <- metID(compMS2_object, "dbProb", minTimesId=1)

## predict Phase II metabolites from SMILES codes
compMS2_object <- metID(compMS2_object, "predSMILES")


#```{r metFrag, eval = TRUE}
###from notes: in silico fragmentation using command line version of MetFrag CL2.3
##turning this off as well because of the adduct issue with one file
## metFrag in silico fragmentation.
compMS2_object <- metID(compMS2_object, "metFrag",frag_mzabs=0.5)

#```{r CFM, eval = TRUE}
###CFM is a Wishart tool, described in 2015 in Metabolomics: Allen, F., Greiner, R., Wishart, D., Competitive Fragmentation Modeling of ESI-MS/MS spectra for putative metabolite identification, Metabolomics, 11:1, pp 98-110, 2015.
## CFM in silico fragmentation.
compMS2_object <- metID(compMS2_object, "CFM")


##export the pieces: annoying way to export two things from one function in R
fullOutput <- list("oneComp" = compMS2_object,"peakTable" = peakTable)

}
