### Portfolio Script to investigate range of simulation outcomes
###Kiri Daust, Nov 2020

require(foreach)
require(reticulate)
require(Rcpp)
library(gridExtra)
library(data.table)
library(scales)
library(tidyr)
library(magrittr)
library(ggthemes)
require(rgdal)
require(tidyverse)
require(Rmisc)

##Set drive with cloud data
if(dir.exists("C:/users/whmacken/Sync")){
  cloud_dir <- "C:/users/whmacken/Sync/CCISS_data/"
}else{
  cloud_dir <- "./inputs/"
}

## source functions
sourceCpp("./CppFunctions/SimGrowth.cpp")
source_python("./PythonFns/PortfolioOptimisation.py")
source("CCISS_Fns.R")

###load data
inputDatName <- "PortPoints_Quesnel_90 GCMsMSY.csv"
SuitTable <- fread("InputsGit/Feasibility_v11_21.csv") ##tree spp suitability
SuitTable <- unique(SuitTable)
SuitNew <- fread("InputsGit/Feas_toAdd.csv")

colnames(SuitTable)[2:4] <- c("SS_NoSpace","Spp","Suitability")
SuitTable <- SuitTable[,c("BGC","SS_NoSpace","Spp","Suitability")]
SuitTable <- rbind(SuitTable, SuitNew)

SIBEC <- fread("InputsGit/PredSI_May2020.csv") 
SIBECnew <- fread("InputsGit/SI_to_add.csv")
SIBEC <- rbind(SIBEC, SIBECnew)
###import SI data (currently from BART)
SIBEC <- SIBEC[,c("SS_NoSpace","Spp","SIPred")] %>% set_colnames(c("SS_NoSpace","TreeSpp","MeanPlotSiteIndex"))

eda <- fread("InputsGit/Edatopic_v11_20.csv")
eda <- unique(eda[is.na(Special),.(BGC,SS_NoSpace, Edatopic)])

###run cciss predict
load(paste0(cloud_dir, "WNAv11_35_VAR_SubZone_ranger.Rdata"))##BGC model
Edatope <- fread("./InputsGit/Edatopic_v11_20.csv",data.table = T)
rawDat <- fread(paste0(cloud_dir,inputDatName),data.table = T)
CCISSPred <- CCISS_Spp(Y1 = rawDat,BGCmodel = BGCmodel,E1 = as.data.table(Edatope))

###rename and cleanup
SSPredOrig <- as.data.table(CCISSPred[[1]])
SSPredOrig[,allOverlap := NULL]
setnames(SSPredOrig, old = c("BGC","SS_NoSpace","SS.pred"), 
         new = c("MergedBGC","SSCurrent","SS_NoSpace"))
SSPredOrig <- SSPredOrig[,.(MergedBGC,SS_NoSpace,SSratio,SSprob,SSCurrent,FuturePeriod,SiteNo)]

cleanData <- function(SSPredAll,SIBEC,SuitTable,SNum,Trees,timePer,selectBGC){
  SSPred <- SSPredAll[SiteNo == SNum,] ###subset
  
  ##Merge SIBEC data
  SIBEC <- SIBEC[TreeSpp %in% Trees,]
  SSPred <- SSPred[SSPred$FuturePeriod %in% timePer,]
  SSPred <- SSPred[,.(FuturePeriod, SSCurrent,SS_NoSpace, SSprob)]
  SSPred <- merge(SSPred, SIBEC, by = "SS_NoSpace", all.x = TRUE)
  
  
  ###Add rows for species with missing SI - mostly US units here
  add <- foreach(Year = unique(SSPred$FuturePeriod), .combine = rbind) %do%{
    byYear <- SSPred[SSPred$FuturePeriod == Year,]
    foreach(SS = unique(byYear$SS_NoSpace), .combine = rbind) %do%{
      bySS <- byYear[byYear$SS_NoSpace == SS,]
      missing <- Trees[!Trees %in% bySS$TreeSpp]
      new <- bySS[rep(1,length(missing)),]
      new$TreeSpp <- missing
      new
    }
  }
  if(!is.null(add)){
    if(nrow(add) > 0){
      add$MeanPlotSiteIndex <- 5 ##Set missing SI
      SSPred <- rbind(SSPred, add)
    }
  } 
  
  
  SSPred <- SSPred[!is.na(SSPred$TreeSpp),]
  setnames(SSPred,old = "TreeSpp", new = "Spp")
  
  ##Add suitability
  SSPred <- merge(SSPred, SuitTable, by = c("SS_NoSpace","Spp"), all.x = TRUE)
  SSPred$Suitability[is.na(SSPred$Suitability)] <- 5
  
  temp <- SIBEC[SIBEC$SS_NoSpace == selectBGC,]
  if(nrow(temp) == 0){
    return(NULL)
  }
  ###Create current data
  current <- temp %>% 
    merge(SuitTable, by.x = c("TreeSpp","SS_NoSpace"), by.y = c("Spp","SS_NoSpace"), all.x = TRUE) %>%
    unique()
  
  ###check that there aren't errors in the table
  temp <- aggregate(SS_NoSpace ~ TreeSpp, current, FUN = length)
  if(any(temp$SS_NoSpace > 1)){
    stop("There are partial duplicates in the suitablity table. Please fix them. :)")
  }
  
  current <- data.frame(Spp = current$TreeSpp, FuturePeriod = 2000, 
                        MeanSI = current$MeanPlotSiteIndex, MeanSuit = current$Suitability)
  
  missing <- Trees[!Trees %in% current$Spp]
  if(length(missing) > 0){
    new <- current[rep(1,length(missing)),]
    new$Spp <- missing
    new$MeanSI <- 10
    new$MeanSuit <- 5
    current <- rbind(current, new)
  }
  
  ##Summarise data- average SI and Suit weighted by SSProb
  SS.sum <- SSPred[,.(MeanSI = sum(MeanPlotSiteIndex*(SSprob/sum(SSprob))),
                      MeanSuit = round(sum(Suitability*(SSprob/sum(SSprob))), digits = 0)),
                   by = .(Spp,FuturePeriod)]
  
  SS.sum <- rbind(SS.sum, current)
  SS.sum <- SS.sum[SS.sum$FuturePeriod %in% timePer,]
  ###not sure what we were doing here?
  SS.sum$MeanSI[SS.sum$MeanSuit == 4] <- 5
  SS.sum$MeanSI[SS.sum$MeanSuit == 5] <- 0
  SS.sum$MeanSuit[SS.sum$MeanSuit == 5] <- 4
  SS.sum <- unique(SS.sum)
  SS.sum <- SS.sum[order(SS.sum$Spp,SS.sum$FuturePeriod),]
  #SS.sum$MeanSuit[is.na(SS.sum$MeanSuit)] <- 3
  return(SS.sum)
}

edatopicSubset <- function(SSPredOrig, eda, pos = "Zonal"){
  if(pos == "Zonal"){
    SSPredFull <- SSPredOrig[grep("01",SSPredOrig$SSCurrent),]
  }else{
    edaSub <- eda[Edatopic == pos,]
    SSPredFull <- SSPredOrig[SSCurrent %in% edaSub$SS_NoSpace,]
  }
  SSPredFull$BGC_analysis <- gsub("/.*","", SSPredFull$SSCurrent)
  
  SSPredFull <- SSPredFull[,c("MergedBGC", "SS_NoSpace", "SSprob", "SSCurrent", 
                              "FuturePeriod", "SiteNo","BGC_analysis")]
  
  ##remove cases where not all timeperiods available
  SSPredFull <- as.data.table(SSPredFull)
  temp <- SSPredFull[,.(Num = length(unique(FuturePeriod))), by = c("SiteNo","SSCurrent")]
  temp <- temp[Num == 3,-c("Num")]
  SSPredFull <- SSPredFull[temp,on = c("SiteNo","SSCurrent")]
  SSPredFull[,SiteNo := as.numeric(SiteNo)]
}

SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.35,1,1.8,3), 
                       "NoMort" = c(70,60,50,30), "RuinSeverity" = c(0.4,0.5,0.7,0.8)) ####ProbDead- out of 100 trees, how many will die each year at each suitability. NoMort- Percent of time no mortality

minAccept <- 0.01 ##min acceptable weight in porfolio - if lower, will remove and re-optimize

cols <- fread("./InputsGit/PortfolioSppColours.csv")
cols <- cols[HexColour != "",]
myPal <- cols$HexColour
names(myPal) <- cols$TreeCode
Trees <- unique(SIBEC$TreeSpp)
myColours <- data.table(TreeCode = Trees)
myColours <- cols[myColours, on = "TreeCode"]
myColours <- myColours[!is.na(HexColour),]
pal <- myColours$HexColour
names(pal) <- myColours$TreeCode
colScale <- scale_fill_manual(name = "variable", values = pal)
Trees <- myColours$TreeCode

##set min and max weights for each species
boundDat <- data.table(Spp = Trees)
boundDat[,`:=`(minWt = 0, maxWt = 1)]

timePeriods = c(2000,2025,2055)
returnValue = 0.9

###subset for zonal or other
SSPredFull <- edatopicSubset(SSPredOrig,eda,pos = "Zonal")

treeList <- Trees
nSpp <- length(Trees)
Units <- unique(SSPredFull$BGC_analysis)

###SELECT BGC
BGC = Units[2]
  
SSPredBGC <- SSPredFull[BGC_analysis == BGC,-("BGC_analysis")]
SSList <- unique(SSPredBGC$SSCurrent)
selectBGC = SSList[1]
SSPredAll <- SSPredBGC[SSPredBGC$SSCurrent == selectBGC,]

SiteList <- unique(SSPredAll$SiteNo)
#SiteList <- rep(SiteList, each = round(15/length(SiteList)))
SSPredAll <- SSPredAll[SSPredAll$SiteNo %in% SiteList & !is.na(SSPredAll$SSprob),]

SL <- SiteList
allSitesSpp <- foreach(SNum = SL, .combine = rbind, 
                       .packages = c("reshape2","Rcpp","magrittr","reticulate",
                                     "data.table","foreach","ggplot2","ggthemes","scales"), 
                       .noexport = c("gs2gw", "simGrowthCBST","simGrowthCpp"),
                       .export = c("cleanData","loopCombine")) %do% {
                         reticulate::source_python("./PythonFns/PortfolioOptimisation.py") 
                         
                         cat("Optimising site",SNum,"...\n")
                         SS.sum <- cleanData(SSPredAll,SIBEC,SuitTable,SNum, Trees, 
                                             timePer = timePeriods,selectBGC = selectBGC)
                         if(any(is.na(SS.sum$MeanSuit))){
                           warning("Missing Suitability in unit ",
                                   BGC,", sitenumber ",SNum," for ",
                                   SS.sum$Spp[is.na(SS.sum$MeanSuit)], ": They will be filled with suit = 4")
                           SS.sum$MeanSuit[is.na(SS.sum$MeanSuit)] <- 4
                         }
                         SS.sum$FuturePeriod <- as.numeric(SS.sum$FuturePeriod)
                         if(length(timePeriods) == 1){
                           temp <- SS.sum
                           temp$FuturePeriod <- SS.sum$FuturePeriod[1]+85
                           SS.sum <- rbind(SS.sum, temp)
                         }
                         
                         if(!is.null(SS.sum)){
                           annualDat <- data.frame("Year" = seq(2000,2100,1))
                           output <- data.frame("year" = annualDat$Year)
                           
                           for (k in 1:nSpp){ ##for each tree
                             DatSpp <- SS.sum[SS.sum$Spp == treeList[k],]
                             dat <- data.frame("Period" = rescale(as.numeric(DatSpp$FuturePeriod), 
                                                                  to = c(2000,2085)), 
                                               "SIBEC" = DatSpp$MeanSI, "Suit" = DatSpp$MeanSuit)
                             dat$SIBEC <- dat$SIBEC/50 ##for mean annual increment
                             dat <- merge(dat, SuitProb, by = "Suit")
                             s <- approx(dat$Period, dat$SIBEC, n = 101) ##Smooth SI
                             p <- approx(dat$Period, dat$ProbDead, n = 101) ###Smooth Prob Dead
                             m <- approx(dat$Period, dat$NoMort, n = 101) ##Smooth No Mort
                             ##r <- approx(dat$Period, dat$RuinSeverity, n = 101)
                             
                             ###data frame of annual data
                             annualDat <- data.frame("Year" = seq(2000,2100,1), "Growth" = s[["y"]], 
                                                     "MeanDead" = p[["y"]], "NoMort" = m[["y"]]) ##create working data
                             
                             Returns <- simGrowthCpp(DF = annualDat)
                             tmpR <- c(0,Returns)
                             assets <- Returns - tmpR[-length(tmpR)]
                             temp <- data.frame(Spp = rep(treeList[k],101), 
                                                Year = 1:101, Returns = Returns)
                             output <- cbind(output, assets)
                           } ## for each tree species
                           
                           colnames(output) <- c("Year", treeList)
                           
                           ####Portfolio#######################################
                           returns <- output
                           rownames(returns) <- returns[,1]
                           returns <- returns[,-1]
                           ###only include species with mean return > 1 in portfolio
                           use <- colnames(returns)[colMeans(returns) > 1] ###should probably be higher
                           returns <- returns[,use]
                           sigma2 <- as.data.frame(cor(returns)) ###to create cov mat from returns
                           
                           ef <- ef_weights_v2(returns, sigma2, boundDat,minAccept) 
                           ef_w <- ef[[1]]
                           ef_w$Sd <- ef[[2]]
                           ef_w$Return <- 1:20
                           ef_w$RealRet <- ef[[3]]
                           ef_w$Sharpe <- ef[[4]]
                           
                           eff_front2 <- as.data.table(ef_w)
                           eff_front2[,RealRet := RealRet/max(RealRet)]
                            
                           ##now use this portfolio for 100 simulations to get volume
                           risk <- eff_front2[20,!c("Sd","Return","Sharpe","RealRet")] 
                           safe <- eff_front2[Sharpe == max(Sharpe),!c("Sd","Return","Sharpe","RealRet")]
                           v90 <- eff_front2[which.min(abs(RealRet - 0.9)),!c("Sd","Return","Sharpe","RealRet")]
                           currTrees <- colnames(safe)
                           
                           simVolume <- foreach(i = 1:100, .combine = rbind) %do% {
                             output <- data.table(Spp = character(),Vol = numeric(),nTree = numeric())
                             for(k in 1:length(currTrees)) { ##for each tree
                               DatSpp <- SS.sum[SS.sum$Spp == currTrees[k],]
                               dat <- data.frame("Period" = rescale(as.numeric(DatSpp$FuturePeriod), 
                                                                    to = c(2000,2085)), 
                                                 "SIBEC" = DatSpp$MeanSI, "Suit" = DatSpp$MeanSuit)
                               dat$SIBEC <- dat$SIBEC/50 ##for mean annual increment
                               dat <- merge(dat, SuitProb, by = "Suit")
                               s <- approx(dat$Period, dat$SIBEC, n = 101) ##Smooth SI
                               p <- approx(dat$Period, dat$ProbDead, n = 101) ###Smooth Prob Dead
                               m <- approx(dat$Period, dat$NoMort, n = 101) ##Smooth No Mort
                               r <- approx(dat$Period, dat$RuinSeverity, n = 101)
                               
                               ###data frame of annual data
                               annualDat <- data.frame("Year" = seq(2000,2100,1), "Growth" = s[["y"]], 
                                                       "MeanDead" = p[["y"]], "NoMort" = m[["y"]], "RuinSeverity" = r[["y"]]) ##create working data
                               simRes <- SimGrowth_Volume(DF = annualDat)
                               out <- data.table(Spp = currTrees[k],Vol = simRes[1], nTree = simRes[2])
                               output <- rbind(output,out)
                             } 
                             data.table(it = c(i,i),
                                        Stat = c("Vol","nTree"),
                                        MaxReturn = c(sum(output$Vol*t(risk)),sum(output$nTree*t(risk))), 
                                        Return90 = c(sum(output$Vol*t(v90)),sum(output$nTree*t(v90))), 
                                        Sharpe = c(sum(output$Vol*t(safe)),sum(output$nTree*t(safe))))
                            }
                           
                           #simVolume <- simVolume[Stat == "nTree",]
                           simVolume <- melt(simVolume, id.vars = c("it","Stat"))
                           #simVolume <- simVolume[value < 8000,] ##something weird happens occasionally - need to sort this out
                           # ggplot(simVolume, aes(x = variable, y = value))+
                           #   geom_violin()+
                           #   labs(x = "Portfolio Choice", y = "Volume of Stand")
                           
                           simVolume$SiteNo <- SNum
                           simVolume
                         }else{
                           NULL
                         }
                         
                       }

allSitesSpp <- as.data.table(allSitesSpp)
allSitesSpp <- allSitesSpp[value < 5000,]
allSitesSpp[value < 0, value := 0]
simVol <- allSitesSpp[Stat == "Vol",]
simNumTree <- allSitesSpp[Stat == "nTree",]
simNumTree <- simNumTree[value < 100,]

ggplot(simVol, aes(x = variable, y = value))+
  geom_violin(draw_quantiles = 0.5)+
  labs(x = "Portfolio Choice", y = "Volume of Stand")+
  ggtitle(BGC)

ggplot(simNumTree, aes(x = variable, y = value))+
  geom_violin(draw_quantiles = 0.5)+
  labs(x = "Portfolio Choice", y = "# of Trees")+
  ggtitle(BGC)
    