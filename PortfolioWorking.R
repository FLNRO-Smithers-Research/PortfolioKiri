### Portfolio Working Script

cleanData <- function(SSPredAll,SIBEC,SuitTable,SNum,Trees,timePer,selectBGC){
  SSPred <- SSPredAll[SiteNo == SNum,] ###subset
  
  ##Merge SIBEC data
  SIBEC <- SIBEC[TreeSpp %in% Trees,]
  SSPred <- SSPred[SSPred$FuturePeriod %in% timePer,]
  SSPred <- SSPred[,c(6,3,4)]
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
  colnames(SSPred)[4] <- "Spp"
  
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

combineList <- function(...) {##Combine multiple dataframe in foreach loop
  mapply(FUN = rbind, ..., SIMPLIFY=FALSE)
}

edatopicSubset <- function(SSPredOrig, eda, pos = "Zonal"){
  if(pos == "Zonal"){
    SSPredFull <- SSPredOrig[grep("01",SSPredOrig$SSCurrent),]
  }else{
    edaSub <- eda[Edatopic == pos,]
    SSPredFull <- SSPredOrig[SSCurrent %in% edaSub$SS_NoSpace,]
  }
  SSPredFull$BGC_analysis <- gsub("/.*","", SSPredFull$SSCurrent)
  
  SSPredFull <- SSPredFull[,c("MergedBGC", "Source", "SS_NoSpace", "SSprob", "SSCurrent", 
                              "FuturePeriod", "SiteNo","BGC_analysis")]
  
  ##remove cases where not all timeperiods available
  SSPredFull <- as.data.table(SSPredFull)
  temp <- SSPredFull[,.(Num = length(unique(FuturePeriod))), by = c("SiteNo","SSCurrent")]
  temp <- temp[Num == 3,-c("Num")]
  SSPredFull <- SSPredFull[temp,on = c("SiteNo","SSCurrent")]
  SSPredFull[,SiteNo := as.numeric(SiteNo)]
}

####ProbDead- out of 100 trees, how many will die each year at each suitability. NoMort- Percent of time no mortality
#SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.35,1,1.8,4), "NoMort" = c(70,60,50,30), "RuinSeverity" = c(0.4,0.5,0.7,0.8)) 
SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.35,1,1.8,4), "NoMort" = c(70,60,50,30), "RuinSeverity" = c(0.7,0.7,0.7,0.7))                        




minAccept <- 0.01 ##min acceptable weight in porfolio - if lower, will remove and re-optimize

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

##adjust min and max weights (or remove species)
## see below examples
# boundDat[Spp == "Hw",maxWt := 0] ##remove Hw
# boundDat[Spp == "Pl",minWt := 0.1] ## set min weight for Pl to 0.1

boundDat2 <-  data.table::transpose(boundDat) %>% row_to_names (row_number = 1)

timePeriods = c(2000,2025,2055)
returnValue = 0.9
SSPredFull <- edatopicSubset(SSPredOrig,eda,pos = "Zonal")

treeList <- Trees
nSpp <- length(Trees)
Units <- unique(SSPredFull$BGC_analysis)
BGC = Units[1]##select BGC
  
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
                                 r <- approx(dat$Period, dat$RuinSeverity, n = 101)
                                 
                                 ###data frame of annual data
                                 annualDat <- data.frame("Year" = seq(2000,2100,1), "Growth" = s[["y"]], 
                                                         "MeanDead" = p[["y"]], "NoMort" = m[["y"]], "RuinSeverity" = r[["y"]]) ##create working data
                                 
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
                                 annualDat <- data.frame("Year" = seq(2000,2100,1))
                                 output <- data.frame("year" = annualDat$Year)
                                 for (k in 1:length(currTrees)){ ##for each tree
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
                                   Returns <- simGrowthCpp(DF = annualDat)
                                   output <- cbind(output, Returns)
                                 } ## for each 
                                 output <- as.data.table(output)
                                 setnames(output, c("Year",names(safe)))
                                 totVol <- output[nrow(output), !"Year"]
                                 data.table(it = i, MaxReturn = sum(totVol*risk), Return90 = sum(totVol*v90), Sharpe = sum(totVol*safe))
                                }
                               
                               simVolume <- melt(simVolume, id.vars = "it")
                               simVolume <- simVolume[value < 8000,] ##something weird happens occasionally - need to sort this out
                               # ggplot(simVolume, aes(x = variable, y = value))+
                               #   geom_violin()+
                               #   labs(x = "Portfolio Choice", y = "Volume of Stand")
                               
                               simVolume$SiteNo <- SNum
                               simVolume
                             }else{
                               NULL
                             }
                             
                           }

    allSitesSpp[value < 0, value := 0]
    ggplot(allSitesSpp, aes(x = variable, y = value))+
      geom_violin()+
      labs(x = "Portfolio Choice", y = "Volume of Stand")    
    
    