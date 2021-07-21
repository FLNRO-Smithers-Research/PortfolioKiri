### Portfolio Script to investigate range of simulation outcomes
###Kiri Daust, Nov 2020

require(foreach)
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
library(tictoc)

##Set drive with cloud data
if(dir.exists("C:/users/whmacken/Sync")){
  cloud_dir <- "C:/users/whmacken/Sync/Portfolio_Work/"
}else{
  cloud_dir <- "./PortfolioData/"
}

## source functions
sourceCpp("./CppFunctions/SimGrowth.cpp")
source("./PortfolioOpt_R.R")
#source("CCISS_Fns.R")

###load data - in package data
SuitTable <- fread("~/CommonTables/Feasibility_v12_7.csv") ##tree spp suitability
SuitTable[,Confirmed := NULL]
SuitTable[,SppVar := substr(SppVar,1,2)]
SuitTable <- unique(SuitTable)
setnames(SuitTable,c("BGC","SS_NoSpace","Spp","Suitability"))

SIBEC <- fread("InputsGit/PredSI_May2020.csv") 
SIBECnew <- fread("InputsGit/SI_to_add.csv")
SIBEC <- rbind(SIBEC, SIBECnew)
###import SI data (currently from BART)
SIBEC <- SIBEC[,c("SS_NoSpace","Spp","SIPred")] %>% set_colnames(c("SS_NoSpace","TreeSpp","MeanPlotSiteIndex"))

eda <- fread("~/CommonTables/Edatopic_v12_5.csv")
eda <- unique(eda[,.(BGC,SS_NoSpace, Edatopic)])

## get cciss data
siteids <- c(6487982,6484391,6484900,6485410,6485920)
library(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user = "postgres", password = "postgres", host = "138.197.168.220", 
                 port = 5432, dbname = "cciss")
library(ccissdev)
bgcDat <- dbGetCCISS(con,siteids,avg = F, scn = "ssp370")
sspreds <- edatopicOverlap(bgcDat,Edatope = E1) ##reuse from uData$sspreds
###rename and cleanup
SSPredOrig <- sspreds
SSPredOrig[,allOverlap := NULL]
setnames(SSPredOrig, old = c("BGC","SiteRef"), new = c("MergedBGC","SiteNo"))
SSPredOrig <- SSPredOrig[,.(MergedBGC,SS_NoSpace,SSratio,SSprob,SS.pred,FuturePeriod,SiteNo)]

##this can be package function
cleanData <- function(SSPredAll,SIBEC,SuitTable,SNum,Trees,timePer,selectBGC){
  SSPred <- SSPredAll[SiteNo == SNum,] ###subset
  
  ##Merge SIBEC data
  SIBEC <- SIBEC[TreeSpp %in% Trees,]
  SSPred <- SSPred[SSPred$FuturePeriod %in% timePer,]
  SSPred <- SSPred[,.(FuturePeriod, SS_NoSpace,SS.pred, SSprob)]
  SSPred <- SIBEC[SSPred, on = c(SS_NoSpace = "SS.pred")]
  setnames(SSPred, old = c("SS_NoSpace","i.SS_NoSpace"), new = c("SS.pred","SS_NoSpace"))
  
  ###Add rows for species with missing SI - mostly US units here
  add <- foreach(Year = unique(SSPred$FuturePeriod), .combine = rbind) %do%{
    byYear <- SSPred[SSPred$FuturePeriod == Year,]
    foreach(SS = unique(byYear$SS.pred), .combine = rbind) %do%{
      bySS <- byYear[byYear$SS.pred == SS,]
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
  SSPred[SuitTable, Suitability := i.Suitability, on = c(SS.pred = "SS_NoSpace", "Spp")]
  SSPred$Suitability[is.na(SSPred$Suitability)] <- 5
  
  temp <- SIBEC[SS_NoSpace == selectBGC,]
  if(nrow(temp) == 0){
    return(NULL)
  }

  ##Summarise data- average SI and Suit weighted by SSProb
  SS.sum <- SSPred[,.(MeanSI = sum(MeanPlotSiteIndex*(SSprob/sum(SSprob))),
                      MeanSuit = round(sum(Suitability*(SSprob/sum(SSprob))), digits = 0)),
                   by = .(Spp,FuturePeriod)]
  
  SS.sum <- SS.sum[FuturePeriod %in% timePer,]
  ###not sure what we were doing here?
  SS.sum[MeanSuit == 4, MeanSI := 5]
  SS.sum[MeanSuit == 5, MeanSI := 0]
  SS.sum[MeanSuit == 5, MeanSuit := 4]
  SS.sum <- unique(SS.sum)
  setorder(SS.sum,Spp,FuturePeriod)
  return(SS.sum)
}

##also package function
edatopicSubset <- function(SSPredOrig, eda, pos = "Zonal"){
  if(pos == "Zonal"){
    SSPredFull <- SSPredOrig[grep("01",SSPredOrig$SS_NoSpace),]
  }else{
    edaSub <- eda[Edatopic == pos,]
    SSPredFull <- SSPredOrig[SS_NoSpace %in% edaSub$SS.pred,]
  }
  SSPredFull[,BGC_analysis := gsub("/.*","", SS_NoSpace)]
  
  SSPredFull <- SSPredFull[,c("MergedBGC", "SS.pred", "SSprob", "SS_NoSpace", 
                              "FuturePeriod", "SiteNo","BGC_analysis")]
  return(SSPredFull)
  ##remove cases where not all timeperiods available
  # temp <- SSPredFull[,.(Num = length(unique(FuturePeriod))), by = c("SiteNo","SS_NoSpace")]
  # temp <- temp[Num == 6,-c("Num")]
  # SSPredFull <- SSPredFull[temp,on = c("SiteNo","SS_NoSpace")]
  # SSPredFull[,SiteNo := as.numeric(SiteNo)]
}


####ProbDead- out of 100 trees, how many will die each year at each suitability. NoMort- Percent of time no mortality
#SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.35,1,1.8,4), "NoMort" = c(70,60,50,30), "RuinSeverity" = c(0.4,0.5,0.7,0.8)) 
#SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.1,0.2,0.3, 1), "NoMort" = c(95,85,75,50), "RuinSeverity" = c(0.5,0.5,0.5,1))                        
#SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.1,0.2,0.3, 1), "NoMort" = c(95,85,75,50))#, "RuinSeverity" = c(1,2,3,4))
#SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.35,1,1.8,4)), "NoMort" = c(95,85,75,50) 
#ProbPest <- 1/100 ## annual probability of an outbreak for a species

#SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.1,0.2,0.3, 1), "NoMort" = c(95,85,75,50), "RuinSeverity" = c(0.3,0.35,0.4,.8))
SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.1,0.5,1,4), "NoMort" = c(95,85,75,50))
ProbPest <- 1/1000 ## annual probability of an outbreak for a species

minAccept <- 0.01 ##min acceptable weight in porfolio - if lower, will remove and re-optimize

cols <- fread("./InputsGit/PortfolioSppColours.csv") ##in package data
cols <- cols[HexColour != "",]
myPal <- cols$HexColour
names(myPal) <- cols$TreeCode
Trees <- c("Py","Fd","At","Pl","Sx","Bl","Cw","Hw","Pw","Ss","Lw","Ba","Hm","Dr","Mb")
#Trees <- unique(SIBEC$TreeSpp)
myColours <- data.table(TreeCode = Trees)
myColours <- cols[myColours, on = "TreeCode"]
myColours <- myColours[!is.na(HexColour),]
pal <- myColours$HexColour
names(pal) <- myColours$TreeCode
colScale <- scale_fill_manual(name = "variable", values = pal)
#Trees <- myColours$TreeCode

##set min and max weights for each species
boundDat <- data.table(Spp = Trees)
boundDat[,`:=`(minWt = 0, maxWt = 1)]

unique(SSPredOrig$FuturePeriod)

timePeriods = c(1961,1991,2021,2041,2061)
returnValue = 0.9

###subset for zonal or other
SSPredFull <- edatopicSubset(SSPredOrig,eda,pos = "Zonal")

treeList <- Trees
#treeList <- c("Py","Fd","At","Pl","Sx","Bl","Cw","Hw","Pw","Ss","Lw","Ba","Hm","Dr","Mb")
nSpp <- length(treeList)
Units <- unique(SSPredFull$BGC_analysis)

###SELECT BGC
BGC = Units[1]
  
SSPredBGC <- SSPredFull[BGC_analysis == BGC,-("BGC_analysis")]
SSList <- unique(SSPredBGC$SS_NoSpace)
selectBGC = SSList[1]
SSPredAll <- SSPredBGC[SSPredBGC$SS_NoSpace == selectBGC,]

SiteList <- unique(SSPredAll$SiteNo)
#SiteList <- rep(SiteList, each = round(15/length(SiteList)))
SSPredAll <- SSPredAll[SiteNo %in% SiteList & !is.na(SSprob),]

#################################
###model climate variability
##summarise by fp to get mean
# dat2 <- rawDat[ID2 == BGC,]
# dat2[,c("GCM","Scenario","FuturePeriod") := tstrsplit(Year, "_")]
# dat2 <- dat2[,.(GCM,FuturePeriod,Scenario,ID1,CMD,Tmin_sp,Tmin_sm,Tmax_sp,Tmax_sm)]
# dat2[,FuturePeriod := as.numeric(gsub(".gcm","",FuturePeriod))]
# dat2 <- dat2[,lapply(.SD,mean), by = .(ID1,FuturePeriod), .SDcols = -c("GCM","Scenario")]

##annual data to get variance - could get from summaries?
# annDat <- fread("PortPoints_Quesnel_1960-2019MSY.csv")
# ann2 <- annDat[ID2 == "SBSmw",.(Year,ID1,CMD,Tmin_sp,Tmax_sm)]
# library(fitdistrplus)
# ##fit normal distribution
# f1 <- fitdist(ann2$CMD, "norm")
# denscomp(f1)

# library(RPostgreSQL)
# drv <- dbDriver("PostgreSQL")
dbDisconnect(con)
con <- dbConnect(drv, user = "postgres", password = "postgres", host = "138.197.168.220", 
                 port = 5432, dbname = "bgc_climate_data") ##connect to climate summaries
#dbDisconnect(con)


climVarFut <- dbGetQuery(con, paste0("select bgc,period,stat,climvar,value from szsum_fut where bgc in ('"
                                  ,BGC,"') and period in ('2021-2040','2041-2060','2061-2080') 
                                  and stat = 'mean'
                                  and climvar in ('CMD','Tmin_sp','Tmax_sm') and scenario = 'ssp370'"))
climVarFut <- as.data.table(climVarFut)

climVarCurr <- dbGetQuery(con, paste0("select bgc,period,stat,climvar,value from szsum_curr where bgc in ('"
                                     ,BGC,"') and period in ('1991 - 2020') 
                                  and stat in ('st.dev.Ann','mean') 
                                  and climvar in ('CMD','Tmin_sp','Tmax_sm')")) %>% as.data.table()
climVarCurr[stat == 'st.dev.Ann',stat := "stdev"]
climVarSD <- climVarCurr[stat == "stdev",]
climVarCurr <- climVarCurr[stat != "stdev",]
climVar <- rbind(climVarCurr,climVarFut)
climVar[,period := as.numeric(substr(period,1,4))]

simulateClimate <- function(climVar){ ##package function
  climParams <- list()
  simResults <- data.table()
  
  for(cvar in c("CMD","Tmin_sp","Tmax_sm")){
    climSub <- climVar[climvar == cvar,.(value = mean(value)), by = .(period)]
    climSD <- climVarSD[climvar == cvar,value]
    ##table of means by period
    dat <- data.table(Year = c(2000,2025,2055,2085),Mean = climSub$value)
    s <- approx(dat$Year, dat$Mean, n = 101) ##smooth
    
    ##simulate using mean and variance
    res <- numeric()
    for(i in 1:101){
      res[i] <- rnorm(1,mean = s$y[i],sd = climSD)
    }
    temp <- data.table(Year = 2000:2100, Value = res)
    temp[,Var := cvar]
    simResults <- rbind(simResults,temp,fill = T)
  }
  simResults <- dcast(simResults,Year ~ Var, value.var = "Value")
  return(simResults)
}
##tmin

##find limits for each species

sppLimits <- foreach(spp = Trees, .combine = rbind) %do% {##package function
  temp <- SuitTable[Suitability == 1 & Spp == spp,] ##what units is Fd 1?
  sppUnits <- unique(temp$BGC)
  
  climSum <- dbGetQuery(con, paste0("select bgc,period,stat,climvar,value from szsum_curr where bgc in ('"
                                    ,paste(sppUnits,collapse = "','"),"') and period = '1991 - 2020' 
                                    and climvar in ('CMD','Tmin_sp','Tmax_sm')"))
  climSum <- as.data.table(climSum)
  climSum2 <- climSum[,.(Min = min(value),Max = max(value)),
                      by = .(stat,climvar)]
  #climSum2 <- climSum2[var == "mean",]
  # climSum3 <- data.table(CMDMin = climSum2[stat == "10%" & climvar == "CMD",Min],
  #                        CMDMax = climSum2[stat == "90%" & climvar == "CMD",Max],
  #                        Tlow = climSum2[stat == "10%" & climvar == "Tmin_sp",Min],
  #                        Thigh = climSum2[stat == "90%" & climvar == "Tmax_sm",Max])
  climSum3 <- data.table(CMDMin = climSum2[stat == "mean" & climvar == "CMD",Min],
                         CMDMax = climSum2[stat == "mean" & climvar == "CMD",Max],
                         Tlow = climSum2[stat == "mean" & climvar == "Tmin_sp",Min],
                         Thigh = climSum2[stat == "mean" & climvar == "Tmax_sm",Max])
  climSum3[,Spp := spp]
  climSum3
}


# ggcmd <- ggplot(data = simResults[Var == "CMD",], aes(x = Year, y = Value))+
#   geom_line()+
#   geom_hline(yintercept = climSum2$CMDMax, col = "red")+
#   geom_hline(yintercept = climSum2$CMDMin, col = "blue")+
#   ggtitle("CMD")
# 
# ggtmin <- ggplot(data = simResults[Var == "Tmin_sp",], aes(x = Year, y = Value))+
#   geom_line()+
#   geom_hline(yintercept = climSum2$Tlow, col = "blue")+
#   ggtitle("Tmin_sp")
# 
# ggtmax <- ggplot(data = simResults[Var == "Tmax_sm",], aes(x = Year, y = Value))+
#   geom_line()+
#   geom_hline(yintercept = climSum2$Thigh, col = "red")+
#   ggtitle("Tmax_sm")
# 
# library(gridExtra)
# grid.arrange(ggcmd,ggtmin,ggtmax, ncol = 3)

##non-package function
SL <- SiteList
SL <- rep(SL, each = 5)
allSitesSpp <- foreach(SNum = SL, .combine = rbind, 
                       .packages = c("Rcpp","magrittr","reticulate",
                                     "data.table","foreach","ggplot2","ggthemes","scales"), 
                       .noexport = c("gs2gw", "simGrowthCBST","simGrowthCpp"),
                       .export = c("cleanData","loopCombine")) %do% {
                         
                         cat("Optimising site",SNum,"...\n")
                         ##simulate climate
                         simResults <- simulateClimate(climVar)
                         SS.sum <- cleanData(SSPredAll,SIBEC,SuitTable,SNum, Trees, 
                                             timePer = timePeriods,selectBGC = selectBGC)
                         if(any(is.na(SS.sum$MeanSuit))){
                           warning("Missing Suitability in unit ",
                                   BGC,", sitenumber ",SNum," for ",
                                   SS.sum$Spp[is.na(SS.sum$MeanSuit)], ": They will be filled with suit = 4")
                           SS.sum$MeanSuit[is.na(SS.sum$MeanSuit)] <- 4
                         }
                         SS.sum[,FuturePeriod := as.numeric(FuturePeriod)]
                         if(length(timePeriods) == 1){
                           temp <- SS.sum
                           temp$FuturePeriod <- SS.sum$FuturePeriod[1]+85
                           SS.sum <- rbind(SS.sum, temp)
                         }
                         
                         if(!is.null(SS.sum)){
                           output <- data.table("year" = seq(2000,2100,1))
                           
                           for (k in 1:nSpp){ ##for each tree
                             DatSpp <- SS.sum[Spp == treeList[k],]
                             dat <- data.table("Period" = rescale(as.numeric(DatSpp$FuturePeriod), 
                                                                  to = c(2000,2085)), 
                                               "SIBEC" = DatSpp$MeanSI/50, "Suit" = DatSpp$MeanSuit)
                             
                             dat <- merge(dat, SuitProb, by = "Suit")
                             s <- approx(dat$Period, dat$SIBEC, n = 101) ##Smooth SI
                             p <- approx(dat$Period, dat$ProbDead, n = 101) ###Smooth Prob Dead
                             m <- approx(dat$Period, dat$NoMort, n = 101) ##Smooth No Mort
                             r <- approx(dat$Period, dat$Suit, n = 101)
                             ##r <- approx(dat$Period, dat$RuinSeverity, n = 101)
                             
                             ###data frame of annual data
                             annualDat <- data.table("Growth" = s[["y"]], "MeanDead" = p[["y"]], "NoMort" = m[["y"]], "Suit" = r[["y"]]) ##create working data
                             annualDat <- cbind(simResults,annualDat)
                             limits <- sppLimits[Spp == treeList[k],]
                             Returns <- SimGrowth_v2(DF = annualDat,ProbPest = 0.005,
                                                     cmdMin = limits[[1]],cmdMax = limits[[2]],
                                                     tempMin = limits[[3]],tempMax = limits[[4]],climLoss = 0.08)
                             tmpR <- c(0,Returns)
                             assets <- Returns - tmpR[-length(tmpR)]
                             temp <- data.frame(Spp = rep(treeList[k],101), 
                                                Year = 1:101, Returns = Returns)
                             output <- cbind(output, assets)
                           } ## for each tree species
                           
                           colnames(output) <- c("Year", treeList)
                           
                           ####Portfolio#######################################
                           returns <- output
                           returns[,Year := NULL]
                           ###only include species with mean return > 1 in portfolio
                           use <- colnames(returns)[colMeans(returns) > 1] ###should probably be higher
                           returns <- returns[,..use]
                           sigma2 <- cor(returns) ###to create cov mat from returns
                           
                           ef <- run_opt(returns, sigma2, boundDat,minAccept) 
                           setnames(ef,old = c("frontier_sd","return","sharpe"),
                                    new = c("Sd","RealRet","Sharpe"))
                           ef[,Return := 1:20]
                           
                           eff_front2 <- ef
                           eff_front2[,RealRet := RealRet/max(RealRet)]
                           eff_front2[,SiteNo := SNum]
                           melt(eff_front2,id.vars = c("SiteNo", "Return"),variable.name = "Spp")
                         }else{NULL}
                       }
##preprocess for plotting
efAll <- allSitesSpp
efAll <- dcast(efAll,Return ~ Spp, fun.aggregate = function(x){sum(x)/(length(SL))})
efAll <- na.omit(efAll)
#efAll$RealRet <- efAll$RealRet/max(efAll$RealRet) ##standardise return
RetCurve <- approx(efAll$RealRet,efAll$Sd,xout = returnValue)
ret90 <- RetCurve$y
maxSharpe <- efAll[Sharpe == max(Sharpe),!c("Return","Sharpe")]
maxSPos <- maxSharpe$Sd
maxSharpe <- t(maxSharpe) %>% as.data.frame() %>% 
  mutate(Spp = rownames(.)) %>% set_colnames(c("Sharpe_Opt","Spp"))
ret90Props <- efAll[which.min(abs(RealRet - returnValue)),-c("Return","Sharpe")]
ret90Props <- t(ret90Props) %>% as.data.frame() %>% 
  mutate(Spp = rownames(.)) %>% set_colnames(c("Set_Return","Spp"))
maxSharpe$SSCurrent <- selectBGC
maxSharpe$Unit <- BGC
maxSharpe$Set_Return <- ret90Props$Set_Return
maxSharpe$Set_Return[maxSharpe$Spp == "Sd"] <- ret90
efAll <- efAll[,-c("Return","Sharpe")]
efAll <- melt(efAll, id.vars = "Sd")
efAll$Unit <- BGC

ef_plot <- function(efAll,intDat){
  # efAll <- outAll$GraphDat
  # intDat <- outAll$MaxS
  efAll$variable <- factor(efAll$variable, levels = sort(unique(as.character(efAll$variable))))
  ggplot(efAll[efAll$variable != "RealRet",],aes(x = Sd, y = value,group = variable))+
    geom_area(aes(fill = variable), size = 0.00001, col = "grey50", stat = "identity")+
    colScale +
    geom_vline(data = intDat[intDat$Spp == "Sd",], aes(xintercept = Sharpe_Opt,colour = "blue"), 
               linetype = "twodash", size = .75)+
    geom_vline(data = intDat[intDat$Spp == "Sd",], aes(xintercept = Set_Return,colour = "grey52"),
               linetype = "dashed", size = .75)+
    geom_line(data = efAll[efAll$variable == "RealRet",], 
              aes(x = Sd, y = value,colour = "black"),linetype = "F1",size = .75)+
    scale_colour_identity(name = "", guide = 'legend', labels = c("Return","MaxSharpe","90%"))+
    scale_x_reverse() +
    xlab("Max Return --> Minimized Risk")+
    ylab("Portfolio Ratio")+
    guides(fill=guide_legend("Species"))+
    theme_few()+
    facet_wrap(.~Unit, scales = "free_x")
}
ef_plot(efAll,maxSharpe)

##table
temp <- as.data.table(maxSharpe)
temp <- temp[!Spp %chin% c("RealRet","Sd"),.(Spp,Sharpe_Opt,SetRet)]
knitr::kable(temp, digits = 2)
###################################################################################################




##volume simulations
SL <- SiteList
allSitesSpp <- foreach(SNum = SL, .combine = rbind, 
                       .packages = c("Rcpp","magrittr","reticulate",
                                     "data.table","foreach","ggplot2","ggthemes","scales"), 
                       .noexport = c("gs2gw", "simGrowthCBST","simGrowthCpp"),
                       .export = c("cleanData","loopCombine")) %do% {
                         reticulate::source_python("./PythonFns/PortfolioOptimisation.py") 
                         
                         cat("Optimising site",SNum,"...\n")
                         ##simulate climate
                         simResults <- simulateClimate(climVar)
                         SS.sum <- cleanData(SSPredAll,SIBEC,SuitTable,SNum, Trees, 
                                             timePer = timePeriods,selectBGC = selectBGC)
                         if(any(is.na(SS.sum$MeanSuit))){
                           warning("Missing Suitability in unit ",
                                   BGC,", sitenumber ",SNum," for ",
                                   SS.sum$Spp[is.na(SS.sum$MeanSuit)], ": They will be filled with suit = 4")
                           SS.sum$MeanSuit[is.na(SS.sum$MeanSuit)] <- 4
                         }
                         SS.sum[,FuturePeriod := as.numeric(FuturePeriod)]
                         if(length(timePeriods) == 1){
                           temp <- SS.sum
                           temp$FuturePeriod <- SS.sum$FuturePeriod[1]+85
                           SS.sum <- rbind(SS.sum, temp)
                         }
                         
                         if(!is.null(SS.sum)){
                           output <- data.table("year" = seq(2000,2100,1))
                           
                           for (k in 1:nSpp){ ##for each tree
                             DatSpp <- SS.sum[Spp == treeList[k],]
                             dat <- data.table("Period" = rescale(as.numeric(DatSpp$FuturePeriod), 
                                                                  to = c(2000,2085)), 
                                               "SIBEC" = DatSpp$MeanSI/50, "Suit" = DatSpp$MeanSuit)
 
                             dat <- merge(dat, SuitProb, by = "Suit")
                             s <- approx(dat$Period, dat$SIBEC, n = 101) ##Smooth SI
                             p <- approx(dat$Period, dat$ProbDead, n = 101) ###Smooth Prob Dead
                             m <- approx(dat$Period, dat$NoMort, n = 101) ##Smooth No Mort
                             r <- approx(dat$Period, dat$Suit, n = 101)
                             ##r <- approx(dat$Period, dat$RuinSeverity, n = 101)
                             
                             ###data frame of annual data
                             annualDat <- data.table("Growth" = s[["y"]], "MeanDead" = p[["y"]], "NoMort" = m[["y"]], "Suit" = r[["y"]]) ##create working data
                             annualDat <- cbind(simResults,annualDat)
                             limits <- sppLimits[Spp == treeList[k],]
                             Returns <- SimGrowth_v2(DF = annualDat,ProbPest = 0.005,
                                                     cmdMin = limits[[1]],cmdMax = limits[[2]],
                                                     tempMin = limits[[3]],tempMax = limits[[4]])
                             tmpR <- c(0,Returns)
                             assets <- Returns - tmpR[-length(tmpR)]
                             temp <- data.frame(Spp = rep(treeList[k],101), 
                                                Year = 1:101, Returns = Returns)
                             output <- cbind(output, assets)
                           } ## for each tree species
                           
                           colnames(output) <- c("Year", treeList)
                           
                           ####Portfolio#######################################
                           returns <- output
                           returns[,Year := NULL]
                           ###only include species with mean return > 1 in portfolio
                           use <- colnames(returns)[colMeans(returns) > 1] ###should probably be higher
                           returns <- returns[,..use]
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
                             climSim <- simulateClimate(climVar)
                             output <- data.table(Spp = character(),Vol = numeric())
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
                               r <- approx(dat$Period, dat$Suit, n = 101)
                               
                               ###data frame of annual data
                               annualDat <- data.table("Growth" = s[["y"]], "MeanDead" = p[["y"]], "NoMort" = m[["y"]], "Suit" = r[["y"]]) ##create working data
                               annualDat <- cbind(climSim,annualDat)
                               limits <- sppLimits[Spp == treeList[k],]
                               simRes <- SimGrowth_v2(DF = annualDat,ProbPest = 0.005,
                                                       cmdMin = limits[[1]],cmdMax = limits[[2]],
                                                       tempMin = limits[[3]],tempMax = limits[[4]])
                               out <- data.table(Spp = currTrees[k],Vol = simRes[length(simRes)])
                               output <- rbind(output,out)
                             } 
                             data.table(it = i,
                                        Stat = c("Vol"),
                                        MaxReturn = c(sum(output$Vol*t(risk))), 
                                        Return90 = c(sum(output$Vol*t(v90))), 
                                        Sharpe = c(sum(output$Vol*t(safe))))
                            }
                           
                           #simVolume <- simVolume[Stat == "nTree",]
                           simVolume <- melt(simVolume, id.vars = c("it","Stat"))
                           # simVolume <- simVolume[value < 8000,] ##something weird happens occasionally - need to sort this out
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

ggplot(simVol, aes(x = variable, y = value))+
  geom_violin(draw_quantiles = 0.5,scale = "width")+
  labs(x = "Portfolio Choice", y = "Volume of Stand")+
  ggtitle(BGC)

ggplot(simNumTree, aes(x = variable, y = value))+
  geom_violin(draw_quantiles = 0.5, scale = "width")+
  labs(x = "Portfolio Choice", y = "# of Trees")+
  ggtitle(BGC)


##investigate multiple GCM runs
tic()
tempclm <- climatebc_mult(inputFile = "ClimBCTest.csv",vip = 1,ysm = "YSM",period = "ACCESS1-0_rcp85_2025.gcm")
toc()

fname <- "./TimeSeries/PemSamplePts_CanESM2_RCP45_r"
dat <- foreach(it = 1:5, .combine = rbind) %do%{
  temp <- fread(paste0(fname,it,"1i1p1_2020-2100MSY.csv"))
  temp[,ModRun := it]
  temp
}

temp <- data.table(Year = 2020:2100, Chunk = c(rep(1:8, each = 10),8))
# dat3[temp, Cunk := i.Chunk, on = "Year"]
# statMean <- dat3[,lapply(.SD, mean),by = .(Cunk,ModRun), .SDcols = -"Year"]
# statMean[,Stat := "Mean"]
# statVar <- dat3[,lapply(.SD, var),by = .(Cunk,ModRun), .SDcols = -"Year"]
# statVar[,Stat := "Var"]
# statByRun <- rbind(statMean,statVar)
# fwrite(statByRun,"TenYearStatsByRun.csv")
# 
# statMean2 <- dat3[,lapply(.SD, mean),by = .(Cunk), .SDcols = -c("Year","ModRun")]
# statMean2[,Stat := "Mean"]
# statVar2 <- dat3[,lapply(.SD, var),by = .(Cunk), .SDcols = -c("Year","ModRun")]
# statVar2[,Stat := "Var"]
# statCombined2 <- rbind(statMean2,statVar2)
# fwrite(statCombined2,"TenYearStats.csv")
# 
# dat3[,Cunk := as.factor(Cunk)]
# ggplot(dat3, aes(x = Cunk, y = CMD, col = ModRun))+
#   geom_boxplot()

dat3 <- dat[,.(Year,ModRun,ID2,CMD,DD5,Tmin_wt,MCMT,EXT,MSP)]
dat3[temp, Chunk := i.Chunk, on = "Year"]
dat3[,ModRun := as.factor(ModRun)]
dat3[,Chunk := as.factor(Chunk)]
datLong <- melt(dat3, id.vars = c("Chunk","Year","ModRun","ID2"), value.name = "Value",variable.name = "ClimVar")

ggplot(datLong, aes(x = Chunk, y = Value, col = ModRun))+
  geom_boxplot()+
  facet_wrap(.~ClimVar, scales = "free_y")

# ###Create current data
# current <- temp %>% 
#   merge(SuitTable, by.x = c("TreeSpp","SS.pred"), by.y = c("Spp","SS.pred"), all.x = TRUE) %>%
#   unique()
# 
# ###check that there aren't errors in the table
# temp <- aggregate(SS.pred ~ TreeSpp, current, FUN = length)
# if(any(temp$SS.pred > 1)){
#   stop("There are partial duplicates in the suitablity table. Please fix them. :)")
# }
# 
# current <- data.frame(Spp = current$TreeSpp, FuturePeriod = 2000, 
#                       MeanSI = current$MeanPlotSiteIndex, MeanSuit = current$Suitability)
# 
# missing <- Trees[!Trees %in% current$Spp]
# if(length(missing) > 0){
#   new <- current[rep(1,length(missing)),]
#   new$Spp <- missing
#   new$MeanSI <- 10
#   new$MeanSuit <- 5
#   current <- rbind(current, new)
# }