##This script uses a Monte Carlo approach to model tree growth wiht SI numbers and predicted suitability 
##It then uses the Markowitz portfolio method to calculate the optimal mix of species.
##Kiri Daust, June 2018


##.libPaths("E:/R packages351")

require(stats)
require(rgl)
require(tcltk)
require(bindr)
require(lattice)
require(dplyr)
require(ggplot2)
require(MonteCarlo)
require(openxlsx)
require(MASS)
require(xts)
require(PortfolioAnalytics)
require(ROI.plugin.quadprog)
require(ROI)
require(MASS)
require(pcaPP)
require(tseries)
require(magrittr)
require(foreach)
require(reshape2)
require(doParallel)

rm(list=ls())
wd=tk_choose.dir(); setwd(wd)
setwd("C:/Users/Kiri Daust/Desktop/PortfolioKiri")

###OPTIONAL: Set up to run loops in parallel###
require(doParallel)
set.seed(123321)
coreNum <- as.numeric(detectCores()-2)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)
clusterEvalQ(coreNo, .libPaths("E:/R packages351"))

repeat.before = function(x) {   # repeats the last non NA value. Keeps leading NA
  ind = which(!is.na(x))      # get positions of nonmissing values
  if(is.na(x[1]))             # if it begins with a missing, add the 
    ind = c(1,ind)        # first position to the indices
  rep(x[ind], times = diff(   # repeat the values at these indices
    c(ind, length(x) + 1) )) # diffing the indices + length yields how often 
}

sigma <- read.csv("CovarianceMatrix_Full.csv")
rownames(sigma) <- sigma[,1]
sigma <- sigma[,-1]
Trees <- c("Bl","Cw","Fd","Hw","Lw","Pl","Py","Sx") ##set species to use in portfolio
nSpp <- length(Trees)
treeList <- Trees
sigma <- sigma[Trees,Trees]
sigma <- as.matrix(sigma)
momentargs <- list()
momentargs$sigma <- sigma ##set up to use in PortfolioAnalystics 

nsuitF <- file.choose() ##Import suitability table
SuitTable <- read.csv(nsuitF, stringsAsFactors = FALSE)
SuitTable <- unique(SuitTable)

colnames(SuitTable)[2] <- "SS_NoSpace"

SIBEC <- read.csv("PredSIforPort_Sept20.csv", stringsAsFactors = FALSE)
SIBEC <- SIBEC[,-5]
colnames(SIBEC)[c(1,3,5)] <- c("SS_NoSpace", "MeanPlotSiteIndex","TreeSpp")


SSPredAll <- read.csv(file.choose(), stringsAsFactors = FALSE) ##Import SS predictions from CCISS tool: must have columns MergedBGC, Source, SS_NoSpace, SSprob, SSCurrent, FuturePeriod, SiteNo

selectBGC <- select.list(choices = sort(unique(SSPredAll$SSCurrent)), graphics = TRUE) ###Select BGC to run for
SSPredAll <- SSPredAll[SSPredAll$SSCurrent == selectBGC,]

#####Randomly select 100 sites############
sites <- as.numeric(as.character(unique(SSPredAll$SiteNo)))
SiteList <- sample(sites, 100, replace = FALSE)
SSPredAll <- SSPredAll[SSPredAll$SiteNo %in% SiteList,]

SSPredAll <- SSPredAll[!is.na(SSPredAll$SSprob),]
########################

combineList <- function(...) {##Combine multiple dataframe in foreach loop
  mapply(FUN = rbind, ..., SIMPLIFY=FALSE)
}

###foreach site
allSitesSpp <- foreach(SNum = unique(SSPredAll$SiteNo), .combine = combineList, .packages = c("foreach","reshape2","dplyr","magrittr","PortfolioAnalytics")) %dopar% {
    SSPred <- SSPredAll[SSPredAll$SiteNo == SNum,]
    EF.out.all <- rep(1:25, each = nSpp)
    EF.ret.all <- 1:25
    
    ##Merge SIBEC data
    SIBEC <- SIBEC[SIBEC$TreeSpp %in% Trees,]
    SSPred <- SSPred[,c(6,3,4)]
    SSPred <- merge(SSPred, SIBEC[,c(1,3,5)], by = "SS_NoSpace", all.x = TRUE)
    
    ###Add rows for species with missing SI
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
    if(nrow(add) > 0){
      add$MeanPlotSiteIndex <- 10 ##Set missing SI to 10
      SSPred <- rbind(SSPred, add)
    }
    SSPred <- SSPred[!is.na(SSPred$TreeSpp),]
    colnames(SSPred)[5] <- "Spp"
    
    ##Add suitability
    SSPred <- merge(SSPred, SuitTable, by = c("SS_NoSpace","Spp"), all.x = TRUE)
    SSPred$Suitability[is.na(SSPred$Suitability)] <- 5
    
    ###Create current data
    current <- SIBEC[SIBEC$SS_NoSpace == selectBGC,c(1,3,5)]
    current <- merge(current, SuitTable, by.x = c("TreeSpp","SS_NoSpace"), by.y = c("Spp","SS_NoSpace"), all.x = TRUE)
    current <- unique(current)
    
    ###check that there aren't errors in the table
    temp <- aggregate(SS_NoSpace ~ TreeSpp, current, FUN = length)
    if(any(temp$SS_NoSpace > 1)){
      stop("There are partial duplicates in the suitablity table. Please fix them. :)")
    }
    
    current$Suitability[current$TreeSpp == "Fd"] <- 3
    current[is.na(current)] <- 5
    current[current == 0] <- 10
    current <- data.frame(Spp = current$TreeSpp, FuturePeriod = 2000, MeanSI = current$MeanPlotSiteIndex, MeanSuit = current$Suitability)
    
    missing <- Trees[!Trees %in% current$Spp]
    if(length(missing) > 0){
      new <- current[rep(1,length(missing)),]
      new$Spp <- missing
      new$MeanSI <- 10
      new$MeanSuit <- 5
      current <- rbind(current, new)
    }
    
    ##Summarise data- average SI and Suit weighted by SSProb
    SS.sum <- SSPred %>%
      group_by(Spp, FuturePeriod) %>%
      summarise(MeanSI = sum(MeanPlotSiteIndex*(SSprob/sum(SSprob))), MeanSuit = round(sum(Suitability*(SSprob/sum(SSprob))), digits = 0))
    
    SS.sum <- as.data.frame(SS.sum)
    SS.sum <- rbind(SS.sum, current)
    SS.sum$MeanSuit[SS.sum$MeanSuit == 5] <- 4
    SS.sum <- unique(SS.sum)
    SS.sum <- SS.sum[order(SS.sum$Spp,SS.sum$FuturePeriod),]
    
    cols <- rainbow(length(treeList))
    annualDat <- data.frame("Year" = seq(2000,2100,1))

    portOutput <- data.frame("Species" = treeList)###set up plot and output
    plot(0,0, type = "n", xlim = c(1,100), ylim = c(0,3000), xlab = "Year", ylab = "Volume")###plot
    
    for (w in 1:50){ ##number of iterations
      output <- data.frame("year" = annualDat$Year)
      
      for (k in 1:nSpp){ ##for each tree
        DatSpp <- SS.sum[SS.sum$Spp == treeList[k],]
        SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.5,0.5,1,4), 
                               "NoMort" = c(70,60,50,30)) ####ProbDead- out of 100 trees, how many will die each year at each suitability. NoMort- Percent of time no mortality
        dat <- data.frame("Period" = c(2000,2025,2055,2085), "SIBEC" = DatSpp$MeanSI, "Suit" = DatSpp$MeanSuit)
        dat$SIBEC <- dat$SIBEC/50 ##for mean annual increment
        dat <- merge(dat, SuitProb, by = "Suit")
        s <- approx(dat$Period, dat$SIBEC, n = 101) ##Smooth SI
        p <- approx(dat$Period, dat$ProbDead, n = 101) ###Smooth Prob Dead
        m <- approx(dat$Period, dat$NoMort, n = 101) ##Smooth No Mort
        
        ###data frame of annual data
        annualDat <- data.frame("Year" = seq(2000,2100,1), "Growth" = s[["y"]], "MeanDead" = p[["y"]], "NoMort" = m[["y"]]) ##create working data
        
        Returns <- vector(mode = "numeric", length = 101)
        modData <- data.frame("V1"= Returns)
        #probCheck <- vector(mode = "numeric", length = 101)
        
        ##growth simulation
          Returns <- vector(mode = "numeric", length = 86)
          nTrees <- 100 ##number of tree to start
          for (i in 1:101){ ##for each year
            height <- sum(annualDat$Growth[1:i]) ##total height
            Returns[i] <- nTrees*height ##volume
            if(runif(1, min = 0, max = 100) > annualDat$NoMort[i]){##will there be mortality?
              prevTrees <- nTrees
              percentDead <- rgamma(1, shape = 1, scale = annualDat$MeanDead[i])###what percent of trees will die based on gamma distribution?
              numDead <- (percentDead/100)*prevTrees###number of dead trees
              nTrees <- prevTrees - numDead ##update number of trees
            }
          } ##for each year
          modData <- cbind(modData, Returns)
        
        
        modData <- modData[,-1]
        lines(modData, type = 'l', col = cols[k])##plot
        
        assets <- vector("numeric", 101)
        assets[1] <- modData[1]
        for (z in 1:100){ ##convert from cumulative volume to change by year
          assets[z+1] <- modData[z+1] - modData[z]
        }
        output <- cbind(output, assets)
      } ## for each tree species
      
      legend("topleft", treeList, col = cols, pch = 15)
      
      colnames(output) <- c("Year", treeList)
      output$Year <- paste(output$Year,"-01-01", sep = "")
      
      ####Portfolio#######################################
      returns <- output
      rownames(returns) <- returns[,1]
      returns <- returns[,-1]
      returnsTS <- as.xts(returns) ###convert to time series object
      
      ##set portfolio constraints
      init.portfolio <- portfolio.spec(assets = colnames(returnsTS))
      init.portfolio <- add.constraint(portfolio = init.portfolio, type = "weight_sum", min_sum = 0.99, max_sum = 1.01) ###weights should add to about 1
      init.portfolio <- add.constraint(portfolio = init.portfolio, type = "box", min = rep(0, nSpp), max = rep(1, nSpp)) ##set min and max weight for each species
      
      ###loop optimisation with increasing risk aversion to obtain efficient frontier
      for(q in seq(from = 2, to = 50, by = 2)){
        qu <- add.objective(portfolio=init.portfolio, type="return", name="mean")
        qu <- add.objective(portfolio=qu, type="risk", name="var", risk_aversion = q)
        
        minSD.opt <- optimize.portfolio(R = returnsTS, portfolio = qu, optimize_method = "ROI", trace = TRUE, momentargs = momentargs)
        if(q == 2){
          EF.out <-  as.data.frame(minSD.opt[["weights"]])
          EF.ret <- as.data.frame(minSD.opt[["objective_measures"]][["mean"]])
        }else{
          temp <- as.data.frame(minSD.opt[["weights"]])
          EF.out <- rbind(EF.out, temp)
          EF.ret <- rbind(EF.ret, minSD.opt[["objective_measures"]][["mean"]])
        }
        
      }
      
      ##set portfolio objectives (currently set to maximise return given certain risk aversion- Quadratic Utility)
      qu <- add.objective(portfolio=init.portfolio, type="return", name="mean")
      qu <- add.objective(portfolio=qu, type="risk", name="var", risk_aversion = 15) ###set risk_aversion
      
      ##optimise
      minSD.opt <- optimize.portfolio(R = returnsTS, portfolio = qu, optimize_method = "ROI", trace = TRUE, momentargs = momentargs)
      test <- as.data.frame(minSD.opt[["weights"]]) ##save weights
      colnames(test) <- "Run"
      portOutput <- cbind(portOutput, test$Run)
      
      EF.out.all <- cbind(EF.out.all, EF.out)
      EF.ret.all <- cbind(EF.ret.all, EF.ret)
      
    }  
    
    ###calculate average of efficient portfolio weights and their returns for frontier chart
    EF.ret.all$Mean <- apply(EF.ret.all[,-1],1,mean) ###mean of returns
    EF.ret.all$EF.ret.all <- 2*EF.ret.all$EF.ret.all
    EF.ret.all$Mean <- EF.ret.all$Mean/max(EF.ret.all$Mean) ##scale return out of 1
    EF.ret.all <- EF.ret.all[,c(1,52)]
    colnames(EF.ret.all) <- c("Risk","MeanRet")
    EF.out.all$Mean <- apply(EF.out.all[,-1],1,mean) ### mean of weights
    EF.sum <- data.frame(Spp = rep(treeList,25), 
                         Risk = EF.out.all$EF.out.all*2, Weight = EF.out.all$Mean, Site = SNum)
    
    ###weights for violin plots
    rownames(portOutput) <- portOutput$Species
    portOutput <- portOutput[,-1]
    portOut <- t(portOutput)
    portOut <- as.data.frame(portOut)
    portOut$Site <- SNum
    
    ###Combine into list
    outList <- list(Frontier = EF.sum, Weights = portOut, Return = EF.ret.all)
    outList
}

EF.sum <- allSitesSpp$Frontier
EF.sum <- aggregate(Weight ~ Spp + Risk, EF.sum, FUN = mean)
EF.ret.all <- allSitesSpp$Return
EF.ret.all <- aggregate(MeanRet ~ Risk, EF.ret.all, FUN = mean)
myColours <- c("red","pink", "orange","yellow","green","blue","magenta","darkgoldenrod")
names(myColours) <- levels(EF.sum$Spp)
colScale <- scale_fill_manual(name = "Spp", values = myColours)

SpeciesPlot <- ggplot(EF.sum)+
  geom_bar(aes(x = Risk, y = Weight, fill = Spp), stat = "identity")+
  colScale+
  geom_line(data = EF.ret.all, aes(x = Risk, y = MeanRet), size = 2)+
  geom_vline(xintercept = 12)+
  labs(x = "Risk (0 = High Risk, 50 = Low Risk)")+
  ggtitle(("Species Portfolio"))


####CBST Portfolio#############
gs2gw <- function(x){###function to convert gs to gw
  5e-25*exp(x*55.9)
}
gs2gwVec <- Vectorize(gs2gw)

###import genetic matrix
setwd("C:/Users/Kiri Daust/Desktop/PortfolioKiri/CBSTPortfolio")
fullMat <- read.csv("Sx genetic suitability matrix - no assisted migration.csv")
rownames(fullMat) <- fullMat$BECvar
fullMat <- fullMat[,-1]
fullMatSx <- fullMat
fullMat <- read.csv("Pl genetic suitability matrix - no assisted migration.csv")
rownames(fullMat) <- fullMat$BECvar
fullMat <- fullMat[,-1]
fullMatPl <- fullMat

##import genetic list
SuitTable <- read.csv("Sx genetic suitability list - no assisted migration.csv")
colnames(SuitTable)[1:2] <- c("Site","Seed")
SuitTableSx <- SuitTable
SuitTable <- read.csv("Pl genetic suitability list - no assisted migration.csv")
colnames(SuitTable)[1:2] <- c("Site","Seed")
SuitTablePl <- SuitTable

SuitTable <- merge(SuitTablePl, SuitTableSx, by = c("Site","Seed"), all = TRUE)
colnames(SuitTable)[3:4] <- c("Spp1","Spp2")

##import BGC prediction by model
SSPredAll <- read.csv(file.choose())##import BGC predictions from CCISS script
BGC <- gsub("/.*","",selectBGC)
SSPredAll <- SSPredAll[SSPredAll$BGC == BGC,]
SSPredAll <- SSPredAll[SSPredAll$FuturePeriod == 2025,]
SSPredAll <- SSPredAll[SSPredAll$SiteNo %in% SiteList,]
SeedList <- as.character(unique(SuitTable$Seed))
outRaw <- data.frame(Seed = SeedList, Number = 0)


simulateGrowth <- function(x){
  if(any(x < 0) || max(x) < 0.4){return(NULL)}
  
  annualDat <- data.frame("Year" = seq(2000,2040,1))
  
  portOutput <- data.frame("Seed" = SeedList)###set up plot and output
  s <- spline(c(2000,2025), x, n = 41)
  growthRate <- s[["y"]]
  pDead <- 1 - s[["y"]]
  pDead <- rescale(pDead, to = c(0.01,0.1), from = c(0,1))
  
  nTrees <- 100 ##number of tree to start
  Returns <- numeric(length = 41)
  
  ###Simulate Growth
  for (i in 1:41){ ##for each year
    height <- sum(growthRate[1:i]) ##total height
    Returns[i] <- nTrees*height ##volume
    prevTrees <- nTrees
    percentDead <- rgamma(1, shape = 1, scale = pDead[i])###what percent of trees will die based on gamma distribution?
    numDead <- (percentDead/100)*prevTrees###number of dead trees
    nTrees <- prevTrees - numDead ##update number of trees
  } ##for each year
  ##returns <- cbind(returns, Returns)
  ##colnames(returns)[length(returns)] <- seed
  
  assets <- vector("numeric", 41)
  assets[1] <- Returns[1]
  
  for (z in 1:40){ ##convert from cumulative volume to change by year
    assets[z+1] <- Returns[z+1] - Returns[z]
  }
  return(assets)
}

cleanFrontier <- function(EF.out.all, EF.ret.all){
  EF.out.all$NumNA <- apply(is.na(EF.out.all), 1, FUN = sum)##number of models where it didn't show up
  EF.out.all <- EF.out.all[EF.out.all$NumNA < 16,] ##remove seed stock showing up in < 50% of models
  EF.out.all[is.na(EF.out.all)] <- 0
  EF.out.all$Mean <- apply(EF.out.all[,-c(1,2,length(EF.out.all))],1,mean) ### mean of weights
  for(R in EF.out.all$RA){###scale each out of 1
    EF.out.all$Mean[EF.out.all$RA == R] <- EF.out.all$Mean[EF.out.all$RA == R]/sum(EF.out.all$Mean[EF.out.all$RA == R])
  }
  
  
  EF.ret.all$Mean <- apply(EF.ret.all[,-1],1,mean, na.rm = TRUE) ###mean of returns
  EF.ret.all$Mean <- EF.ret.all$Mean/max(EF.ret.all$Mean) ##scale return out of 1
  EF.ret.all <- EF.ret.all[,c(1,length(EF.ret.all))]
  colnames(EF.ret.all) <- c("Risk","MeanRet")
  EF.sum <- EF.out.all[,c("Seed","RA","Mean")]

  if(nrow(EF.sum) > 0){
    EF.sum$Site <- SNum
  }
  EF.ret.all$Site <- SNum
  return(list(EF.sum, EF.ret.all))
}

###foreach site
allSites <- foreach(SNum = SiteList, .combine = combineList,.packages = c("scales", "foreach","reshape2","dplyr","magrittr","PortfolioAnalytics")) %dopar% {
  EF.out.all1 <- data.frame(Seed = "SBSmc2", RA = 2)
  EF.out.all2 <- data.frame(Seed = "SBSmc2", RA = 2)
  EF.ret.all1 <- data.frame(RA = 2)
  EF.ret.all2 <- data.frame(RA = 2)
  SSPred <- SSPredAll[SSPredAll$SiteNo == SNum,]
  currBGC <- as.character(SSPred$BGC[1])
  
  SSPred <- merge(SSPred,SuitTable, by.x = "BGC.pred", by.y = "Site", all.x = TRUE)
  modList <- as.character(unique(SSPred$GCM))
  output <- data.frame(Seed = SeedList)
  
  for(mod in modList){
    SSPredMod <- SSPred[SSPred$GCM == mod,]
    returns <- data.frame(Year = seq(2000,2040,1))
    modData1 <- data.frame(Year = seq(2000,2040,1))
    modData2 <- data.frame(Year = seq(2000,2040,1))
    
    for(seed in SeedList){
      SSPredSd <- SSPredMod[SSPredMod$Seed == seed,]
      SSPredSd <- SSPredSd[order(SSPredSd$GCM,SSPredSd$FuturePeriod),]
      SS.sum <- SSPredSd[,c("FuturePeriod","Spp1","Spp2")]
      curr <- data.frame(FuturePeriod = 2000, 
                         Spp1 = SuitTable$Spp1[SuitTable$Site == BGC & SuitTable$Seed == as.character(seed)],
                         Spp2 = SuitTable$Spp2[SuitTable$Site == BGC & SuitTable$Seed == as.character(seed)])
      SS.sum <- rbind(curr, SS.sum)
      SS.sum$Spp1 <- gs2gwVec(SS.sum$Spp1)
      SS.sum$Spp2 <- gs2gwVec(SS.sum$Spp2)
      SS.sum[SS.sum < 0.00354] <- -1
      SS.sum <- SS.sum[complete.cases(SS.sum),]
      
      Spp1 <- simulateGrowth(SS.sum$Spp1)
      Spp2 <- simulateGrowth(SS.sum$Spp2)
      if(is.null(Spp1) & is.null(Spp2)){next}
      
      if(!is.null(Spp1)){
        modData1 <- cbind(modData1, Spp1)
        colnames(modData1)[length(modData1)] <- seed
      }
      if(!is.null(Spp2)){
        modData2 <- cbind(modData2, Spp2)
        colnames(modData2)[length(modData2)] <- seed
      }
      
    }
    
    datNames <- c("modData1","modData2")
    for(i in 1:length(datNames)){
      modData <- get(datNames[i])
      if(ncol(modData) < 3){next}
      returns <- modData
      rownames(returns) <- paste(returns$Year,"-01-01", sep = "")
      returns <- returns[,-1]
      #returns <- returns[1:76,]
      returnsTS <- as.xts(returns)
      
      init.portfolio <- portfolio.spec(assets = colnames(returnsTS))
      nSpp <- length(colnames(returnsTS))
      init.portfolio <- add.constraint(portfolio = init.portfolio, type = "weight_sum", min_sum = 0.9, max_sum = 1.1) ###weights should add to about 1
      init.portfolio <- add.constraint(portfolio = init.portfolio, type = "box", min = rep(0, nSpp), max = rep(0.95, nSpp)) ##set min and max weight for each species
      
      for(q in seq(from = 2, to = 50, by = 2)){
        qu <- add.objective(portfolio=init.portfolio, type="return", name="mean")
        qu <- add.objective(portfolio=qu, type="risk", name="var", risk_aversion = q)
        
        minSD.opt <- optimize.portfolio(R = returnsTS, portfolio = qu, optimize_method = "ROI", trace = TRUE)
        if(q == 2){
          EF.out <-  as.data.frame(minSD.opt[["weights"]])
          EF.out$Seed <- rownames(EF.out)
          EF.out$RA <- q
          EF.ret <- as.data.frame(minSD.opt[["objective_measures"]][["mean"]])
          EF.ret$RA <- q
        }else{
          temp <- as.data.frame(minSD.opt[["weights"]])
          temp$Seed <- rownames(temp)
          temp$RA <- q
          EF.out <- rbind(EF.out, temp)
          temp <- as.data.frame(minSD.opt[["objective_measures"]][["mean"]])
          temp$RA <- q
          EF.ret <- rbind(EF.ret, temp)
        }
        
      }
      if(i == 1){
        EF.out.all1 <- merge(EF.out.all1, EF.out, by = c("Seed","RA"), all = TRUE)
        EF.ret.all1 <- merge(EF.ret.all1, EF.ret, by = "RA", all = TRUE)
      }else{
        EF.out.all2 <- merge(EF.out.all2, EF.out, by = c("Seed","RA"), all = TRUE)
        EF.ret.all2 <- merge(EF.ret.all2, EF.ret, by = "RA", all = TRUE)
      }
      
    }
    
  }
  Spp1Out <- cleanFrontier(EF.out.all1, EF.ret.all1)
  Spp2Out <- cleanFrontier(EF.out.all2, EF.ret.all2)
  outList <- list(Spp1Out[[1]],Spp1Out[[2]], Spp2Out[[1]],Spp2Out[[2]])
  outList
} 

for(i in c(1,3)){
  EFall <- allSites[[i]]   
  EFall <- aggregate(Mean ~ Seed + RA, data = EFall, FUN = mean)
  maxW <- aggregate(Mean ~ Seed, data = EFall, FUN = max)
  EFall <- EFall[EFall$Seed %in% maxW$Seed[maxW$Mean > 0.05],]
  for(R in unique(EFall$RA)){###scale each out of 1
    EFall$Mean[EFall$RA == R] <- EFall$Mean[EFall$RA == R]/sum(EFall$Mean[EFall$RA == R])
  }
  EFret <- aggregate(MeanRet ~ Risk, data = allSites[[i+1]], FUN = mean)
  
  assign(paste("CBSTp",i, sep = ""),
    ggplot(EFall)+
    geom_bar(aes(x = RA, y = Mean, fill = Seed), stat = "identity")+
    geom_line(data = EFret, aes(x = Risk, y = MeanRet), size = 2)+
    scale_fill_discrete()+
    labs(x = "Risk (0 = High Risk, 50 = Low Risk)"))
}

layoutMat <- rbind(c(1,1,1,1,2,2),
                   c(1,1,1,1,NA,NA),
                   c(1,1,1,1,3,3))
library(gridExtra)

grid.arrange(grobs = list(SpeciesPlot,CBSTp3,CBSTp1), layout_matrix = layoutMat)


###maybe for later: https://stackoverflow.com/questions/35631889/align-multiple-plots-with-varying-spacings-and-add-arrows-between-them/35634129#35634129
####Efficient Frontier



###Violin plot
portOut <- allSites$Weights[,-(nSpp+1)]

portOut <- melt(portOut)
portOut$variable <- as.factor(portOut$variable)


ggplot(portOut)+
  geom_violin(aes(x = variable, y = value),colour = "purple", draw_quantiles = c(0.25,0.5,0.75), scale = "width")+
  xlab("Species")+
  ylab("Weight")

###Old Version
library(vioplot)
vioplot(portOut$Fd,portOut$Lw,portOut$Pl,portOut$Bl,portOut$Sx,portOut$Cw,portOut$Py, names = colnames(portOut), col = "purple") ##have to manually change if adding species


##===================================================================================
##                  OLD CODE
##====================================================================================



y <- dgamma(x, shape = 1, scale = 4)
plot(x,y, type = "l")

#opt.ts <- portfolio.optim(returnsTS, pm = 14,covmat = sigma)

#portOutput <- cbind(portOutput, opt.ts$pw)
#colnames(portOutput)[w+1] <- paste("Run",w,sep = "")
covmat = sigma
print(opt_qu)
portOutput.long <- portOutput
portOutput.long <- portOutput.long[,-1]
portOutput.long <- as.data.frame(t(portOutput.long))
portOutput.long <- portOutput.long[complete.cases(portOutput.long),]
vioplot(portOutput.long$V1,portOutput.long$V2,portOutput.long$V3,portOutput.long$V4,portOutput.long$V5,portOutput.long$V6, names = c("Pl","Sx","Cw","Bl","Lw","Fd"))
colnames(portOutput.long) <- c("Pl","Sx","Cw","Bl","Lw","Fd")
write.csv(portOutput.long, file = "Lw_SI30_S2.csv")

frontier <- extractEfficientFrontier(opt_qu, match.col = "StdDev", risk_aversion = seq(10,70,2))
chart.EfficientFrontier(frontier, match.col = "StdDev")
chart.EF.Weights(frontier, match.col = "StdDev")

colMeans(portOutput.long)
##################################################################################
growForest <- function(growth, meanDead){
  height <- sum(growth)
  nTrees <- 100 - rnorm(1, mean = meanDead, sd = 6)
  Returns <- nTrees*height
  return(list(Returns = Returns))
}

param_list = list("growth" = annualDat$Growth, "meanDead" = annualDat$MeanDead)

MCResult <- MonteCarlo(func = growForest, nrep = 100, param_list = param_list, max_grid = 10000)

MakeTable(MCResult, rows = "meanDead", cols = "growth", digits = 2)
rlnorm(1, meanlog = 5, sdlog = 5)

lines(output$year,output$Run1, col = "black", type = "l")
lines(output$year,output$Run2, col = "purple", type = "l")
plot(output$year,output$Run3, col = "red", type = "l")
legend("topleft", c("Original","Changed SIBEC", "Changed Suit"), col = c("black","purple","red"), pch = 15)

randomCorr <- data.frame("Year" = seq(1,100,1))
prev <- 0
for (j in 1:6){
  pests <- vector(mode = "numeric", length = 100)
  for (i in 1:100){
    pests[i] <- prev - rweibull(1,shape = 1, scale = 0.5)
  }
  randomCorr <- cbind(randomCorr, pests)
}
randomCorr <- randomCorr[,-1]
colnames(randomCorr) <- c("Pl", "Fd", "Cw", "Lw", "Bl", "Sx")
require(pcaPP)
corr <- cor.fk(randomCorr)

x <- seq(0,10,0.1)
y <- dweibull(x,shape = 1.5, scale = 1.5)
plot(x,y,type = "l")


print(minSD.opt)
frontier1 <- create.EfficientFrontier(R = returnsTS, portfolio = minSD.portfolio, type = "mean-StdDev", match.col = "StdDev", n.portfolios = 30)
chart.EF.Weights(frontier1, main= "Risk", match.col = "StdDev")

#####Random Portfolio#############################
init.portfolio <- portfolio.spec(assets = colnames(returnsTS))
init.portfolio <- add.constraint(portfolio = init.portfolio, type = "weight_sum", min_sum = 0.99, max_sum = 1.01)
init.portfolio <- add.constraint(portfolio = init.portfolio, type = "box", min = rep(0.01,6), max = rep(0.9, 6))
maxret.portfolio <- add.objective(portfolio=init.portfolio, type="return", name="mean")

opt.maxRet <- optimize.portfolio(R = returnsTS, portfolio = maxret.portfolio, optimize_method = "DEoptim", trace = TRUE, momentargs = momentargs, search_size = 20000)
maxret <- as.numeric(opt.maxRet[["objective_measures"]])

#init.portfolio <- add.constraint(portfolio = init.portfolio, type = "return", return_target = 0.75*maxret)
good.portfolio <- add.objective(portfolio = init.portfolio, type = "risk", name = "StdDev", target = 0.1)
good.portfolio <- add.objective(portfolio=init.portfolio, type="return", name="mean", target = 0.75*maxret)

#good.opt <- optimize.portfolio(R = returnsTS, portfolio = good.portfolio, optimize_method = "DEoptim", trace = TRUE, momentargs = momentargs, search_size = 20000)
rp <- random_portfolios(init.portfolio, 10000, "sample")
good.opt <- optimize.portfolio(R = returnsTS, portfolio = good.portfolio, optimize_method = "random", rp = rp, trace = TRUE)
test <- as.data.frame(good.opt[["weights"]])
colnames(test) <- "Run"
portOutput <- cbind(portOutput, test$Run)
print(good.opt)

if(sigma > 0 & sigma < 1){
  sigma <- sigma*0.5
}

half <- function(x){
  if(x > 0 & x < 1){
    x <- x/2
  }
  return(x)
}

sigma <- sigmaSave
sigma <- as.data.frame(sigma)
sigma$Pl <- lapply(sigma$Pl, half)
sigma$Sx <- lapply(sigma$Sx, half)
sigma$Cw <- lapply(sigma$Cw, half)
sigma$Bl <- lapply(sigma$Bl, half)
sigma$Lw <- lapply(sigma$Lw, half)
sigma$Fd <- lapply(sigma$Fd, half)

sigma$Pl <- as.numeric(sigma$Pl)
sigma$Sx <- as.numeric(sigma$Sx)
sigma$Cw <- as.numeric(sigma$Cw)
sigma$Bl <- as.numeric(sigma$Bl)
sigma$Lw <- as.numeric(sigma$Lw)
sigma$Fd <- as.numeric(sigma$Fd)
sigma <- as.matrix(sigma)
sigma <- as.numeric(sigma)

x <- seq(from = -2, to = 5, by = 0.05)
y <- dgamma(x, shape = 1, scale = 4)
plot(x,y, type = "l")
 