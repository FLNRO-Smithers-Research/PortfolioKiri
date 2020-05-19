##This script uses a Monte Carlo approach to model tree growth wiht SI numbers and predicted suitability 
##It then uses the Markowitz portfolio method to calculate the optimal mix of species.
##Python script PortfolioOptimisation.py and C++ script SimGrowth.cpp must be sourced.
##Kiri Daust, May 2020

require(tcltk)
require(dplyr)
require(ggplot2)
require(MASS)
require(magrittr)
require(foreach)
require(reshape2)
require(doParallel)
require(reticulate)
require(Rcpp)
library(gridExtra)
library(data.table)
library(scales)
library(tidyr)

rm(list=ls())
##wd=tk_choose.dir(); setwd(wd)
###source C++ and Python Scripts
#setwd("C:/Users/kirid/Desktop/PortfolioKiri")
#setwd(tk_choose.dir())
sourceCpp("CppFunctions/SimGrowth.cpp")
source_python("PythonFns/PortfolioOptimisation.py") ###make sure you have python as a path variable - easiest way is to install Anaconda

############################################################
####If doing CBST, stop here and go to line 350############
##########################################################



###OPTIONAL: Set up to run loops in parallel###
# require(doParallel)
# set.seed(123321)
# coreNum <- as.numeric(detectCores()-2)
# coreNo <- makeCluster(coreNum)
# registerDoParallel(coreNo, cores = coreNum)


sigma <- read.csv("InputsGit/CovarianceMatrix_Full.csv")
rownames(sigma) <- sigma[,1]
sigma <- sigma[,-1]
Trees <- c("Bl","Cw","Fd","Hw","Lw","Pl","Py","Sx") ##set species to use in portfolio
nSpp <- length(Trees)
treeList <- Trees
sigma <- sigma[Trees,Trees]

repeat.before = function(x) {   # repeats the last non NA value. Keeps leading NA
  ind = which(!is.na(x))      # get positions of nonmissing values
  if(is.na(x[1]))             # if it begins with a missing, add the 
    ind = c(1,ind)        # first position to the indices
  rep(x[ind], times = diff(   # repeat the values at these indices
    c(ind, length(x) + 1) )) # diffing the indices + length yields how often 
}

nsuitF <- file.choose() ##Import suitability table
SuitTable <- read.csv(nsuitF, stringsAsFactors = FALSE)
SuitTable <- unique(SuitTable)

colnames(SuitTable)[2:4] <- c("SS_NoSpace","Spp","Suitability")

SIBEC <- read.csv("InputsGit/BartPredSI.csv", stringsAsFactors = FALSE) ###import SI data
SIBEC <- SIBEC[,-4]
colnames(SIBEC) <- c("SS_NoSpace", "TreeSpp","MeanPlotSiteIndex")

SSPredAll <- read.csv(file.choose(), stringsAsFactors = FALSE) ##Import SS predictions from CCISS tool: must have columns MergedBGC, Source, SS_NoSpace, SSprob, SSCurrent, FuturePeriod, SiteNo
SSPredAll <- SSPredAll[,c("MergedBGC", "Source", "SS_NoSpace", "SSprob", "SSCurrent", 
                          "FuturePeriod", "SiteNo")]
selectBGC <- select.list(choices = sort(unique(SSPredAll$SSCurrent)), graphics = TRUE) ###Select BGC to run for
SSPredAll <- SSPredAll[SSPredAll$SSCurrent == selectBGC,]

SSPredSave <- SSPredAll
SSPredAll <- SSPredSave

#####Randomly select 100 sites - if testing to speed up############
# sites <- as.numeric(as.character(unique(SSPredAll$SiteNo)))
# SiteList <- sample(sites, 50, replace = FALSE)
#################################################################

SiteList <- unique(SSPredAll$SiteNo[grep("Ele",SSPredAll$SiteNo)]) ###select Elephant hill fire
##SiteList <- unique(SSPredAll$SiteNo)
SSPredAll <- SSPredAll[SSPredAll$SiteNo %in% SiteList,]
SSPredAll <- SSPredAll[!is.na(SSPredAll$SSprob),]
SList <- unique(SSPredAll$SiteNo)
########################

combineList <- function(...) {##Combine multiple dataframe in foreach loop
  mapply(FUN = rbind, ..., SIMPLIFY=FALSE)
}

Trees <- c("Bl","Cw","Fd","Hw","Lw","Pl","Py","Sx")
minWt <- c(0,0,0,0,0.05,0.2,0,0.1)
maxWt <- c(0.3,0.3,0.8,0.5,0.4,0.6,1,1)
boundDat <- data.frame(Spp = Trees, Min = minWt, Max = maxWt)


###foreach site
allSitesSpp <- foreach(SNum = SList[1:20], .combine = combineList, 
                       .packages = c("foreach","reshape2","dplyr","magrittr","PortfolioAnalytics", "Rcpp"), 
                       .noexport = c("simGrowthCpp")) %do% {
                         SSPred <- SSPredAll[SSPredAll$SiteNo == SNum,] ###subset
                         
                         ##Merge SIBEC data
                         SIBEC <- SIBEC[SIBEC$TreeSpp %in% Trees,]
                         SSPred <- SSPred[,c(6,3,4)]
                         SSPred <- merge(SSPred, SIBEC, by = "SS_NoSpace", all.x = TRUE)
                         
                         ###Add rows for species with missing SI - usually don't need this but safer to keep in
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
                           add$MeanPlotSiteIndex <- 5 ##Set missing SI
                           SSPred <- rbind(SSPred, add)
                         }
                         SSPred <- SSPred[!is.na(SSPred$TreeSpp),]
                         colnames(SSPred)[4] <- "Spp"
                         
                         ##Add suitability
                         SSPred <- merge(SSPred, SuitTable, by = c("SS_NoSpace","Spp"), all.x = TRUE)
                         SSPred$Suitability[is.na(SSPred$Suitability)] <- 5
                         
                         ###Create current data
                         current <- SIBEC[SIBEC$SS_NoSpace == selectBGC,] %>% 
                           merge(SuitTable, by.x = c("TreeSpp","SS_NoSpace"), by.y = c("Spp","SS_NoSpace"), all.x = TRUE) %>%
                           unique()
                         
                         ###check that there aren't errors in the table
                         temp <- aggregate(SS_NoSpace ~ TreeSpp, current, FUN = length)
                         if(any(temp$SS_NoSpace > 1)){
                           stop("There are partial duplicates in the suitablity table. Please fix them. :)")
                         }
                         
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
                         
                         SS.sum <- as.data.frame(SS.sum) %>% 
                           rbind(current)
                         SS.sum$MeanSI[SS.sum$MeanSuit == 4] <- 5
                         SS.sum$MeanSI[SS.sum$MeanSuit == 5] <- 0
                         SS.sum$MeanSuit[SS.sum$MeanSuit == 5] <- 4
                         SS.sum <- unique(SS.sum)
                         SS.sum <- SS.sum[order(SS.sum$Spp,SS.sum$FuturePeriod),]
                         
                         cols <- rainbow(length(treeList))
                         annualDat <- data.frame("Year" = seq(2000,2100,1))
                         
                         plot(0,0, type = "n", xlim = c(1,100), ylim = c(0,3000), xlab = "Year", ylab = "Volume")###plot
                         
                         eff_front <- foreach(w = 1, .combine = combineList) %do% { ##number of iterations
                           output <- data.frame("year" = annualDat$Year)
                           growthSim <- data.frame(Spp = character(), Year = numeric(), Returns = numeric())
                           
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
                             
                             Returns <- simGrowthCpp(DF = annualDat)
                             lines(Returns, type = 'l', col = cols[k])##plot
                             tmpR <- c(0,Returns)
                             assets <- Returns - tmpR[-length(tmpR)]
                             temp <- data.frame(Spp = rep(treeList[k],101), Year = 1:101, Returns = Returns)
                             growthSim <- rbind(growthSim, temp)
                             output <- cbind(output, assets)
                           } ## for each tree species
                           
                           legend("topleft", treeList, col = cols, pch = 15)
                           
                           colnames(output) <- c("Year", treeList)
                           
                           ####Portfolio#######################################
                           returns <- output
                           rownames(returns) <- returns[,1]
                           returns <- returns[,-1]
                           ###only include species with mean return > 1 in portfolio
                           use <- colnames(returns)[colMeans(returns) > 1] ###should probably be higher
                           returns <- returns[,use]
                           sigma2 <- as.data.frame(cor(returns)) ###to create cov mat from returns
                           ####sigma2 <- sigma[use,use] ###to use pre-made cov mat
                           target <- set_target(returns, sigma2) ###find range of efficient frontier
                           
                           temp <- boundDat[boundDat$Spp %in% colnames(returns),]
                           ef <- ef_weights(returns, sigma2, target,0,1,0.1) ###change to "mean" for maximising return - main python function
                           #maxSharpe <- max_sharpe_ratio(returns, sigma2,0)
                           ef_w <- ef[[1]]
                           
                           ###add in spp not used
                           notUse <- setdiff(Trees, colnames(ef_w))
                           temp <- as.data.frame(matrix(data = 0, nrow = length(target), ncol = length(notUse)))
                           colnames(temp) <- notUse
                           ##ef_w <- rbind(ef_w, maxSharpe)
                           ef_w <- cbind(ef_w, temp)
                           ef_w$Sd <- c(ef[[2]])
                           ef_w$Return <- rescale(target, to = c(0,1))
                           ef_w <- ef_w[,c(Trees, "Sd","Return")]
                           
                           # ##to plot individual portfolios
                           # #############################################
                           # efAll <- ef_w[ef_w$Return != "mSharpe",]
                           # efAll$Sd <- efAll$Sd/max(efAll$Sd)
                           # efAll$Return <- as.numeric(efAll$Return)
                           # efAll <- melt(efAll, id.vars = "Return")
                           # print(ggplot(efAll[efAll$variable != "Sd",])+
                           #   geom_bar(aes(x = Return, y = value, fill = variable),size = 0.00001, col = "black", stat = "identity")+
                           #   colScale +
                           #   geom_line(data = efAll[efAll$variable == "Sd",], aes(x = Return, y = value))+
                           #   #theme(legend.position = "none")+
                           #   xlab("Volatility"))
                           # ############################################
                           ef_w$It <- w
                           
                           sigmaOut <- sigma2
                           sigmaOut$Row <- rownames(sigmaOut)
                           sigmaOut <- melt(sigmaOut)
                           colnames(sigmaOut)[2] <- "Column"
                           sigmaOut$It <- 4
                           list(frontier = ef_w, sim = growthSim, sig = sigmaOut)
                         }
                         
                         eff_front2 <- melt(eff_front[['frontier']], id.vars = c("Return","It"))
                         eff_front2 <- dcast(eff_front2, Return ~ variable, fun.aggregate = mean)
                         eff_front2$SiteNo <- SNum
                         growthSim <- eff_front$sim
                         growthSim$SiteNo <- SNum
                         sigmaOut <- eff_front$sig
                         sigmaOut$SiteNo <- SNum
                         list(frontier = eff_front2, sim = growthSim, sig = sigmaOut)
                       }

###The optimisation works much better with set returns and optimising stdev - the opposite of what we're plotting
###We have to therefore round Sd values so that they're not unequally spaced when we combine and average

##############################################
efAll <- allSitesSpp[['frontier']] ###portfolio output
efAll$Return <- round(efAll$Return, digits = 2)
efAll <- efAll[,-length(efAll)]
efAll <- melt(efAll, id.vars = "Return") %>% dcast(Return ~ variable, fun.aggregate = mean)
efAll <- efAll[complete.cases(efAll),]
efAll$Return <- efAll$Return/max(efAll$Return) ##standardise return
ef1 <- c(efAll$Sd,1)
ef2 <- c(0, efAll$Sd)
temp <- ef1 - ef2
temp <- temp[-1]
efAll <- efAll[temp > 0,]
efAll <- melt(efAll, id.vars = "Sd")

myColours <- c("red","pink", "orange","yellow","green","blue","magenta","darkgoldenrod")
names(myColours) <- levels(factor(Trees))
colScale <- scale_fill_manual(name = "variable", values = myColours)

ggplot(efAll[efAll$variable != "Return",])+
  geom_area(aes(x = Sd, y = value, fill = variable),size = 0.00001, col = "black", stat = "identity")+
  colScale +
  geom_line(data = efAll[efAll$variable == "Return",], aes(x = Sd, y = value))+
  scale_x_reverse() +
  xlab("Volatility")+
  ggtitle("EleHill IDFdk3")

#############################################################3
###plot by return#####################
efAll <- allSitesSpp[['frontier']] ###portfolio output
efAll$Return <- round(efAll$Return, digits = 2)
efAll <- efAll[,-length(efAll)]
efAll <- melt(efAll, id.vars = "Return") %>% dcast(Return ~ variable, fun.aggregate = mean)
efAll <- apply(efAll, 2, repeat.before) %>% as.data.frame()
efAll$Return <- efAll$Return/max(efAll$Return) ##standardise return
efAll[is.na(efAll)] <- 0
efAll <- melt(efAll, id.vars = "Return")

myColours <- c("red","pink", "orange","yellow","green","blue","magenta","darkgoldenrod")
names(myColours) <- levels(factor(Trees))
colScale <- scale_fill_manual(name = "variable", values = myColours)

ggplot(efAll[efAll$variable != "Sd",])+
  geom_area(aes(x = Return, y = value, fill = variable),size = 0.00001, col = "black", stat = "identity")+
  colScale +
  geom_line(data = efAll[efAll$variable == "Sd",], aes(x = Return, y = value))+
  scale_x_reverse() +
  xlab("Volatility")+
  ggtitle("Plat IDFdk3")



####if want to add sharpe ratio to plot - will need to uncomment a bunch of stuff above
maxS_plot <- ggplot(sharpe)+
  geom_bar(aes(x = "", y = Weight, fill = Spp),size = 0.00001, col = "black", stat = "identity")+
  colScale+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  xlab("")

layoutMat <- rbind(c(1,1,1,2),
                   c(1,1,1,2),
                   c(1,1,1,2))

grid.arrange(grobs = list(ef_plot,maxS_plot), layout_matrix = layoutMat)

###investigate simulated growth
simGrowth <- allSitesSpp[['sim']]
simGrowth <- dcast(simGrowth, SiteNo+Year ~ Spp, fun.aggregate = mean, value.var = "Returns")
simGrowth <- melt(simGrowth, id.vars = c("SiteNo","Year"))
simGrowth$g <- interaction(simGrowth$SiteNo, simGrowth$variable)

myColours <- c("red","pink", "orange","yellow","green","blue","magenta","darkgoldenrod")
names(myColours) <- levels(factor(Trees))
colScale <- scale_colour_manual(name = "variable", values = myColours)
ggplot(simGrowth, aes(x = as.numeric(Year), y = value, colour = variable, group = g))+
  geom_line(alpha = 0.6)+
  colScale

####mean cov matrix
sig <- allSitesSpp$sig
sig <- dcast(sig, Row ~ Column, value.var = "value", fun.aggregate = mean)


################################################################################################
####CBST Portfolio############
################################################################################################

####function to calculate parameters for gs2gw function given cutoff (default 0.9)
###for CBST Portfolio
setGW <- function(k = 0.9){
  dat <- data.frame(x = c(0.85,0.87,k,0.99,1), y = c(0,0,0.01,0.95,1))
  cont <- nls.control(maxiter = 500, minFactor = 1/5000, warnOnly = TRUE)
  fit <- nls(y ~ a*exp(x*b), data = dat, start = list(a = 5e-25, b = 55.9), control = cont)
  return(coef(fit))
}
params <- setGW(0.96)
####plot gs2gw function to investigate############3
x <- seq(0.85,1,by = 0.001)
y <- gs2gw(x, params[1], params[2])
plot(x,y, type = "l")
#####################################

modelYears <- 60
modPeriod <- c(2025,2055)

###read in CBST data
setwd("./CBSTPortfolio")
Trees <- c("Bl","Cw","Fdi","Hw","Lw","Pli","Py","Sx") ##set species to use in portfolio
dir <- "CBSTNoMigration"
files <- list.files(dir)
allSppDat <- foreach(Spp = Trees, .combine = rbind) %do% {
  fName <- paste(dir, files[grep(paste(Spp,"_",sep = ""),files)], sep = "/")
  temp <- fread(fName, data.table = F)
  temp$Spp <- Spp
  temp
}
colnames(allSppDat) <- c("Site","Seed","Height", "Spp")


##import BGC prediction by model
SSPredAll <- read.csv(file.choose())##import BGC predictions from CCISS script
BGC <- "IDFdk3" ###select BGC
# SSPredAll <- separate(SSPredAll, GCM, c("GCM","Scenario"), sep = "_")
# SSPredAll$GCM <- paste(SSPredAll$GCM, SSPredAll$Scenario, sep = "-")
# SSPredAll <- SSPredAll[,-2]
SSPredAll <- SSPredAll[SSPredAll$BGC == BGC,]
SSPredAll <- SSPredAll[grep("Ele",SSPredAll$SiteNo),]
SSPredAll <- SSPredAll[SSPredAll$FuturePeriod %in% modPeriod,]

# SSPredAll <- SSPredAll[SSPredAll$SiteNo %in% SiteList,]


###foreach spp - not yet working with all species, do each species individually
##allSitesSpp <- foreach(Spp = Trees, .combine = combineList) %do% {

  Spp = "Sx"
  SuitTable <- allSppDat[allSppDat$Spp == Spp,]
  SuitTable <- SuitTable[,-4]
  SSPredAll$SiteNo <- as.character(SSPredAll$SiteNo)
  SSPredAll <- SSPredAll[order(SSPredAll$SiteNo,SSPredAll$GCM),]
  SiteList <- unique(SSPredAll$SiteNo)
  
  
  simulateGrowth <- function(x, nYears = modYears){ ###remove any not suitable in any time period
    if(any(x < 0) || max(x) < 0.6){return(NULL)}
    else{
      s <- spline(c(2000,modPeriod), x, n = nYears+1)
      growthRate <- s[["y"]]
      pDead <- 1 - s[["y"]]
      pDead <- rescale(pDead, to = c(0.01,0.1), from = c(0,1)) %>% multiply_by(100)
      
      annualDat <- data.frame("Year" = seq(2000,2000+nYears,1), "Growth" = growthRate, "MeanDead" = pDead, "NoMort" = rep(25, nYears+1)) ##create working data
      
      Returns <- simGrowthCBST(DF = annualDat)
      
      return(Returns)
    }
    
  }
  
  worker.init <- function(){
    Rcpp::sourceCpp("../CppFunctions/SimGrowth.cpp")
    reticulate::source_python("../PythonFns/PortfolioOptimisation.py")
  }
  
  
  
  require(doParallel)
  cl <- makePSOCKcluster(detectCores()-2)
  clusterCall(cl, worker.init)
  registerDoParallel(cl)
  
  allSites <- foreach(SNum = SiteList[1:5], .combine = rbind, .packages = c("reshape2","Rcpp","magrittr","scales","reticulate"), .noexport = 
                        c("gs2gw", "simGrowthCBST","simGrowthCpp")) %do% {
    
    reticulate::source_python("../PythonFns/PortfolioOptimisation.py")
    
    SSPred <- SSPredAll[SSPredAll$SiteNo == SNum,]
    currBGC <- as.character(SSPred$BGC[1])
    
    SSPred <- merge(SSPred,SuitTable, by.x = "BGC.pred", by.y = "Site", all.x = TRUE)
    ##SSPred <- SSPred[complete.cases(SSPred),]
    
    modList <- as.character(unique(SSPred$GCM))
    output <- data.frame(Return = numeric(), variable = character(), value = numeric(), Model = character())
    temp <- unique(SSPred$Seed)
    #modList <- modList[grep("rcp85",modList)] ##select just rcp8.5 models
    grList <- list()
    
    for(mod in modList){
      cat("Model",mod,"\n")
      SSPredMod <- SSPred[SSPred$GCM == mod,]
      if(any(is.na(SSPredMod$Seed))){next}
      SeedList <- unique(SSPredMod$Seed)
      returns <- data.frame(Year = seq(2000,2000+modelYears,1))
      modData <- data.frame(Seed = character(), Year = numeric(), Returns = numeric())
      
      for(seed in SeedList){
        SSPredSd <- SSPredMod[SSPredMod$Seed == seed,]
        SSPredSd <- SSPredSd[order(SSPredSd$GCM,SSPredSd$FuturePeriod),]
        SS.sum <- SSPredSd[,c("FuturePeriod","Height")]
        curr <- data.frame(FuturePeriod = 2000, 
                           Height = SuitTable$Height[SuitTable$Site == BGC & SuitTable$Seed == as.character(seed)])
        SS.sum <- rbind(curr, SS.sum)
        SS.sum$Height <- gs2gw(SS.sum$Height, as.numeric(params[1]), as.numeric(params[2]))
        ##SS.sum[SS.sum < 0.005] <- -1
        SS.sum <- SS.sum[complete.cases(SS.sum),]
        
        Returns <- simulateGrowth(SS.sum$Height, nYears = modelYears)
        if(is.null(Returns)){next}
        tmpR <- c(0,Returns)
        assets <- Returns - tmpR[-length(tmpR)]
        dat <- data.frame(Seed = rep(seed,modelYears+1), Year = 1:(modelYears+1), Returns = assets)
        modData <- rbind(modData, dat)
      }
      modData <- dcast(modData, Year ~ Seed, value.var = "Returns")
      returns <- modData
      rownames(returns) <- returns[,1]
      returns <- returns[,-1]
      sigma2 <- as.data.frame(cor(returns)) ###to create cov mat from returns
      ##target <- set_target(returns, sigma2) ###find range of efficient frontier
      target <- seq(0.1,0.9,by = 0.02)
      
      ef <- ef_weights_cbst(returns, sigma2, target, 0, 1, 0.01) ###change to "mean" for maximising return - main python function
      ef_w <- as.data.frame(ef[[1]])
      ef_w$Sd <- c(ef[[2]])
      ###########################################
      dat <- ef_w
      dat$Return <- target
      
      #####graph individual model############
      # datRet <- melt(dat, id.vars = "Sd")
      # datRet <- datRet[datRet$variable != "Return",]
      # datRet$Sd <- round(datRet$Sd, digits = 2)
      # datRet <- dcast(datRet, Sd ~ variable, fun.aggregate = mean) %>% melt(id.vars = "Sd")
      # for(R in unique(datRet$Sd)){###scale each out of 1
      #   datRet$value[datRet$Sd == R] <- datRet$value[datRet$Sd == R]/sum(datRet$value[datRet$Sd == R])
      # }
      # grList[[mod]] <- (ggplot(datRet)+
      #   geom_area(aes(x = Sd, y = value, fill = variable),size = 0.00001, col = "black", stat = "identity")+
      #   scale_x_reverse()+
      #   ggtitle(mod))
      #####
      
      dat$Sd <- round(dat$Sd, digits = 2)
      ##dat <- melt(dat, id.vars = "Sd") %>% dcast(Sd ~ variable, fun.aggregate = mean)
      # dat <- merge(temp,dat, by = "Sd", all = T)
      # dat <- apply(dat,2, repeat.before) %>% as.data.frame()
      # dat[is.na(dat)] <- 0
      # dat <- melt(dat, id.vars = "Sd")
      dat <- melt(dat, id.vars = "Return")
      dat$Model <- mod
      output <- rbind(output, dat)
    }
    output$SiteNo <- SNum
    output
  } 
  
 ###To look at each site individually
  # grList <- foreach(SNum = SiteList,.combine = c) %do% {
  #   dat <- allSites[allSites$SiteNo == SNum,]
  #   dat <- dat[!is.nan(dat$value),]
  #   dat <- dcast(dat, Return ~ variable, fun.aggregate = function(x){sum(x)/(15)})
  #   dat <- dat[,colMeans(dat, na.rm = T) > 0.02]
  #   dat$Sd <- round(dat$Sd,digits = 1)
  #   dat <- melt(dat, id.vars = "Sd") %>% dcast(Sd ~ variable, fun.aggregate = mean)
  #   dat <- melt(dat, id.vars = "Sd")
  #   ret <- dat[dat$variable == "Return",]
  #   ret$value <- ret$value/max(ret$value)
  #   dat <- dat[dat$variable != "Return",]
  #   for(R in unique(dat$Sd)){###scale each out of 1
  #     dat$value[dat$Sd == R] <- dat$value[dat$Sd == R]/sum(dat$value[dat$Sd == R])
  #   }
  #   colnames(dat) <- c("Sd","Seed","Proportion")
  #   tmp <- list(ggplot(dat)+
  #     geom_area(aes(x = Sd, y = Proportion, fill = Seed),size = 0.00001, col = "black", stat = "identity")+
  #     geom_line(data = ret, aes(x = Sd, y = value))+
  #     scale_x_reverse())
  #   tmp
  # }
  # 
  # layoutMat <- rbind(c(1,2,3),
  #                    c(4,5,6),
  #                    c(7,8,9),
  #                    c(10,11,12),
  #                    c(13,14,15))
  # 
  # grid.arrange(grobs = grList[1:15], layout_matrix = layoutMat)
  
  output <- allSites[complete.cases(allSites),]
  dat <- output
  dat <- dat[!is.nan(dat$value),]
  dat <- dcast(dat, Return ~ variable, fun.aggregate = function(x){sum(x)/(30*length(SiteList))})
  dat <- dat[,colMeans(dat, na.rm = T) > 0.02]
  dat$Sd <- round(dat$Sd,digits = 1)

  datRet <- melt(dat, id.vars = "Return")
  datRet <- datRet[datRet$variable != "Sd",]
  for(R in unique(datRet$Return)){###scale each out of 1
    datRet$value[datRet$Return == R] <- datRet$value[datRet$Return == R]/sum(datRet$value[datRet$Return == R])
  }
  ggplot(datRet)+
    geom_area(aes(x = Return, y = value, fill = variable),size = 0.00001, col = "black", stat = "identity")+
    ggtitle("CBST by Return")
  
  dat <- melt(dat, id.vars = "Sd") %>% dcast(Sd ~ variable, fun.aggregate = mean)
  dat <- melt(dat, id.vars = "Sd")
  ret <- dat[dat$variable == "Return",]
  ret$value <- ret$value/max(ret$value)
  dat <- dat[dat$variable != "Return",]
  for(R in unique(dat$Sd)){###scale each out of 1
    dat$value[dat$Sd == R] <- dat$value[dat$Sd == R]/sum(dat$value[dat$Sd == R])
  }
  colnames(dat) <- c("Sd","Seed","Proportion")
  ggplot(dat)+
    geom_area(aes(x = Sd, y = Proportion, fill = Seed),size = 0.00001, col = "black", stat = "identity")+
    geom_line(data = ret, aes(x = Sd, y = value))+
    scale_x_reverse()+
    ggtitle(paste("CBST ",Spp, sep = ""))

#}

########################OLD CODE#################################

####################################################################

titles <- c("CBST Pine","", "CBST Spruce")
zoneCols <- read.csv("PortfolioColours.csv")
for(i in c(1,3)){
  EFall <- allSites[[i]]   
  EFall <- aggregate(Mean ~ Seed + RA, data = EFall, FUN = mean)
  maxW <- aggregate(Mean ~ Seed, data = EFall, FUN = max)
  EFall <- EFall[EFall$Seed %in% maxW$Seed[maxW$Mean > 0.06],]
  for(R in unique(EFall$RA)){###scale each out of 1
    EFall$Mean[EFall$RA == R] <- EFall$Mean[EFall$RA == R]/sum(EFall$Mean[EFall$RA == R])
  }
  EFret <- aggregate(MeanRet ~ Risk, data = allSites[[i+1]], FUN = mean)
  EFall <- merge(EFall,zoneCols, by.x = "Seed", by.y = "BGC", all.x = TRUE)
  EFall$Labs <- NA
  EFall <- foreach(x = unique(as.character(EFall$Seed)), .combine = rbind)%do%{
    temp <- EFall[EFall$Seed == x,]
    temp$Labs[temp$Mean == max(temp$Mean)] <- x
    temp
  }
  EFall$Seed <- factor(EFall$Seed, levels = sort(unique(as.character(EFall$Seed))))
  
  
  assign(paste("CBSTp",i, sep = ""),
         ggplot(EFall)+
           geom_bar(aes(x = RA, y = Mean, fill = Colour), size = 0.00001, col = 'black', stat = "identity")+
           scale_fill_identity()+
           geom_line(data = EFret, aes(x = Risk, y = MeanRet), size = 1)+
           geom_label(aes(x = RA, y = Mean, label = Labs, group = Colour), position = position_stack(vjust = 0.5),size =2)+
           labs(x = "Risk (0 = High Risk, 50 = Low Risk)", title = titles[i])
  )
}

layoutMat <- rbind(c(1,2,3),
                   c(4,5,6),
                   c(7,8,9))
layoutMat <- rbind(c(1,2),
                   c(3,4),
                   c(5,6))
library(gridExtra)

pdf(file = "CBSTByMod2.pdf", paper = "letter", pointsize = 8)
grid.arrange(grobs = graphList[7:12], layout_matrix = layoutMat)
dev.off()


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

# testFn <- function(annualDat){
#   Returns <- vector(mode = "numeric", length = 101)
#   nTrees <- 100 ##number of tree to start
#   for (i in 1:101){ ##for each year
#     height <- sum(annualDat$Growth[1:i]) ##total height
#     Returns[i] <- nTrees*height ##volume
#     if(runif(1, min = 0, max = 100) > annualDat$NoMort[i]){##will there be mortality?
#       prevTrees <- nTrees
#       percentDead <- rgamma(1, shape = 1, scale = annualDat$MeanDead[i])###what percent of trees will die based on gamma distribution?
#       numDead <- (percentDead/100)*prevTrees###number of dead trees
#       nTrees <- prevTrees - numDead ##update number of trees
#     }
#   } ##for each year
#   return(Returns)
# }
# 
# cppFunction("double testProb(int n, double shape, double scale){
#               return(Rcpp::rgamma(n,shape,scale)[0]);
#             }")
# cppFunction("double testUnif(int n, double min, double max){
#               return(Rcpp::runif(n,min,max)[0]);
#             }")
# 
# set.seed(42); testFn(annualDat)
# set.seed(42); simGrowthCpp(annualDat)
# 
# 
# set.seed(314152)
# simGrowthCpp(annualDat)
# testFn(annualDat)
if(SNum == SList[1] && w == 1){
  target <- seq(0.1,1.2,by = 0.05) ###if maximising return
  ef <- ef_weights(returns, sigma2, target,0,"mean") ###change to "mean" for maximising return
  temp <- ef[[2]]
  indMax <- which(temp == max(temp))
  temp <- temp[1:indMax]
  d1 <- c(temp,0) - c(0,temp)
  d2 <- c(d1,0) - c(0,d1)
  d2[c(1,length(d2))] <- 0
  indMin <- which(d2 == max(d2))
  target <- seq(0.05*(indMin-1),0.05*indMax,by = 0.03)
}
target <- seq(0.55,0.75,by = 0.02)
target <- seq(0.65, 0.8, by = 0.01)
target <- seq(0.3,0.65, by = 0.02)


