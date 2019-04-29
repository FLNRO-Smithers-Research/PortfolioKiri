##This script uses a Monte Carlo approach to model tree growth wiht SI numbers and predicted suitability 
##It then uses the Markowitz portfolio method to calculate the optimal mix of species.
##Kiri Daust, June 2018

rm(list=ls())

library(stats)
library(rgl)
require(tcltk)
library(bindr)
library(lattice)
require(dplyr)
library(ggplot2)
library(MonteCarlo)
require(openxlsx)
require(MASS)
require(xts)
library(PortfolioAnalytics)
library(ROI.plugin.quadprog)
library(ROI)
library(MASS)
require(pcaPP)
library(tseries)
library(magrittr)
library(foreach)
library(reshape2)
library(scales)

wd=tk_choose.dir()
setwd(wd)

###OPTIONAL: Set up to run loops in parallel###
require(doParallel)
set.seed(123321)
coreNum <- as.numeric(detectCores()-2)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)

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

SSPredAll <- Y3.sub
SSPredAll$GCM <- paste(SSPredAll$GCM,SSPredAll$Scenario, sep = "-")
SSPredAll <- SSPredAll[,-c(2,4,8)]


###import covariance matrix
fullMat <- read.csv("Pl genetic suitability matrix - no assisted migration.csv")
rownames(fullMat) <- fullMat$BECvar
fullMat <- fullMat[,-1]

########################

combineList <- function(...) {##Combine multiple dataframe in foreach loop
  mapply(FUN = rbind, ..., SIMPLIFY=FALSE)
}

SiteList <- as.character(SSPredAll$SiteNo[SSPredAll$BGC == "SBSmc2"])
outRaw <- data.frame(Seed = SeedList, Number = 0)
###foreach site
allSites <- foreach(SNum = SiteList[1:100], .combine = rbind, 
                    .packages = c("foreach","reshape2","dplyr","magrittr","PortfolioAnalytics")) %do% {
    SSPred <- SSPredAll[SSPredAll$SiteNo == SNum,]
    currBGC <- as.character(SSPred$BGC[1])
    
    SSPred <- merge(SSPred,SuitTable, by.x = "BGC.pred", by.y = "Site", all.x = TRUE)
      modList <- as.character(unique(SSPred$GCM))
      SeedList <- unique(SSPred$Seed)
      output <- data.frame(Seed = SeedList)
      
      for(mod in modList){
        SSPredMod <- SSPred[SSPred$GCM == mod,]
        returns <- data.frame(Year = seq(2000,2100,1))
        modData <- data.frame(Year = seq(2000,2100,1))
        plot(0,0, type = "n", xlim = c(1,100), ylim = c(0,50), xlab = "Year", ylab = "Volume")###plot
        
        for(seed in SeedList){
          SSPredSd <- SSPredMod[SSPredMod$Seed == seed,]
          SSPredSd <- SSPredSd[order(SSPredSd$GCM,SSPredSd$FuturePeriod),]
          SS.sum <- SSPredSd[,c("FuturePeriod","HTp_pred")]
          colnames(SS.sum)[2] <- "Growth"
          curr <- data.frame(FuturePeriod = 2000, Growth = SuitTable$HTp_pred[SuitTable$Site == SSPredMod$BGC[1] & SuitTable$Seed == seed])
          SS.sum <- rbind(curr, SS.sum)
          SS.sum$Growth <- (SS.sum$Growth - 0.9)*10
          if(any(SS.sum$Growth < 0)){next}###need to account for different numbers of reps
          outRaw$Number[outRaw$Seed == seed] <- outRaw$Number[outRaw$Seed == seed] + 1
          
          cols <- rainbow(50)
          annualDat <- data.frame("Year" = seq(2000,2100,1))
          
          portOutput <- data.frame("Seed" = SeedList)###set up plot and output
         ## s <- approx(SS.sum$FuturePeriod, SS.sum$Growth, n = 101) ##Smooth SI
          s <- spline(SS.sum$FuturePeriod, SS.sum$Growth, n = 101)
          growthRate <- s[["y"]]
          #lines(growthRate)
          pDead <- 1 - s[["y"]]
          pDead <- rescale(pDead, to = c(0.01,0.1), from = c(0,1))

          nTrees <- 100 ##number of tree to start
          Returns <- numeric(length = 101)

          ###Simulate Growth
          for (i in 1:101){ ##for each year
            height <- sum(growthRate[1:i]) ##total height
            Returns[i] <- nTrees*height ##volume
            prevTrees <- nTrees
            percentDead <- rgamma(1, shape = 1, scale = pDead[i])###what percent of trees will die based on gamma distribution?
            numDead <- (percentDead/100)*prevTrees###number of dead trees
            nTrees <- prevTrees - numDead ##update number of trees
          } ##for each year
          returns <- cbind(returns, Returns)
          colnames(returns)[length(returns)] <- seed

          assets <- vector("numeric", 101)
          assets[1] <- Returns[1]

          for (z in 1:100){ ##convert from cumulative volume to change by year
            assets[z+1] <- Returns[z+1] - Returns[z]
          }
          lines(assets)
          modData <- cbind(modData, assets)
          colnames(modData)[length(modData)] <- seed
          # returns <- cbind(returns, s[["y"]])
          # lines(s[["y"]])
          # colnames(returns)[length(returns)] <- seed
        }
        if(ncol(modData) < 3){next}
        returns <- modData
        rownames(returns) <- paste(returns$Year,"-01-01", sep = "")
        returns <- returns[,-1]
        #returns <- returns[1:76,]
        returnsTS <- as.xts(returns)

        # nSpp <- ncol(returns)
        # sdList <- colnames(returns)
        # sigma <- cor(fullMat, method = "pearson")
        # sigma <- sigma[sdList,sdList]
        # # tempMat <- fullMat[sdList,sdList]
        # # sigma <- cor(tempMat, method = "pearson")
        # sigma <- as.matrix(sigma)
        # momentargs <- list()
        # momentargs$sigma <- sigma ##set up to use in PortfolioAnalystics 
        
        init.portfolio <- portfolio.spec(assets = colnames(returnsTS))
        nSpp <- length(colnames(returnsTS))
        init.portfolio <- add.constraint(portfolio = init.portfolio, type = "weight_sum", min_sum = 0.99, max_sum = 1.01) ###weights should add to about 1
        init.portfolio <- add.constraint(portfolio = init.portfolio, type = "box", min = rep(0, nSpp), max = rep(1, nSpp)) ##set min and max weight for each species
        
        # ret_const <- return_objective(name = "mean")
        # risk_const <- portfolio_risk_objective(name = "var", risk_aversion = 30) ###set risk_aversion
        # qu <- list(ret_const,risk_const)
        qu <- add.objective(portfolio=init.portfolio, type="risk", name="var", risk_aversion = 15) ###set risk_aversion
        qu <- add.objective(portfolio=qu, type="return", name="mean")
        
        ##optimise
        minSD.opt <- optimize.portfolio(R = returnsTS, portfolio = qu, 
                                        optimize_method = "quadprog", trace = TRUE)
        
        test <- as.data.frame(minSD.opt[["weights"]]) ##save weights
        test$Seed <- rownames(test)
        colnames(test)[1] <- mod
        output <- merge(output, test, by = "Seed", all.x = TRUE)
  }
    
  rownames(output) <- output$Seed
  output <- output[,-1]
  output[output < 0.05] <- NA
  output <- output[rowSums(is.na(output)) != ncol(output), ]
  output$SiteNo <- SNum
} 
    
  








  ###foreach model
    
    ###Create current dat
    current <- data.frame(Seed = BGCnoSS, FuturePeriod = 2000, MeanSI = 1)
    
    ##Summarise data- average SI and Suit weighted by SSProb
    SS.sum <- SSPred %>%
      group_by(MergedBGC, FuturePeriod) %>%
      summarise(MeanSI = sum(HTp_pred*(SSprob/sum(SSprob))))
    
    SS.sum <- as.data.frame(SS.sum)
    colnames(SS.sum)[1] <- "Seed"
    SS.sum <- rbind(SS.sum, current)
    SS.sum <- unique(SS.sum)
    SS.sum <- SS.sum[order(SS.sum$Seed,SS.sum$FuturePeriod),]
    
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

###Violin plot
portOut <- allSites$Weights[,-(nSpp+1)]

portOut <- melt(portOut)
portOut$variable <- as.factor(portOut$variable)


ggplot(portOut)+
  geom_violin(aes(x = variable, y = value), draw_quantiles = c(0.25,0.5,0.75), scale = "width")+
  labs(x = "Spp", y = "Weight")

###Old Version
vioplot(portOut$Fd,portOut$Lw,portOut$Pl,portOut$Bl,portOut$Sx,portOut$Cw,portOut$Py, names = colnames(portOut), col = "purple") ##have to manually change if adding species


####Efficient Frontier
EF.sum <- allSites$Frontier
EF.sum <- aggregate(Weight ~ Spp + Risk, EF.sum, FUN = mean)
EF.ret.all <- allSites$Return
EF.ret.all <- aggregate(MeanRet ~ Risk, EF.ret.all, FUN = mean)

ggplot(EF.sum)+
  geom_bar(aes(x = Risk, y = Weight, fill = Spp), stat = "identity")+
  geom_line(data = EF.ret.all, aes(x = Risk, y = MeanRet), size = 2)+
  geom_vline(xintercept = 12)+
  labs(x = "Risk (0 = High Risk, 50 = Low Risk)")+
  ggtitle((""))

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
library(pcaPP)
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
 