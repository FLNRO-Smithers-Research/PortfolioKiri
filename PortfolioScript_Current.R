##This script uses a Monte Carlo approach to model tree growth wiht SI numbers and predicted suitability 
##It then uses the Markowitz portfolio method to calculate the optimal mix of species.

.libPaths("E:/R packages")
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
library(ROI)
library(MASS)
require(pcaPP)
library(tseries)
library(magrittr)
library(foreach)

wd=tk_choose.dir()
setwd(wd)

repeat.before = function(x) {   # repeats the last non NA value. Keeps leading NA
  ind = which(!is.na(x))      # get positions of nonmissing values
  if(is.na(x[1]))             # if it begins with a missing, add the 
    ind = c(1,ind)        # first position to the indices
  rep(x[ind], times = diff(   # repeat the values at these indices
    c(ind, length(x) + 1) )) # diffing the indices + length yields how often 
}

sigma <- read.csv("SuitCovMat.csv")
rownames(sigma) <- sigma[,1]
sigma <- sigma[,-1]
Trees <- c("Fd","Lw","Pl","Py","Sx") ##set species to use in portfolio
treeList <- Trees
sigma <- sigma[Trees,Trees]
sigma <- as.matrix(sigma)
momentargs <- list()
momentargs$sigma <- sigma ##set up to use in PortfolioAnalystics 

nsuitF <- file.choose() ##Import suitability table
SuitTable <- read.csv(nsuitF, stringsAsFactors = FALSE)
SuitTable <- unique(SuitTable)
colnames(SuitTable)[2] <- "SS_NoSpace"

SIBEC <- read.xlsx("SIBEC_for_Portfolio2.xlsx", sheet = 1) ##Import SIBEC data

SSPredAll <- read.csv("CaribooSSIDFdk3_ByMod.csv", stringsAsFactors = FALSE) ##Import SS predictions from CCISS tool
mods <- unique(SSPredAll$Model)
SSPredAll <- SSPredAll[!is.na(SSPredAll$SSprob),]
SSPredAll <- SSPredAll[SSPredAll$SSCurrent == "IDFdk3/01",] ###Current site series

fires <- c("EleHill","Plat","WillLk")

pdf("Vioplots_WL_6.pdf", width = 8, height = 16, paper = "letter")
layout(rbind(c(1,2), c(3,4),c(5,6),c(7,8),c(9,10)),widths=c(1.2,1.2), heights =c(1,1,1,1,1), respect=TRUE) ## set up plot areas
layout(rbind(c(1,2), c(3,4),c(5,6)),widths=c(1,1), heights =c(1,1,1), respect=TRUE)
par(mai = c(0.3, 0.3, 0.3, 0.3))

#for(a in 26:30){ ## portfolio for each individual model
    
  for (x in 1:3){ ##for each location
    #SSPred <- SSPredAll[SSPredAll$Model == mods[a],-8]
    SSPred <- SSPredAll
    EF.out.all <- rep(1:25, each = 5)
    EF.ret.all <- 1:25
    
    ##Merge SIBEC data
    SIBEC <- SIBEC[SIBEC$TreeSpp %in% Trees,]
    SSPred <- SSPred[SSPred$SiteNo == x,]
    SSPred <- SSPred[,c(6,3,4)]
    SSPred <- merge(SSPred, SIBEC[,c(9,7,5)], by = "SS_NoSpace", all.x = TRUE)
    
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
    
    add$MeanPlotSiteIndex <- 10 ##Set missing SI to 10
    SSPred <- rbind(SSPred, add)
    SSPred <- SSPred[!is.na(SSPred$TreeSpp),]
    colnames(SSPred)[5] <- "Spp"
    
    ##Add suitability
    SSPred <- merge(SSPred, SuitTable, by = c("SS_NoSpace","Spp"), all.x = TRUE)
    SSPred$Suitability[is.na(SSPred$Suitability)] <- 5
    
    ###Create current data
    current <- SIBEC[SIBEC$SS_NoSpace == "IDFdk3/01",c(9,7,5)]
    current <- merge(current, SuitTable, by.x = c("TreeSpp","SS_NoSpace"), by.y = c("Spp","SS_NoSpace"), all.x = TRUE)
    current[is.na(current)] <- 5
    current[current == 0] <- 10
    current <- data.frame(Spp = current$TreeSpp, FuturePeriod = 2000, MeanSI = current$MeanPlotSiteIndex, MeanSuit = current$Suitability)
    
    missing <- Trees[!Trees %in% current$Spp]
    new <- current[rep(1,length(missing)),]
    new$Spp <- missing
    new$MeanSI <- 10
    new$MeanSuit <- 5
    current <- rbind(current, new)
    
    ##Summarise data- average SI and Suit weighted by SSProb
    SS.sum <- SSPred %>%
      group_by(Spp, FuturePeriod) %>%
      summarise(MeanSI = sum(MeanPlotSiteIndex*(SSprob/sum(SSprob))), MeanSuit = round(sum(Suitability*(SSprob/sum(SSprob))), digits = 0))
    
    SS.sum <- as.data.frame(SS.sum)
    SS.sum <- rbind(SS.sum, current)
    SS.sum$MeanSuit[SS.sum$MeanSuit == 5] <- 4
    SS.sum <- unique(SS.sum)
    SS.sum <- SS.sum[order(SS.sum$Spp,SS.sum$FuturePeriod),]
    
    cols <- rainbow(5)
    annualDat <- data.frame("Year" = seq(2000,2100,1))

    portOutput <- data.frame("Species" = c("Fd", "Lw", "Pl", "Py", "Sx"))
    plot(0,0, type = "n", xlim = c(1,100), ylim = c(0,3000), xlab = "Year", ylab = "Volume")
    
    for (w in 1:50){ ##number of iterations
      output <- data.frame("year" = annualDat$Year)
      
      for (k in 1:length(Trees)){ ##for each tree
        DatSpp <- SS.sum[SS.sum$Spp == treeList[k],]
        SuitProb <- data.frame("Suit" = c(1,2,3,4), "ProbDead" = c(0.5,0.5,1,4), 
                               "NoMort" = c(70,60,50,30)) ####ProbDead- out of 100 trees, how many will die each year at each suitability. NoMort- Percent of time no mortality
        dat <- data.frame("Period" = c(2000,2025,2055,2085), "SIBEC" = DatSpp$MeanSI, "Suit" = DatSpp$MeanSuit)
        dat$SIBEC <- dat$SIBEC/50 ##for mean annual increment
        dat <- merge(dat, SuitProb, by = "Suit")
        s <- approx(dat$Period, dat$SIBEC, n = 101) ##Smooth SI
        p <- approx(dat$Period, dat$ProbDead, n = 101) ###Smooth Prob Dead
        m <- approx(dat$Period, dat$NoMort, n = 101) ##Smooth No Mort
        
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
              percentDead <- rgamma(1, shape = 1, scale = annualDat$MeanDead[i])###what percent of trees will die?
              numDead <- (percentDead/100)*prevTrees
              nTrees <- prevTrees - numDead ##update number of trees
            }
          } ##for each year
          modData <- cbind(modData, Returns)
        
        
        modData <- modData[,-1]
        lines(modData, type = 'l', col = cols[k])##graphs
        
        assets <- vector("numeric", 101)
        assets[1] <- modData[1]
        for (z in 1:100){ ##convert from cumulative volume to change by year
          assets[z+1] <- modData[z+1] - modData[z]
        }
        output <- cbind(output, assets)
      } ## for each tree species
      
      legend("topleft", treeList, col = cols, pch = 15)
      
      colnames(output) <- c("Year", "Fd", "Lw", "Pl", "Py", "Sx")
      output$Year <- paste(output$Year,"-01-01", sep = "")
      
      ####Portfolio#######################################
      returns <- output
      rownames(returns) <- returns[,1]
      returns <- returns[,-1]
      returnsTS <- as.xts(returns) ###convert to time series object
      
      ##set portfolio constraints
      init.portfolio <- portfolio.spec(assets = colnames(returnsTS))
      init.portfolio <- add.constraint(portfolio = init.portfolio, type = "weight_sum", min_sum = 0.99, max_sum = 1.01) ###weights should add to about 1
      init.portfolio <- add.constraint(portfolio = init.portfolio, type = "box", min = c(0,0,0,0,0), max = c(1,1,1,1,1)) ##set min and max weight for each species
      
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
      qu <- add.objective(portfolio=qu, type="risk", name="var", risk_aversion = 10)
      
      ##optimise
      minSD.opt <- optimize.portfolio(R = returnsTS, portfolio = qu, optimize_method = "ROI", trace = TRUE, momentargs = momentargs)
      test <- as.data.frame(minSD.opt[["weights"]]) ##save weights
      colnames(test) <- "Run"
      portOutput <- cbind(portOutput, test$Run)
      
      EF.out.all <- cbind(EF.out.all, EF.out)
      EF.ret.all <- cbind(EF.ret.all, EF.ret)
      
    }  
    
    ###calculate average of efficient portfolio weights and their returns
    EF.ret.all$Mean <- apply(EF.ret.all[,-1],1,mean)
    EF.ret.all$EF.ret.all <- 2*EF.ret.all$EF.ret.all
    EF.ret.all$Mean <- EF.ret.all$Mean/max(EF.ret.all$Mean)
    EF.ret.all <- EF.ret.all[,c(1,52)]
    colnames(EF.ret.all) <- c("Risk","MeanRet")
    EF.out.all$Mean <- apply(EF.out.all[,-1],1,mean)
    EF.sum <- data.frame(Spp = rep(c("Fd","Lw","Pl","Py","Sx"),25), 
                         Risk = EF.out.all$EF.out.all*2, Weight = EF.out.all$Mean)
    
    ##graph efficient frontier weights
    assign(paste("m",x,"EleHill", sep = ""), ggplot(EF.sum)+
             geom_bar(aes(x = Risk, y = Weight, fill = Spp), stat = "identity")+
             geom_line(data = EF.ret.all, aes(x = Risk, y = MeanRet), size = 2)+
             geom_vline(xintercept = 12)+
             labs(x = "Risk (0 = High Risk, 50 = Low Risk)")+
             ggtitle((fires[x])))

    ###Create violin plots
    rownames(portOutput) <- portOutput$Species
    portOutput <- portOutput[,-1]
    portOut <- t(portOutput)
    portOut <- as.data.frame(portOut)
    
    vioplot(portOut$Fd,portOut$Lw,portOut$Pl,portOut$Py,portOut$Sx,names = colnames(portOut), col = "purple")
    title(main = mods[a])
    
  }
#}
dev.off()

##plot frontiers
require(gridExtra)
ordLay <- rbind(c(1,2),
                c(3,NA))
grid.arrange(m1f30,m2f30,m3f30, layout_matrix=ordLay)




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
 