##This script uses a Monte Carlo approach to model tree growth wiht SI numbers and predicted suitability 
##It then uses the Markowitz portfolio method to calculate the optimal mix of species.
##Kiri Daust, June 2018

rm(list=ls())
.libPaths("E:/R packages351")

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
require(doParallel)
require(corrplot)

wd=tk_choose.dir(); setwd(wd)

###OPTIONAL: Set up to run loops in parallel###
require(doParallel)
set.seed(123321)
coreNum <- as.numeric(detectCores()-2)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)
clusterEvalQ(coreNo, .libPaths("E:/R packages351"))



nsuitF <- read.csv ("TreeSppSuit_v10.10BulkleyTest.csv")
#nsuitF <- read.csv ("TreeSppSuit_v10.10.csv")

SuitTable <- unique(nsuitF)

colnames(SuitTable)[2] <- "SS_NoSpace"

SIBEC <- read.csv("PredSIforPort_Sept20.csv", stringsAsFactors = FALSE)
#SIBEC <- read.csv("PredSIforPort_Sept20_Testing.csv", stringsAsFactors = FALSE)

SIBEC <- SIBEC[,-5]
colnames(SIBEC)[c(1,3,5)] <- c("SS_NoSpace", "MeanPlotSiteIndex","TreeSpp")

SSPredAll <- read.csv("BulkleyTSA_SSpredicted.csv", stringsAsFactors = FALSE) ##Import SS predictions from CCISS tool: must have columns MergedBGC, Source, SS_NoSpace, SSprob, SSCurrent, FuturePeriod, SiteNo
#SSPredAll <- read.csv("BulkleyTSA_SSpredictedNOCHANGE.csv", stringsAsFactors = FALSE) ##Import SS predictions from CCISS tool: must have columns MergedBGC, Source, SS_NoSpace, SSprob, SSCurrent, FuturePeriod, SiteNo
##############Select Site Series to Run
selectBGC <- select.list(choices = sort(unique(SSPredAll$SSCurrent)), graphics = TRUE) ###Select BGC to run for
SSPredAll <- SSPredAll[SSPredAll$SSCurrent == selectBGC,]
sites <- as.numeric(as.character(unique(SSPredAll$SiteNo)))

##Build Correlation matrix for species of interest
treeSuitraw <- read.csv ("TreeSppSuit_v10.10BulkleyTest.csv")# some suitability adjustments
treeSuitraw <- treeSuitraw[treeSuitraw$Unit %in% SSPredAll$SS_NoSpace,]
treeSuitraw$Suitability <- ifelse(treeSuitraw$Suitability %in% c('3'), 3, 
                                  ifelse(treeSuitraw$Suitability %in% c('2'), 2, 
                                         ifelse(treeSuitraw$Suitability %in% c('1'), 1,treeSuitraw$Suitability )))
treeSuitraw$All <- paste(treeSuitraw$Unit, treeSuitraw$Spp)
treeSuitremove <- treeSuitraw[duplicated(treeSuitraw$All),]# remove duplicate entries 
treeSuitnew <-treeSuitraw[!duplicated(treeSuitraw$All),]
treeSuit <- treeSuitraw[,-c(5)]
#Trees <- c("Bl","Pl","Sx")#SBSmc2 and ESSFmc current
#Trees <- c("Bl","Cw","Fd","Lw","Pl","Sx")#SBSmc2 futures "Hw",
#Trees <- c("Ba", "Bl","Hw","Pl","Sx") #ICHmc1Current
Trees <- c("Ba", "Bl","Cw","Fd","Hw","Lw", "Pl","Py","Sx","Ss") #
treeSuit <- treeSuit[treeSuit$Spp %in% Trees,]
treeSuitMatrix <- dcast(treeSuit, Unit ~ Spp, mean)
treeSuitMatrix[is.na(treeSuitMatrix)] <- 0
#sigma <- treeSuitMatrix
rownames(treeSuitMatrix) <- treeSuitMatrix[,1]
treeSuitMatrix <- treeSuitMatrix[,-1]
##create correlation matrix
sigma <- cor(treeSuitMatrix, method=c("spearman"))
write.csv(sigma, "CorrelationMatrix1.csv")
##Graphical output1
corrplot(sigma,type = "lower", order = "hclust", tl.col = "black", tl.srt = 45, tl.cex = .5, is.corr = TRUE)
##Graphical output correlations limited species
chart.Correlation(sigma, histogram = TRUE, pch = 19)


############USE THIS ONLY IF NOT GENERATING ON THE FLY
#sigma <- read.csv("CorrelationMatrix1_modified.csv")
#rownames(sigma) <- sigma[,1]
#sigma <- sigma[,-1]
###########
nSpp <- length(Trees)
treeList <- Trees
sigma <- sigma[Trees,Trees]
sigma <- as.matrix(sigma)
momentargs <- list()
momentargs$sigma <- sigma ##set up to use in PortfolioAnalystics 
#####Randomly select 100 sites############

#SSPredAll <- SSPredAll[SSPredAll$SiteNo %in% sample(sites, 10, replace = FALSE),]
#SSPredAll <- SSPredAll[!is.na(SSPredAll$SSprob),]
########################

combineList <- function(...) {##Combine multiple dataframe in foreach loop
  mapply(FUN = rbind, ..., SIMPLIFY=FALSE)
}

###foreach site
allSites <- foreach(SNum = unique(SSPredAll$SiteNo), .combine = combineList, .packages = c("foreach","reshape2","dplyr","magrittr","PortfolioAnalytics")) %dopar% {
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
    
    ###check that there aren't errors in the table
    temp <- aggregate(SS_NoSpace ~ TreeSpp, current, FUN = length)
    if(any(temp$SS_NoSpace > 1)){
      stop("There are partial duplicates in the suitablity table. Please fix them. :)")
    }
    
    #current$Suitability[current$TreeSpp == "Fd"] <- 3
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
    plot(0,0, type = "n", xlim = c(1,100), ylim = c(0,3000), xlab = "Year", ylab = "Volume")
    
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
     # init.portfolio <- add.constraint(portfolio = init.portfolio, type = "box", min = c(0,0,0,0,0,0,0,0,0,0), max = c(1,1,1,1,0,1,1,1,1,1)) ##set min and max weight for each species
      init.portfolio <- add.constraint(portfolio = init.portfolio, type = "box", min = rep(0, nSpp),max = rep(1, nSpp))  ##set min and max weight for each species
      
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
      qu <- add.objective(portfolio=qu, type="risk", name="var", risk_aversion = 25) ###set risk_aversion
      
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

violin <-ggplot(portOut)+
  geom_violin(aes(x = variable, y = value), draw_quantiles = c(0.25,0.5,0.75), scale = "width")+
  labs(x = "Spp", y = "Weight")
selectBGC <- gsub("/", "_", selectBGC)
name1 = paste(selectBGC,"_Bulkley_TSA_CC_Violin",".pdf", sep = "")
plot(violin)
   ggsave(name1, violin, device="pdf")
dev.off ()
####Efficient Frontier
EF.sum <- allSites$Frontier
EF.sum <- aggregate(Weight ~ Spp + Risk, EF.sum, FUN = mean)
EF.ret.all <- allSites$Return
EF.ret.all <- aggregate(MeanRet ~ Risk, EF.ret.all, FUN = mean)
write.csv (EF.sum, paste(selectBGC,"_CC_Portfolio_Weights",".csv", sep = ""))
portfolio <- ggplot(EF.sum)+
  geom_bar(aes(x = Risk, y = Weight, fill = Spp), stat = "identity")+
  geom_line(data = EF.ret.all, aes(x = Risk, y = MeanRet), size = 2)+
  geom_vline(xintercept = 25)+
  labs(x = "Risk (0 = High Risk, 50 = Low Risk)", y= "% of species in portfolio")+
  scale_fill_brewer(palette = "Paired")+
  ggtitle(selectBGC, "CC_Portfolio_NoHw")
plot(portfolio)
selectBGC <- gsub("/", "_", selectBGC)
name2 = paste(selectBGC,"Bulkley_TSA_CC",".pdf", sep = "")
ggsave(name2, portfolio, device="pdf")
dev.off ()
 
