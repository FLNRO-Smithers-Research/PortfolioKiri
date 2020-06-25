require(data.table)
require(doBy)
require(randomForest)
require(foreach)
require(dplyr)
require(reshape2)
library(doParallel)
library(tidyr)
require(ranger)

addVars <- function(dat){
  dat$PPT_MJ <- dat$PPT05 + dat$PPT06  # MaY/June precip
  dat$PPT_JAS <- dat$PPT07 + dat$PPT08 + dat$PPT09  # July/Aug/Sept precip
  dat$PPT.dormant <- dat$PPT_at + dat$PPT_wt  # for calculating spring deficit
  dat$CMD.def <- 500 - (dat$PPT.dormant)  # start of growing season deficit original value was 400 but 500 seems better
  dat$CMD.def[dat$CMD.def < 0] <- 0  #negative values set to zero = no deficit
  dat$CMDMax <- dat$CMD07
  dat$CMD.total <- dat$CMD.def + dat$CMD
  dat$CMD.grow <- dat$CMD05 + dat$CMD06 +dat$CMD07 +dat$CMD08 +dat$CMD09
  dat$DD5.grow <- dat$DD5_05 + dat$DD5_06 + dat$DD5_07 + dat$DD5_08 + dat$DD5_09
  dat$CMDMax <- dat$CMD07 # add in so not removed below
  dat$DDgood <- dat$DD5 - dat$DD18
  dat$DDnew <- (dat$DD5_05 + dat$DD5_06 +dat$DD5_07  + dat$DD5_08)  - (dat$DD18_05 + dat$DD18_06 +dat$DD18_07 +dat$DD18_08)
  dat$TmaxJuly <- dat$Tmax07
  return(dat)
}

CCISS_cbst <- function(Y1,BGCmodel){
  Y1 <- addVars(Y1)
  
  vars <- BGCmodel[["forest"]][["independent.variable.names"]]
  varList = c("Model", "SiteNo", "BGC", vars)
  colnames (Y1) [1:3] = c("Model", "SiteNo", "BGC")
  Y1=Y1[,varList]
  
  Y1$BGC <- as.factor(Y1$BGC)
  
  ##Predict future subzones######
  Y1$BGC.pred <- predict(BGCmodel, Y1[,-c(1:3)])[['predictions']]
  
  Y1 <- separate(Y1, Model, into = c("Model","Scenario","FuturePeriod"), sep = "_", remove = T)
  Y1$FuturePeriod <- gsub(".gcm","",Y1$FuturePeriod)
  CBSTdat <- Y1[,c("Model","Scenario","FuturePeriod","SiteNo","BGC","BGC.pred")]
  CBSTdat$Model <- paste(CBSTdat$Model,CBSTdat$Scenario,sep = "_")
  CBSTdat <- CBSTdat[,-2]
  colnames(CBSTdat)[1] <- "GCM"
  return(CBSTdat)
}

CCISS_Spp <- function(Y1,BGCmodel,E1){
  Y1 <- addVars(Y1)
  
  vars <- BGCmodel[["forest"]][["independent.variable.names"]]
  varList = c("Model", "SiteNo", "BGC", vars)
  colnames (Y1) [1:3] = c("Model", "SiteNo", "BGC")
  Y1=Y1[,varList]
  
  Y1$BGC <- as.factor(Y1$BGC)
  
  ##Predict future subzones######
  Y1$BGC.pred <- predict(BGCmodel, Y1[,-c(1:3)])[['predictions']]
  
  Y1 <- separate(Y1, Model, into = c("GCM","Scenario","FuturePeriod"), sep = "_", remove = T)
  Y1$FuturePeriod <- gsub(".gcm","",Y1$FuturePeriod)
  Y1 <- Y1[,c("GCM","Scenario","FuturePeriod","SiteNo","BGC","BGC.pred")]
  
  #E1 <-fread("./InputsGit/Edatopic_v11_7.csv")
  
  E1 <- E1[,-c(7:12)]
  E1 <- unique(E1)
  
  ###create list of focal BGCs & edatopic space
  e1 <- unique(Y1$BGC)
  edatopic1 <- E1[E1$BGC %in% e1,]
  edatopic1$Codes[edatopic1$Codes == ""] <- NA
  
  ###create list of predicted BGCs and edatopic space
  e2 <- unique(Y1$BGC.pred)
  edatopic2 <- E1[E1$BGC %in% e2,]
  edatopic2$Codes[edatopic2$Codes == ""] <- NA
  
  ##calculate BGC proportion
  Y3.sub1 <-Y1
  Y3.sub1 <- Y3.sub1[,c("SiteNo","FuturePeriod","BGC","BGC.pred")]
  Y3.sub1 <- Y3.sub1[order(Y3.sub1$SiteNo, Y3.sub1$FuturePeriod, Y3.sub1$BGC,Y3.sub1$BGC.pred),]
  Y3.sub1 <- mutate(Y3.sub1, BGC = as.character(BGC),BGC.pred = as.character(BGC.pred))
  Y3.sub1$BGC.len <- ave(Y3.sub1$BGC, Y3.sub1$SiteNo, Y3.sub1$FuturePeriod, Y3.sub1$BGC, FUN = length)
  Pred.len <- aggregate(BGC.len ~ SiteNo + FuturePeriod + BGC + BGC.pred, Y3.sub1, FUN = length)
  colnames(Pred.len)[5] <- "Pred.len"
  Y3.sub1 <- merge(Y3.sub1, Pred.len, by = c("SiteNo","FuturePeriod","BGC","BGC.pred"), all = TRUE)
  Y3.sub1$BGC.len <- as.numeric(Y3.sub1$BGC.len); Y3.sub1$Pred.len <- as.numeric(Y3.sub1$Pred.len)
  Y3.sub1$BGC.prop <- Y3.sub1$Pred.len/Y3.sub1$BGC.len
  Y3.sub1 <- Y3.sub1[order(Y3.sub1$SiteNo,Y3.sub1$FuturePeriod,Y3.sub1$BGC,Y3.sub1$BGC.pred),]
  
  
  BGClist = unique(Y3.sub1$BGC)
  FuturePeriod.list <- unique(Y3.sub1$FuturePeriod)
  BGCfutures.list <- unique(Y3.sub1$BGC.pred) ### to use later on to limit the site series
  BGCfocalE <- edatopic1[edatopic1$BGC %in% Y3.sub1$BGC  ,] ### extracts edatopic space for BGC focal of SiteNo
  BGCfutureE <- edatopic2[edatopic2$BGC %in% Y3.sub1$BGC.pred  ,] #extracts edatopic info only for predicted BGCs
  Y3.sub1$SiteNo <- as.character(Y3.sub1$SiteNo)
  SiteNo.list = unique(Y3.sub1$SiteNo)
  Y3.sub1$BGC <- as.character(Y3.sub1$BGC)
  Y3.sub1$BGC.pred <- as.character(Y3.sub1$BGC.pred)
  
  gc()
  
  coreNum <- as.numeric(detectCores()-1)
  coreNo <- makeCluster(coreNum)
  registerDoParallel(coreNo, cores = coreNum)
  
  edaOverlap <- function(dat1, dat2){
    return(length(dat1[dat1 %in% dat2])/length(dat2))
  }
  
  SiteNo.suit <-  foreach(SNL = SiteNo.list, .combine = rbind, .packages = c("doBy","foreach")) %dopar% {##  for each SiteNo in the data
    
    SiteFuture.suit <- foreach(i = FuturePeriod.list, .combine = rbind)  %do% {
      
      Y3.each <- Y3.sub1[Y3.sub1$SiteNo %in% SNL ,] ## extracts data for each site
      Y3.each <- Y3.each[Y3.each$FuturePeriod %in% i,] ##extracts data for each time period
      Y3.BGC <- unique(Y3.each$BGC)
      Y3.BGC.pred<- unique(Y3.each$BGC.pred)
      BGCfocalE <- edatopic1[edatopic1$BGC %in% Y3.BGC  , ] ### extracts edatopic space for BGC focal of SiteNo
      BGCfutureE <- edatopic2[edatopic2$BGC %in% Y3.BGC.pred  , ] #extracts edatopic info only for predicted BGCs
      
      Y3.SSlist = as.list(unique(BGCfocalE$SS_NoSpace))
      
      FTS2 <-  foreach(SS = Y3.SSlist, .combine =rbind, .packages = c("doBy","foreach")) %do% {      ##  for each site series for a SiteNo BGC
        
        SSfocal <- BGCfocalE[BGCfocalE$SS_NoSpace %in% SS ,] ###find focal site series cells
        SSfocalE <- as.list(unique(SSfocal$Edatopic))
        
        ##select site series only with some edatopic overlap with SSfocal
        
        SSfutureE <- BGCfutureE[BGCfutureE$Edatopic %in% SSfocalE,]
        futureZones <- unique(SSfutureE$BGC)
        futureSS.names <- unique(SSfutureE$SS_NoSpace)
        if(length(SSfutureE$Edatopic) > 0){
          
          SSfocal <- BGCfocalE[BGCfocalE$SS_NoSpace %in% SS,]
          SSfuture <- BGCfutureE[BGCfutureE$SS_NoSpace %in% futureSS.names,]
          
          ##match site series within each projected subzone
          futureSS <- foreach(futSS = futureZones, .combine = rbind) %do% {
            dat <- SSfuture[SSfuture$BGC == futSS,]
            fut <- dat[dat$Edatopic %in% SSfocal$Edatopic,]
            if(any(fut$Codes[!is.na(fut$Codes)] %in% SSfocal$Codes[!is.na(SSfocal$Codes)]) & (length(fut$Codes[!is.na(fut$Codes)])/length(fut$Edatopic)) > 0.1){###Are there matchin special edatopic cells?
              oldSp <- unique(SSfocal$Codes[!is.na(SSfocal$Codes)])
              newSp <- unique(fut$SS_NoSpace[match(oldSp, fut$Codes)])
              dat <- unique(dat[dat$SS_NoSpace %in% newSp, -c(4:5)]) #which ones have the new special edatope
              dat$alloverlap <- 1/length(dat$SS_NoSpace)
              
            }else{
              if(all(is.na(SSfocal$Codes))){
                dat <- dat[is.na(dat$Codes),] ####remove special edatopes
              }
              dat <- foreach(x = unique(as.character(dat$SS_NoSpace)), .combine = rbind) %do% {
                dat1 <- dat[dat$SS_NoSpace == x,]
                Overlap <- edaOverlap(SSfocal$Edatopic, dat1$Edatopic)
                Revoverlap <- edaOverlap(dat1$Edatopic, SSfocal$Edatopic)
                dat1 <- unique(dat1[-c(4:5)])
                dat1$alloverlap <- Overlap*Revoverlap
                dat1 <- as.data.frame(dat1)
              }
              ##if not special, loop through each SS to calculate overlap
              
            }
            dat
          }
          
          if(length(futureSS) > 0){
            futureSS <- futureSS[futureSS$alloverlap > 0,] ###Adjust this to remove low overlap SS
            Y3.each <- unique(Y3.each)
            Y3.each <- Y3.each[Y3.each$BGC.pred %in% futureSS$BGC,] ###Remove predictions with no overlap and recalculate
            Y3.each$BGC.prop <- Y3.each$BGC.prop/sum(Y3.each$BGC.prop)
            futureSS <- merge(futureSS, Y3.each[,c("BGC.pred","BGC.prop")], by.x = "BGC", by.y = "BGC.pred", all.x = TRUE)
            
            
            ####Calculate the SS ratio
            SSoverlap <- summaryBy(alloverlap~BGC, data=futureSS, id = 'SS_NoSpace', FUN=c(sum))
            futureSS$overlaptot<- SSoverlap$alloverlap.sum[match(futureSS$BGC, SSoverlap$BGC )]
            futureSS$SSratio <- futureSS$alloverlap/futureSS$overlaptot
    
            ####Calculated the overall site series probability
            futureSS$SSprob <- (futureSS$BGC.prop * futureSS$SSratio)
            sum(futureSS$SSprob)
            
            futureSS$SSCurrent <- rep(SS,length(futureSS$SSprob))
            futureSS <- futureSS[,-c(4:7)]
            
            futureSS$FuturePeriod <- as.character(i)  
            futureSS$SiteNo <- as.character(SNL)
            
            futureSS <- as.data.frame(futureSS)
          }else{
            warning("No matching edatopes - some predictions have been removed")
          }
        }
        
      } #For each Site
    } #For each year
    
  } # for all
  stopCluster(coreNo)
  return(list(SiteNo.suit[!is.na(SiteNo.suit$SSprob),],Pred.len))
}

