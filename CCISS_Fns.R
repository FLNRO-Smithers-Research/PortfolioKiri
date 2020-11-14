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

edatopicOverlap <- function(dat,Edatope){
  SS <- Edatope[is.na(Special),.(BGC,SS_NoSpace,Edatopic)]
  SS <- unique(SS)
  SS <- SS[complete.cases(SS),]
  BGC <- unique(dat)
  BGC <- BGC[complete.cases(BGC),]
  SSsp <- Edatope[!is.na(Codes),.(BGC,SS_NoSpace,Codes)]
  SSsp <- unique(SSsp)
  
  ##Special edatopes
  CurrBGC <- SSsp[BGC, on = "BGC", allow.cartesian = T]
  setkey(BGC, BGC.pred)
  setkey(SSsp, BGC)
  FutBGC <- SSsp[BGC, allow.cartesian = T]
  setnames(FutBGC, old = c("BGC","SS_NoSpace","i.BGC"), 
           new = c("BGC.pred","SS.pred","BGC"))
  FutBGC <- FutBGC[!is.na(SS.pred),]
  setkey(FutBGC, SiteNo, FuturePeriod, BGC,BGC.pred, Codes)
  setkey(CurrBGC,SiteNo,FuturePeriod, BGC,BGC.pred, Codes)
  new <- CurrBGC[FutBGC]
  SSsp.out <- new[,.(allOverlap = 1/.N,SS.pred,BGC.prop), keyby = .(SiteNo,FuturePeriod,BGC,BGC.pred,SS_NoSpace)]
  
  ##regular
  CurrBGC <- SS[BGC, on = "BGC", allow.cartesian = T]
  setkey(BGC, BGC.pred)
  setkey(SS, BGC)
  FutBGC <- SS[BGC, allow.cartesian = T]
  setnames(FutBGC, old = c("BGC","SS_NoSpace","i.BGC"), 
           new = c("BGC.pred","SS.pred","BGC"))
  
  setkey(FutBGC, SiteNo, FuturePeriod, BGC,BGC.pred, Edatopic)
  FutBGC[,BGC.prop := NULL]
  setkey(CurrBGC,SiteNo,FuturePeriod, BGC,BGC.pred, Edatopic)
  new <- merge(CurrBGC,FutBGC, all = T)
  setkey(new, SiteNo,FuturePeriod,BGC,BGC.pred,SS_NoSpace,SS.pred)
  ##new <- new[complete.cases(new),]
  
  ###forwards overlap
  SS.out <- new[,.(SS.prob = .N), 
                keyby = .(SiteNo,FuturePeriod,BGC,BGC.pred,SS_NoSpace,SS.pred)]
  SS.out2 <- new[,.(SS.Curr = length(unique(Edatopic)), BGC = unique(BGC.prop)), 
                 keyby = .(SiteNo,FuturePeriod,BGC,BGC.pred,SS_NoSpace)]
  comb <- SS.out2[SS.out]
  
  ###reverse overlap
  SS.out.rev <- new[,.(SS.prob = .N), 
                    keyby = .(SiteNo,FuturePeriod,BGC,BGC.pred,SS.pred,SS_NoSpace)]
  SS.out2.rev <- new[,.(SS.Curr = length(unique(Edatopic)), BGC = unique(BGC.prop)), 
                     keyby = .(SiteNo,FuturePeriod,BGC,BGC.pred,SS.pred)]
  SS.out2.rev <- SS.out2.rev[complete.cases(SS.out2.rev),]
  combRev <- SS.out2.rev[SS.out.rev]
  
  ##combine them
  comb[,SSProb := SS.prob/SS.Curr]
  combRev[,SSProbRev := SS.prob/SS.Curr]
  combAll <- merge(comb,combRev,by = c("SiteNo","FuturePeriod","BGC","BGC.pred","SS_NoSpace","SS.pred"))
  combAll[,allOverlap := SSProb*SSProbRev]
  setnames(combAll, old = "BGC.1.x",new = "BGC.prop")
  combAll <- combAll[,.(SiteNo, FuturePeriod, BGC, BGC.pred, SS_NoSpace, 
                        allOverlap, SS.pred, BGC.prop)]
  combAll <- rbind(combAll, SSsp.out)
  combAll <- combAll[!is.na(SS_NoSpace),]
  
  ##add in BGC probability
  combAll <- combAll[complete.cases(combAll),]
  combAll[,SSratio := allOverlap/sum(allOverlap), by = .(SiteNo, FuturePeriod, BGC, BGC.pred,SS_NoSpace)] ##should check this?
  setorder(combAll, SiteNo, FuturePeriod, BGC, BGC.pred, SS_NoSpace)
  
  combAll <- unique(combAll)
  setkey(combAll, SiteNo, FuturePeriod, BGC,BGC.pred)
  temp <- unique(combAll[,.(SiteNo,FuturePeriod,BGC,BGC.pred,BGC.prop)])
  temp[,BGC.prop := BGC.prop/sum(BGC.prop), by = .(SiteNo,FuturePeriod,BGC)]
  temp <- unique(temp)
  combAll[,BGC.prop := NULL]
  combAll <- temp[combAll]
  combAll[,SSprob := SSratio*BGC.prop]
  
  return(combAll)
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
  Y1=Y1[,..varList]
  
  Y1$BGC <- as.factor(Y1$BGC)
  
  ##Predict future subzones######
  Y1$BGC.pred <- predict(BGCmodel, Y1[,-c(1:3)])[['predictions']]
  
  Y1 <- separate(Y1, Model, into = c("GCM","Scenario","FuturePeriod"), sep = "_", remove = T)
  Y1$FuturePeriod <- gsub(".gcm","",Y1$FuturePeriod)
  Y1 <- Y1[,c("GCM","Scenario","FuturePeriod","SiteNo","BGC","BGC.pred")]
  
  
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
  
  Y3.sub1 <- as.data.table(Y3.sub1)
  Y3.sub1[,`:=`(BGC.len = NULL, Pred.len = NULL)]
  
  SiteNo.suit <- edatopicOverlap(Y3.sub1, E1)
  
  return(list(SiteNo.suit[!is.na(SiteNo.suit$SSprob),],Pred.len))
}

