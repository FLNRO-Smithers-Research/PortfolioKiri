require(dplyr)
require(ggplot2)
require(MASS)
require(magrittr)
require(foreach)
require(reshape2)
require(reticulate)
require(Rcpp)
library(gridExtra)
library(data.table)
library(scales)
library(tidyr)
library(magrittr)
library(ggthemes)
library(flextable)
require(rgdal)
require(tmap)
library(sf)
require(tufte)
require(Hmisc)
require(janitor)
require(tidyverse)
require(Rmisc)
require(bookdown)


cloud_dir <- "C:/Users/kirid/Sync/CCISS_data/"
climDat <- fread(paste0(cloud_dir,"/WNA_4k_HexPts_BGC_Normal_1961_1990MSY.csv"))
climDat <- climDat[Eref_sm > -1000,]
climDat <- climDat[,c("BGC","CMD","FFP","Eref_sm","MCMT","SHM","TD","PAS","DD5_sp","MSP")]
colnames(climDat)[1] <- "BGC"
climAve <- climDat[,lapply(.SD, mean), by = "BGC"]

Trees <- c("Ba", "Bl","Cw","Dr", "Fdi","Hw","Lw","Pli","Py","Sx") ##set species to use in portfolio
dir <- "CBSTNoMigration"
files <- list.files(dir)
allSppDat <- foreach(Spp = Trees, .combine = rbind) %do% {
  fName <- paste(dir, files[grep(paste(Spp,"_",sep = ""),files)], sep = "/")
  temp <- fread(fName)
  temp$Spp <- Spp
  temp
}
colnames(allSppDat) <- c("Site","Seed","Height", "Spp")
allSppDat[Spp == "Fdi",Spp := "Fd"][Spp == "Pli",Spp := "Pl"]
Units <- unique(allSppDat$Site)

library(StatMatch)
test <- mahalanobis.dist(climAve[,-1])

Trees <- c("Ba", "Bl","Cw","Dr", "Fd","Hw","Lw","Pl","Py","Sx")
out <- foreach(SppCurr = Trees, .combine = rbind) %do% {
  cat("Processing",SppCurr,"\n")
  SppDat <- allSppDat[Spp == SppCurr,]
  t1 <- foreach(SiteCurr = Units, .combine = rbind) %do% {
    cat("Processing",SiteCurr,"\n")
    dat <- SppDat[Site == SiteCurr,]
    trainDat <- climAve[dat, on = c(BGC = "Seed")]
    trainDat <- trainDat[!is.na(CMD), -c("Site","Spp")]
    rfMod <- randomForest(Height ~ ., data = trainDat[,-1], importance = T)
    predDat <- climAve[!BGC %chin% trainDat$BGC,]
    predDat$Height <- predict(rfMod, newdata = predDat[,-1])
    t2 <- data.table(Site = rep(SiteCurr,nrow(predDat)),Seed = predDat$BGC, Height = predDat$Height)
    t3 <- data.table(Site = predDat$BGC, Seed = rep(SiteCurr,nrow(predDat)), Height = predDat$Height)
    rbind(t2,t3)
  }
  t1$Spp <- SppCurr
  t1
}

allSppDat <- rbind(allSppDat, out)
fwrite(allSppDat, "FilledCBSTNoMigration.csv")
fill1 <- fread("FilledCBSTNoMigration.csv")
fill1 <- fill1[!Site %chin% allSppDat$Site,]

Trees <- c("Ba", "Bl","Cw","Dr", "Fd","Hw","Lw","Pl","Py","Sx")
Units <- unique(SppDat$Site)
add <- foreach(SppCurr = Trees, .combine = rbind) %do% {
  cat("Processing",SppCurr,"\n")
  SppDat <- fill1[Spp == SppCurr,]
  t1 <- foreach(SiteCurr = Units, .combine = rbind) %do% {
    cat("Processing",SiteCurr,"\n")
    dat <- SppDat[Site == SiteCurr,]
    temp <- data.table(Site = SiteCurr, Seed = SiteCurr, Height = 1, Spp = SppCurr)
    dat <- rbind(dat,temp)
    trainDat <- climAve[dat, on = c(BGC = "Seed")]
    trainDat <- trainDat[!is.na(CMD), -c("Site","Spp")]
    rfMod <- randomForest(Height ~ ., data = trainDat[,-1], importance = T)
    predDat <- climAve[!BGC %chin% trainDat$BGC,]
    predDat$Height <- predict(rfMod, newdata = predDat[,-1])
    t2 <- data.table(Site = rep(SiteCurr,nrow(predDat)),Seed = predDat$BGC, Height = predDat$Height)
    t2 <- rbind(t2,temp[,-4])
    t2
  }
  t1$Spp <- SppCurr
  t1
}

fill2 <- fread("FilledCBSTNoMigration.csv")
fill2 <- rbind(fill2,add)
fwrite(fill2,"FilledCBSTNoMigration_v2.csv")

t1 <- dat[BECvar_site == "CWHvm1",]
trainDat <- climAve[t1, on = c(BGC = "BECvar_seed")]
trainDat <- trainDat[!is.na(CMD), -c("BECvar_site")]
library(randomForest)
rfMod <- randomForest(HTp_pred ~ ., data = trainDat[,-1], importance = T)
#trainDat$rfPred <- predict(rfMod, newdata = trainDat[,-c("BGC","HTp_pred")])
