BGC <- as.data.table(Y3.sub1)
Edatope <- as.data.table(Edatope)
SS <- as.data.table(Edatope)
SS <- SS[is.na(Special),.(BGC,SS_NoSpace,Edatopic)]
SS <- unique(SS)
BGC <- unique(BGC)
BGC <- BGC[,-c("BGC.len","Pred.len")]
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
setkey(CurrBGC,SiteNo,FuturePeriod, BGC,BGC.pred, Edatopic)
new <- CurrBGC[FutBGC]
setkey(new, SiteNo,FuturePeriod,BGC,BGC.pred,SS_NoSpace,SS.pred)
new <- new[complete.cases(new),]

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
combRev <- SS.out2.rev[SS.out.rev]

comb[,SSProb := SS.prob/SS.Curr]
combRev[,SSProbRev := SS.prob/SS.Curr]
combAll <- merge(comb,combRev,by = c("SiteNo","FuturePeriod","BGC","BGC.pred","SS_NoSpace","SS.pred"))
combAll[,allOverlap := SSProb*SSProbRev]
setnames(combAll, old = "BGC.1.x",new = "BGC.prop")
combAll <- combAll[,.(SiteNo, FuturePeriod, BGC, BGC.pred, SS_NoSpace, 
                      allOverlap, SS.pred, BGC.prop)]
combAll <- rbind(combAll, SSsp.out)

##maybe need to recalculate BGC prob?
combAll[,SSratio := allOverlap/sum(allOverlap), by = .(SiteNo, FuturePeriod, BGC, BGC.pred,SS_NoSpace)]
temp <- unique(combAll[,.(SiteNo,FuturePeriod,BGC,BGC.pred,BGC.prop)])
temp[,BGC.prop := BGC.prop/sum(BGC.prop), by = .(SiteNo,FuturePeriod,BGC)]
combAll[,BGC.prop := NULL]
combAll <- temp[combAll, on = c("SiteNo","FuturePeriod","BGC","BGC.pred")]
combAll[,SSprob := SSratio*BGC.prop]
