BGC <- as.data.table(Y3.sub1)
SS <- as.data.table(Edatope)
SS <- SS[,.(BGC,SS_NoSpace,Edatopic)]
SS <- unique(SS)
BGC <- unique(BGC)
BGC <- BGC[,-c("BGC.len","Pred.len")]

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
SS.out <- new[,.(SS.prob = .N), 
              keyby = .(SiteNo,FuturePeriod,BGC,BGC.pred,SS_NoSpace,SS.pred)]
SS.out2 <- new[,.(SS.Curr = .N, BGC = unique(BGC.prop)), keyby = .(SiteNo,FuturePeriod,BGC,BGC.pred,SS_NoSpace)]
comb <- SS.out2[SS.out]
comb[,SSProb := (SS.prob/SS.Curr)*BGC.1]
