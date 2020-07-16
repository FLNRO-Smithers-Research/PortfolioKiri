BGC <- as.data.table(Y3.sub1)
SS <- as.data.table(Edatope)
SS <- SS[is.na(Special),.(BGC,SS_NoSpace,Edatopic)]
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

##maybe need to recalculate BGC prob?
combAll[,SSratio := allOverlap/sum(allOverlap), by = .(SiteNo, FuturePeriod, BGC, BGC.pred,SS_NoSpace)]


dat <- fread(file.choose())
dat[,ID := seq(nrow(dat))]
dat <- data.table::melt(dat, id.vars = "ID")
setnames(dat, c("ID","Model","Response"))
dat3 <- dat[,.(NumMods = length(unique(Response))), by = ID]

dat2 <- dat[,.(Num = .N), by = .(ID, Response)]
dat2 <- dat2[,.(MaxNum = max(Num)), by = ID]
dat2[,PropSame := MaxNum/4]
