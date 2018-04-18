library(plyr)
library(mgcv)
library(tidyr)
datadir<-'C:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/haddock'
datadir<-'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/data/haddock'
setwd(datadir)

#READ IN AND STACK DATA: SHORT FROMAT
#LENGTH AT AGE
ldat<-read.csv('length_atage_RV_2003_assess.csv',header=TRUE)
ldat<-ldat %>% gather(year,len,X1970:X2003)
ldat$year<-as.numeric(substr(ldat$year,2,5))
ldat<-subset(ldat,is.na(Age)==FALSE)
names(ldat)<-tolower(names(ldat))

#NUMBERS AT AGE
ndat<-read.csv('numbers_atage_VPA_2011_assess.csv',header=TRUE)
ndat<-ndat %>% gather(year,n,X1970:X2010)
ndat$year<-as.numeric(substr(ndat$year,2,5))
names(ndat)<-tolower(names(ndat))

#WEIGHT AT AGE
wdat<-read.csv('weight_atage_RV_2003_assess.csv',header=TRUE)
wdat<-wdat %>% gather(year,wgt,X1970:X2003)
wdat$year<-as.numeric(substr(wdat$year,2,5))

wldat<-merge(wdat,ldat,by=c('year','age'),all=TRUE)

#STANDARDIZED LENGTH BINS TO USE THROUGHOUT
len<-seq(min(wldat$len,na.rm=TRUE),max(wldat$len,na.rm=TRUE),1)



#CALCULATE ABUNDANCE AT LENGTH FROM ABUNDANCE AT AGE AND LENGTH AT AGE DATA
f<-function(d){
print(unique(d$year))
#MODEL TO PREDICT AGE FROM STANDARDIZED LENGTHS
mod<-gam(age~s(len,k=5),data=d,gamma=1.4)
pdat<-data.frame(len=len)
p<-predict(mod,newdata=pdat)
pdat$page<-p
#plot(d$len,d$age,pch=15)
#lines(pdat$len,pdat$page)

#GET NUMBERS AT AGE
d2<-subset(ndat,year==unique(d$year))
d2<-merge(d2,d,by=c('year','age'),all.x=TRUE,all.y=FALSE)
#MODEL OF NUMBER AS A FUNCTION OF LENGTH
mod2<-gam(n~s(len,k=10),data=d2,gamma=.8)
#USE THE INITIAL MODEL TO PREDICT THE NUMBER AS A FUNCTION OF ALL LENGTHS
pdat$pnum<-predict(mod2,newdata=pdat)
#PREDICTIONS ARE NEGATIVE, SET TO 0
pdat$pnum<-ifelse(pdat$pnum<0,0,pdat$pnum)
#GET AS PERCENTAGES
pdat$pnum<-pdat$pnum/sum(pdat$pnum)
#NOW MULTIPLY BY TOTAL OBSERVED TO GET ON CORRECT SCALE
pdat$pnum<-pdat$pnum*sum(d2$n,na.rm=TRUE)

return(pdat)
}
NL<-ddply(ldat,.(year),.fun=f,.progress='text')


#CALCULATE DAILY RATION FROM LENGTH-WEIGHT AND PUBLISHED RATION EQN
f<-function(d){
print(unique(d$year))
mod<-gam(wgt~s(len,k=4),data=d,gamma=1.4)
pdat<-data.frame(len=len)
p<-predict(mod,newdata=pdat)
pdat$pwgt<-p
print(length(unique(pdat$len)))

#PREDICTED DAILY RATION FROM PUBLISHED LIT.
pdat$dr<-exp((-0.36*log(pdat$pwgt))-.101)
pdat$drtot<-pdat$dr*pdat$pwgt
#plot(pdat$len,pdat$drtot)
#plot(log(pdat$pwgt),pdat$dr)
return(pdat)
}
PDR<-ddply(wldat,.(year),.fun=f,.progress='text')

#NOW SUM THE PRODUCT OF PDR AND NL FOR EACH YEAR TO GET PREDATION INTENSITY
dat<-merge(PDR,NL,by=c('year','len'),all=TRUE)
dat$pi.nl<-dat$drtot*dat$pnum

dat<-data.frame(year=sort(unique(dat$year)),
           pi=tapply(dat$pi.nl,dat$year,function(x) sum(x,na.rm=TRUE)))
plot(dat$year,dat$pi,pch=15,type='b')
names(dat)[2]<-c('had.pi')

setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
setwd(datadir)
save(dat,file='haddock_predation_index.RData')

#EXPLORE CORRELATION TO HADDOCK SSB: HIGH (R=0.75)
tdat<-data.frame(year=sort(unique(ndat$year)),
                 n=tapply(ndat$n,ndat$year,function(x) sum(x,na.rm=TRUE)))
plot(data$year,data$pi,pch=15,type='b',xlim=c(1970,2010))
par(new=TRUE)
plot(tdat$year,tdat$n,col='red',type='b',xlim=c(1970,2010))
a<-merge(data,tdat,by=c('year'),all=TRUE)
plot(a$pi,a$n,pch=15)
cor(a$pi,a$n,use='pairwise.complete.obs')
