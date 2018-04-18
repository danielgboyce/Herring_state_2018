library(RColorBrewer)
library(segmented)
library(splines)
library(strucchange)
library(data.table)
library(psych)
library(reshape2)
library(gplots)
library(forecast)
library(cluster)
library(vegan)
library(ggplot2)
library(hybridHclust)
library(raster)
library(fields)
library(gridExtra)
library(colorRamps)
library(mapdata)
library(scales)
library(MASS)
library(mgcv)
library(maps)
library(plyr)
library(plotrix)
library(lubridate)
library(fossil)

datadir1<-'N://cluster_2017//scratch//spera//data//stagingdat'
datadir<-'N://cluster_2017//scratch//spera//data//finaldat_v2'
figsdir<-'C://Users//copepod//Documents//aalldocuments//literature//research//active//SPERA//Figures'
figsdir<-'C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures'



#LOADS PATHFINDER DATA TO EXAMINE PHENOLOGY
setwd(datadir)
load("pathfinder_sst_daily_1981_2012_spawnar.RData")
dat$cell<-gsub(' ','',paste(dat$lon,'_',dat$lat))


#ESTIMATE AVERAGE SST DURING FALL FOR EACH YEAR AND CELL
f<-function(d){
d<-subset(d,month>=8)
#mod<-gam(sst~s(year,k=4) + s(day,bs='cc',k=6),data=d,gamma=1.4)
mod<-gam(sst~year,data=d,gamma=1.4)
pdat<-data.frame(year=seq(min(d$year),max(d$year),1),
                 day=250)
p<-predict(mod,newdata=data.frame(year=pdat$year,day=pdat$day),se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit

modfac<-gam(sst~as.factor(year),data=d,gamma=1.4)
pdatfac<-data.frame(year=sort(unique(d$year)))
pf<-predict(modfac,newdata=data.frame(year=pdatfac$year),se.fit=TRUE)
pdatfac$p<-pf$fit
pdatfac$se<-pf$se.fit
s<-summary(mod)
dout<-data.frame(lon=unique(d$lon),
                 lat=unique(d$lat),
                 beta=s$p.table[2,1],
                 se=s$p.table[2,2],
                 pv=s$p.table[2,4],
                 r2=s$r.sq)
#plot(pdatfac$year,pdatfac$p)
#lines(pdat$year,pdat$p)
return(dout)
}
sst.trnd<-ddply(dat,.(cell),.fun=f,.progress='text')

map('world',xlim=c(-70,-65),ylim=c(42,46))
points(sst.trnd$lon,sst.trnd$lat,pch=16,cex=rescale(sst.trnd$beta,newrange=c(.2,3)),col=alpha(ifelse(sst.trnd$beta>0,'red','blue'),.8))



######################################################
#ESTIMATE PHENOLOGY CYCLE FOR EACH GRID CELL AND YEAR; THEN EXTRACT PROPERTIES: SEASONAL MAX, SEASONAL MIN, TIMING OF MAX/MIN, AMPLITUDE, DURATION WHEN TEMP IS >15
f<-function(d){
if(length(unique(d$month))>=10){
mod<-gam(sst~s(day,bs='cc',k=6),data=d,gamma=1.4)
pdat<-data.frame(day=seq(1,365,1))

p<-predict(mod,newdata=data.frame(day=pdat$day),se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit

#TEMPERATURES ABOVE 18
pdat0<-subset(pdat,p>=12.3)
dur12<-dim(pdat0)[1]
tim12<-min(pdat0$day)

#LOOK AT FALL ONLY
pdatfall<-subset(pdat,day>200 & day<300)
mnfall<-mean(pdatfall$p)
sdfall<-sd(pdatfall$p)
mxfall<-max(pdatfall$p)

mx<-subset(pdat,p==max(pdat$p))
mn<-subset(pdat,p==min(pdat$p))

    out<-data.frame(tim=mx$day,
                    tim2=mn$day,
                    amp=diff(range(pdat$p)),
                    mx=mx$p,
                    mn=mn$p,
                    dur12=dur12,
                    tim12=tim12,
                    mnfall=mnfall,
                    mxfall=mxfall,
                    sdfall=sdfall,
                    lon=unique(d$lon),
                    lat=unique(d$lat))
return(out)
} else NULL
}
phen<-ddply(dat,.(cell,year),.fun=f,.progress='text')

par(mfrow=c(2,3))
f<-function(d){
ylb<-names(d)[1]
names(d)[1]<-'y'
plot(d$year,d$y,las=1,pch=15,col=alpha('dodgerblue4',.3),ylab=ylb,xlab='Year')
}
f(subset(phen,select=c('mxfall','year')))
f(subset(phen,select=c('dur20','year')))
f(subset(phen,select=c('sdfall','year')))
f(subset(phen,select=c('tim','year')))
f(subset(phen,select=c('tim2','year')))
f(subset(phen,select=c('amp','year')))
f(subset(phen,select=c('mn','year')))
f(subset(phen,select=c('mx','year')))

##LINEAR TREND IN EACH CELL
dm<-data.frame(cell=sort(unique(phen$cell)),
               nyear=tapply(phen$tim,phen$cell,length))
dm<-subset(dm,nyear>=20)

f2<-function(d){

f3<-function(d2){
nm<-names(d2)[1]
names(d2)[1]<-'y'
mod<-gam(y~year,data=d2)
mod2<-gam(y~s(year),data=d2,gamma=1.4)
s<-summary(mod)
pdat<-data.frame(year=c(1985,2010))
pdat$p<-predict(mod2,newdata=pdat)
out<-data.frame(beta=s$p.table[2,1],
                 se=s$p.table[2,2],
                 pv=s$p.table[2,4],
                 r2=s$r.sq,
                delta=diff(pdat$p))
names(out)<-gsub(' ','',paste(names(out),'.',nm))
return(out)
}
dout<-cbind(f3(subset(d,select=c('mxfall','year'))),
             f3(subset(d,select=c('sdfall','year'))),
             f3(subset(d,select=c('tim','year'))),
             f3(subset(d,select=c('tim2','year'))),
             f3(subset(d,select=c('dur20','year'))),
             f3(subset(d,select=c('mn','year'))),
             f3(subset(d,select=c('mx','year'))),
             f3(subset(d,select=c('amp','year'))))
dout$lon<-unique(d$lon)
dout$lat<-unique(d$lat)
return(dout)
}
phen2<-ddply(subset(phen,cell %in% dm$cell),.(cell),.fun=f2,.progress='text')
save(phen2,file='phen2.RData')



#################################################
#ESTIMATE PHENOLOGY FOR PRE/POST 1988
brks<-seq(1980,2015,5)
dat$tbin5<-cut(dat$year,breaks=brks)
brks7<-seq(1980,2015,7)
dat$tbin7<-cut(dat$year,breaks=brks7)
dat$tbinpp<-ifelse(dat$year<1988,'pre','post')

f<-function(d){
#mod<-gam(sst~s(day,bs='cc',k=6) + s(lon,lat,k=20),data=d,gamma=1)
mod<-gam(sst~s(day,bs='cc',k=7) + s(lon,lat,k=20),data=d,gamma=1)
pdat<-data.frame(lon=-67,
                 lat=43.5,
                 day=seq(1,365,1),
                 year=mean(d$year))
 p<-predict(mod,newdata=pdat,se.fit=TRUE)
 pdat$p<-p$fit
 pdat$se<-p$se.fit
 return(pdat)
}
phen3<-ddply(dat,.(tbinpp),.fun=f,.progress='text')

cl1<-'dodgerblue4'
cl2<-'firebrick4'
cl<-data.frame(tbinpp=sort(unique(phen3$tbinpp)),
           cl=c(cl1,cl2))
phen3<-merge(phen3,cl,by=c('tbinpp'),all.x=TRUE,all.y=FALSE)
save(phen3,file='phen3.RData')



brks<-seq(1980,2015,5)
dat$tbin5<-cut(dat$year,breaks=brks)
brks7<-seq(1980,2015,7)
dat$tbin7<-cut(dat$year,breaks=brks7)
dat$tbinpp<-ifelse(dat$year<1988,'pre','post')

f<-function(d){
    return(data.frame(n=dim(d)[1],
                      nday=length(unique(d$day)),
                      nmonth=length(unique(d$month))))
}
ot<-ddply(dat,.(year),.fun=f)
ot<-subset(ot,nday>1)

f<-function(d){
mod<-gam(sst~s(day,bs='cc',k=5) + s(lon,lat,k=10),data=d,gamma=1)
pdat<-data.frame(lon=-67,
                 lat=43.5,
                 day=seq(1,365,1),
                 year=mean(d$year))
 p<-predict(mod,newdata=pdat,se.fit=TRUE)
 pdat$p<-p$fit
 pdat$se<-p$se.fit
 return(pdat)
}
pdat<-ddply(subset(dat,year %in% ot$year & year>1981),.(year),.fun=f,.progress='text')

f<-function(d){
    return(data.frame(tim=subset(d,p==max(d$p)[1])$day,
                      amp=max(d$p)-min(d$p)))
}
pdat2<-ddply(pdat,.(year),.fun=f)
plot(pdat2$year,pdat2$tim)
plot(pdat2$year,pdat2$amp)


modt<-lm(tim~year,data=pdat2)
moda<-lm(amp~year,data=pdat2)

p1<-predict(modt,newdata=data.frame(year=1982))
p2<-predict(modt,newdata=data.frame(year=2012))
p1<-predict(moda,newdata=data.frame(year=1982))
p2<-predict(moda,newdata=data.frame(year=2012))
