library(lattice)
library(plyr)
library(mgcv)
library(devtools)
install_github('dankelley/oce',ref='develop')
install_github('denkelley/ocedata',ref='master')
library(lubridate)
library(oce)
datadir<-'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/data/abiotic/temperature'
datadir<-'C:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/abiotic/temperature'
setwd(datadir)
hal<-read.odf('Halifax_SST.odf')
hal[['temperature']]
h<-read.odf('BoothBay_SST.odf')

lur<-read.table('Lurcher_Shoals_0m.odf',header=F,col.names=c('date','sst','nsst','sal','nsal'),sep="",skip=140)
lur$lon<--66.18
lur$lat<-43.5

bbay<-read.table('BoothBay_SST.odf',header=F,col.names=c('date','sst'),sep="",skip=100)
bbay$lon<--69.64
bbay$lat<-43.84

prin<-read.table('Prince5_Monthly_Means.odf',header=F,col.names=c('date','depth','sst','sal','nsal'),sep="",skip=130)
prin$lon<--66.81
prin$lat<-44.94
prin<-subset(prin,depth<10)

grg<-read.table('Eastern_Georges_Bank_0m.odf',header=F,col.names=c('date','sst','nsst','sal','nsal'),sep="",skip=140)
grg$lon<--66.45
grg$lat<-41.85

sta<-read.table('Standrew_SST.odf',header=F,col.names=c('date','sst'),sep="",skip=80)
sta$lon<--67.09
sta$lat<-45.08

hal<-read.table('Halifax_SST.odf',header=F,col.names=c('date','sst'),sep="",skip=90)
hal$lon<--63.57
hal$lat<-44.65



dfun<-function(d,nm){
d$date2<-strptime(substr(d$date,1,11),'%d-%b-%Y')
d$year<-as.numeric(substr(d$date2,1,4))
d$month<-as.numeric(substr(d$date2,6,7))
d$day<-yday(d$date2)
d<-subset(d,sst> -3)
d2<-subset(d,month>=8 & month<=11)
    ds1<-data.frame(year=sort(unique(d$year)),
     nmonth=tapply(d$month,d$year,function(x) length(unique(x))))
    ds2<-data.frame(year=sort(unique(d2$year)),
     nmonthf=tapply(d2$month,d2$year,function(x) length(unique(x))))
ds<-merge(ds1,ds2,by=c('year'),all=FALSE)
ds<-subset(ds,nmonthf>=2 & nmonth>=6)
d<-subset(d,year %in% ds$year)
d$date<-NULL
d$date2<-NULL

f2<-function(dd){
    mod<-gam(sst ~ s(day,bs='cs',k=6),data=dd, gamma=1)
    pdat<-data.frame(year=unique(dd$year),
                     day=250)
    p<-predict(mod,newdata=pdat,se.fit=TRUE)
    pdat$sst<-p$fit
    pdat$sst.se<-p$se.fit
#    pdat<-subset(pdat,select=c('year','sst','sst.se'))
    pdat<-subset(pdat,select=c('year','sst'))
    return(pdat)
}
dout<-ddply(d,.(year),.fun=f2)
#names(dout)[2:3]<-gsub(' ','',paste(nm,'.',names(dout)[2:3]))
names(dout)[2]<-gsub(' ','',paste(names(dout)[2],'.',nm))
return(dout)
}
hal2<-dfun(hal,'hal')
sta2<-dfun(sta,'stalt')
grg2<-dfun(grg,'grg')
prin2<-dfun(prin,'prin')
bbay2<-dfun(bbay,'bbay')
lur2<-dfun(lur,'lur')

a<-merge(hal2,sta2,by=c('year'),all=TRUE)
a<-merge(a,grg2,by=c('year'),all=TRUE)
a<-merge(a,prin2,by=c('year'),all=TRUE)
a<-merge(a,bbay2,by=c('year'),all=TRUE)
a<-merge(a,lur2,by=c('year'),all=TRUE)
a<-a[order(a$year),]
plot(a,pch=16)
cr<-round(cor(a[,2:dim(a)[2]],use='pairwise.complete.obs'),digits=2)
rowMeans(cr)

plot(a$year,a$hal.sst,type='l',ylim=c(9,18),lwd=2)
lines(a$year,a$sta.sst,type='l',col='red',lwd=2)
lines(a$year,a$grg.sst,type='l',col='darkblue',lwd=2)
lines(a$year,a$prin.sst,type='l',col='green',lwd=2)
lines(a$year,a$bbay.sst,type='l',col='magenta',lwd=2)
lines(a$year,a$lur.sst,type='l',col='dodgerblue',lwd=2)
adat<-data.frame(year=sort(unique(a$year)),
               sst=rowMeans(a[,2:7],na.rm=TRUE))
plot(adat$year,adat$sst,type='l',ylim=c(0,16))
points(adat$year,adat$sst,pch=16)



setwd('N:/cluster_2017/scratch/spera/data/finaldat_v2')
save(a,file='azmp_temp_series.RData')


#dat<-data.frame(temp=d[['data']]$temperature,
#                salinity=d[['data']]$salinity)















dfun<-function(d,nm){
d$date2<-strptime(substr(d$date,1,11),'%d-%b-%Y')
d$year<-as.numeric(substr(d$date2,1,4))
d$month<-as.numeric(substr(d$date2,6,7))
d$day<-yday(d$date2)
d<-subset(d,sst> -3,select=c('sst','year','day'))
d$sst<-(d$sst-mean(d$sst,na.rm=TRUE))/sd(d$sst,na.rm=TRUE)
d$stn<-nm
return(d)
}
hal2<-dfun(hal,'hal')
sta2<-dfun(sta,'stalt')
grg2<-dfun(grg,'grg')
prin2<-dfun(prin,'prin')
bbay2<-dfun(bbay,'bbay')
lur2<-dfun(lur,'lur')

a<-data.frame(do.call('rbind',list(hal2,sta2,grg2,prin2,bbay2,lur2)))
xyplot(sst~day | stn,data=a,pch=15)

#GET AVERAGE SST FOR EACH DAY OF SERIES
f<-function(d){
    return(d2<-data.frame(sst=mean(d$sst,na.rm=TRUE)))
}
aa<-ddply(a,.(year,day),.fun=f,.progress='text')
bb<-ddply(a,.(stn,year,day),.fun=f,.progress='text')

plot(aa$day,aa$sst,pch=15)




centerfun<-function(d){
    #DETERMINE DAY OF MAX
    mod<-gam(sst~s(day,bs='cs',k=5),data=d,gamma=1.4)
    pdat<-data.frame(xx=seq(min(d$day),max(d$day),1))
    pdat$p<-predict(mod,newdata=data.frame(day=pdat$xx))
    mday<-subset(pdat,p==max(pdat$p))$xx[1]
    #CENTER OF TIMING OF MAX
    shft<-(365/2)-mday
    d$day2<-d$day+shft
    d$day2<-ifelse(d$day2>365,d$day2-365,d$day2)
    d$day2<-ifelse(d$day2<0,d$day2+365,d$day2)
    d$day2<-ceiling(d$day2)
    return(d)
}
aa<-centerfun(aa)
bb<-ddply(bb,.(stn),.fun=centerfun,.progress='text')

plot(aa$day,aa$sst)
plot(aa$day2,aa$sst)
plot(bb$day,bb$sst)
plot(bb$day2,bb$sst)

#CTFUNCTION
fn<-function(dd){
dd$sst<-dd$sst+3
dum<-data.frame(day=sort(unique(dd$day2)),
                sump=tapply(dd$sst,dd$day2,mean))
dum$ct1<-dum$day*dum$sump
return(data.frame(ct=sum(dum$ct1)/sum(dum$sump)))
}
ctdat<-ddply(aa,.(year),.fun=fn,.progress='text')
ctdat2<-ddply(bb,.(stn,year),.fun=fn,.progress='text')
plot(ctdat$year,ctdat$ct,ylim=c(170,200))

xyplot(ct~year |stn,data=ctdat2,pch=15)

library(tidyr)
ctdat3<-ctdat2 %>% spread(stn, ct)
plot(subset(ctdat3,select=c('bbay','grg','hal','lur','prin','stalt')),pch=15)
