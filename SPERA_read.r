library(gmt)
library(data.table)
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
library(rgdal)


#SPECIFY DIRECTORIES
datadir<-'M:/data/chl_phenology/data'
figsdir<-'C:/Users/sailfish/Documents/aalldocuments/literature/postdoc_2013/chl_phenology/figures'
setwd(datadir)




#CPR DATA
setwd('M:/data/iDiv_data/CPR/CPR_GoM')
cprdat2<-read.csv('NOAACPRDATA.csv',header=TRUE,skip=6,col.names=c('cruise','sample','year','month','mday','hour','minute','lat','lon','chl','marmap.taxcode','abundance','sciname'))
cprdat2$lon<-cprdat2$lon*-1
cprdat2<-subset(cprdat2,is.na(chl)==FALSE)
cprdat2$date<-strptime(gsub(' ','',paste(cprdat2$mday,'/',cprdat2$month,'/',cprdat2$year)),"%d/%m/%Y")
cprdat2$day<-yday(cprdat2$date)
cprdat2<-subset(cprdat2,select=c('lon','lat','chl','year','month','day','date'))
cprdat2$db<-'GOM'

setwd(datadir)
cprdat<-read.csv('CPR.2015.csv',header=TRUE)
names(cprdat)<-tolower(names(cprdat))
names(cprdat)[3:7]<-c('time','ampm','lat','lon','chl')
cprdat$date<-strptime(cprdat$sampledate,"%d/%m/%Y")
cprdat$year<-as.numeric(substr(cprdat$date,1,4))
cprdat$day<-yday(cprdat$date)
cprdat$month<-month(cprdat$date)
cprdat<-subset(cprdat,select=c('lon','lat','chl','year','month','day','date'))
cprdat$db<-'SAFOS'

cprdat<-rbind(cprdat,cprdat2)
cprdat<-subset(cprdat,lon>= -66.7 & lon<= -66 & lat>= 43 & lat<= 44)

map('world',xlim=c(-70,-65),ylim=c(42,45))
points(cprdat$lon,cprdat$lat,pch=16,cex=.5,col='red')
map.axes()

dum<-data.frame(year=sort(unique(cprdat$year)),
   chl=tapply(cprdat$chl,cprdat$year,mean),
   nmonths=tapply(cprdat$month,cprdat$year,function(x) length(unique(x))))
plot(dum$year,dum$chl,type='l')

mod<-gam(chl~s(day,bs='cc',k=5) + factor(year),data=cprdat)
dayp<-182
yearp<-sort(unique(cprdat$year))
p<-predict(mod,newdata=data.frame(day=dayp,year=yearp))
pred<-data.frame(year=yearp,
                 p=p)
pred<-subset(pred,year>=1998)
plot(pred$year,pred$p,pch=16,col='red')
lines(pred$year,pred$p,col='red')
plot(cprdat$date,cprdat$chl,pch=16,cex=.5)
plot(cprdat$month,cprdat$chl)

setwd('M:/data/Spera/chl')
rs<-read.csv('chl_sw_mod_lurcher.csv',header=TRUE)
sw<-rs[,1:5];names(sw)<-c('schl','schl2','mday','month','year')
md<-rs[,6:10];names(md)<-c('mchl','mchl2','mday','month','year')
sw$date<-strptime(gsub(' ','',paste(floor(sw$mday),'/',sw$month,'/',sw$year)),"%d/%m/%Y")
sw$day<-yday(sw$date)
md$date<-strptime(gsub(' ','',paste(floor(md$mday),'/',md$month,'/',md$year)),"%d/%m/%Y")
md$day<-yday(md$date)

dum<-merge(cprdat,sw,by=c('year','month'),all=FALSE)
dum<-merge(dum,md,by=c('year','month'),all=FALSE)


cor(dum$chl,dum$schl,use='pairwise.complete.obs')
cor(dum$chl,dum$mchl,use='pairwise.complete.obs')
cor(dum$schl,dum$mchl,use='pairwise.complete.obs')

#FIGURE SHOWING RELATIONSHIPS
pdf('lurcher_modis_seawifs_cpr_chl.pdf',height=10,width=10)
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(dum$chl,dum$schl,col=alpha('black',.25),pch=16,las=1,xlab='CPR',ylab='SWFS')
t<-data.frame(chl=sort(unique(dum$chl)),
              chl=tapply(dum$schl,dum$chl,function(x) mean(x,na.rm=TRUE)))
points(t$chl,t$chl.1,pch=16,col=alpha('red3',.75),cex=2)
legend('topright',legend=c('r = 0.31'),bty='n')

plot(dum$chl,dum$mchl,col=alpha('black',.25),pch=16,las=1,xlab='CPR',ylab='MODIS')
t<-data.frame(chl=sort(unique(dum$chl)),
              chl=tapply(dum$mchl,dum$chl,function(x) mean(x,na.rm=TRUE)))
points(t$chl,t$chl.1,pch=16,col=alpha('red3',.75),cex=2)
legend('topright',legend=c('r = 0.30'),bty='n')

plot(dum$mchl,dum$schl,pch=16,xlab='MODIS',ylab='SWFS',col=alpha('darkred',.5),las=1)
abline(a=0,b=1,lty=2)
legend('topleft',legend=c('r = 0.62'),bty='n')

f<-function(d,cl){
names(d)[1]<-'y'
mod<-gam(y~s(day,bs='cc',k=5) + factor(year),data=d)
dayp<-182
yearp<-sort(unique(d$year))
p<-predict(mod,newdata=data.frame(day=dayp,year=yearp))
pred<-data.frame(year=yearp,
                 p=p)
pred$p<-(pred$p-mean(pred$p))/sd(pred$p)
lines(pred$year,pred$p,col=cl)
}
plot(0,0,xlim=c(1997,2015),ylim=c(-3,3),xlab='Year',ylab='Z-score',las=1)
abline(h=0,lty=2)
f(subset(cprdat,select=c('chl','year','day')),'green')
f(subset(sw,select=c('schl','year','day')),'red3')
f(subset(md,select=c('mchl','year','day')),'blue3')
legend('bottomright',legend=c('CPR','SWFS','MODIS'),col=c('green','red3','blue3'),lwd=2,bty='n')

dev.off()




plot(sw$year,sw$schl)
pred<-merge(pred,swf,by=c('year'),all.x=TRUE,all.y=FALSE)
plot(pred$p,pred$schl,pch=15)
cor(pred$p,pred$schl)
plot(pred$year,pred$p,type='l',col='green')
lines(pred$year,pred$schl,type='l',col='red')
