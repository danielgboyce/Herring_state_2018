
n<-number of sets
k<-outcome interested in testing
p<-probability of outcome (catchability)
d<-rbinom(n,k,prob=p)
hist(d,breaks=100)

p<-.02#probability of capture
1-dbinom(0,size=30,prob=p)#probability of detecting a herring if it is present

p<-0.02

1-dbinom(0,size=230,prob=p)#probability of detecting a herring if it is present



library(rgdal)
library(maptools)
library(segmented)
library(tidyr)
library(DataCombine)
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



setwd('C:/Users/copepod/Documents/aalldocuments/literature/research/active/ESS_trophic_control/data')
setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/ESS_trophic_control/data')
plg<-readShapePoly('polygons_ecnasap.shp')#COMMAND TO READ BACK IN
plg<-subset(plg,region=='NS')
plot(plg)
text(plg,plg$stratum)
map('world',add=TRUE,col='gray',fill=TRUE)


mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast<-readOGR('N://data/shapefiles/naturalearthdata_ne_10m_land_poly',layer='ne_10m_land')#works for Alfcoast<-fortify(coast)
coast.mc<-crop(coast,extent(-180,180,-90,90),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_blank(),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position="right",
                         plot.title = element_text(size=16)))





##############################################################

##############################################################
setwd(datadir)
rvl<-read.csv("herring_lengths_RV_survey_spera_spawnar.csv",header=TRUE)
rvl$flen<-ifelse(rvl$mission=='NED2016016',rvl$flen/10,rvl$flen)
rvl$flen<-ifelse(rvl$mission!='NED2016016',(rvl$flen*1.0866)+0.95632,rvl$flen)
#plot(log10(rvl$flen),log10(rvl$fwt))
rvl<-subset(rvl,log10(rvl$flen)>=.5 & month %in% c(6,7,8))#REMOVE OUTLIERS

rvl$lonc<-round(rvl$lon,digits=0)
rvl$lonc<-ifelse(rvl$lon<=rvl$lonc,rvl$lonc-.25,rvl$lonc+.25)
rvl$latc<-round(rvl$lat,digits=0)
rvl$latc<-ifelse(rvl$lat>=rvl$latc,rvl$latc+.25,rvl$latc-.25)
rvl$cell<-gsub(' ','',paste(rvl$lonc,'_',rvl$latc))


plot(log10(rvl$flen),log10(rvl$fwt))
a<-rvl
a$lfwt<-log10(a$fwt)
a$lflen<-log10(a$flen)
mod<-lm(lfwt ~ lflen,data=a)
abline(mod,col='red')
rvl$fwtest<-10^predict(mod,newdata=data.frame(lflen=a$lflen))
rvl$fwtest<-ifelse(is.na(rvl$fwt)==TRUE,rvl$fwtest,rvl$fwt)


rvl2<-rvl
rvl2$jcat<-ifelse(rvl$flen<=24.75,'J','A')
rvl2$id<-gsub(' ','',paste(rvl2$mission,'_',rvl2$setno))
rvl2$id1<-gsub(' ','',paste(rvl2$mission,'_',rvl2$setno,'_',rvl2$fshno))

a<-subset(rvl2,jcat=='J')
b<-subset(rvl2,jcat=='A')
sum(a$fwtest,na.rm=TRUE)/sum(rvl2$fwtest)
sum(b$fwtest,na.rm=TRUE)/sum(rvl2$fwtest)


f<-function(d){
return(data.frame(year=unique(d$year),
                  lon=unique(d$lon),
                  lat=unique(d$lat),
                  tbin5=unique(d$tbin5),
                  strat=unique(d$strat),
                  time=unique(d$time),
                  no=sum(d$clen,na.rm=TRUE),
                  wt=sum(d$fwtest/1000,na.rm=TRUE)))
}
rvll<-ddply(rvl2,.(id,jcat),.fun=f,.progress='text')


#ADDS 0'S FOR EACH LENGTH CATEGORY
dt<-unique(subset(rvll,select=c('jcat')))
f<-function(d){
d2<-merge(dt,d,by=c('jcat'),all.x=TRUE,all.y=FALSE)
#d2<-rbind.fill(dt,d)
d2$id<-unique(d$id)
d2$tbin5<-unique(d$tbin5)
d2$year<-unique(d$year)
d2$lon<-unique(d$lon)
d2$lat<-unique(d$lat)
d2$strat<-unique(d$strat)
d2$time<-unique(d$time)
d2$no<-ifelse(is.na(d2$no)==TRUE,0,d2$no)
d2$wt<-ifelse(is.na(d2$wt)==TRUE,0,d2$wt)
return(d2)
}
rvll2<-ddply(rvll,.(id),.fun=f,.progress='text')



f<-function(d){
gblon<--66.3
gblat<-43.3
if(length(unique(d$year))>=5 & length(unique(d$lat))>10 & length(unique(d$time))>10 & length(unique(d$lat))>10){

#PREDICT ANNUAL AVERAGE VALUES
modw<-gam(wt ~ as.factor(year) + s(lon,lat,k=30) + s(time,bs='cc',k=5),data=d,family='nb'(link='log'))
modn<-gam(no ~ as.factor(year) + s(lon,lat,k=30) + s(time,bs='cc',k=5),data=d,family='nb'(link='log'))
#modw<-gam(wt ~ as.factor(year)+s(time,bs='cc',k=5),data=d,family='nb'(link='log'))
#modn<-gam(no ~ as.factor(year)+s(time,bs='cc',k=5),data=d,family='nb'(link='log'))
    pdat<-data.frame(year=sort(unique(d$year)),
                     lon=gblon,
                     lat=gblat,
                     time=1200)
    pdat$pn<-predict(modn,newdata=pdat,type='response')
    pdat$pw<-predict(modw,newdata=pdat,type='response')
return(pdat)
} else NULL
}
ot<-ddply(rvll2,.(jcat),.fun=f,.progress='text')

a<-subset(ot,jcat=='A')
plot(a$year,a$pw,pch=15,las=1,xlab='Year',ylab='Predicted weight per tow [Adults]')
plot(a$year,a$pn,pch=15,las=1,xlab='Year',ylab='Predicted numbers per tow [Adults]')
b<-subset(ot,jcat=='J')
plot(b$year,b$pw,pch=15,las=1,xlab='Year',ylab='Predicted weight per tow [Juveniles]')
plot(b$year,b$pn,pch=15,las=1,xlab='Year',ylab='Predicted numbers per tow [Juveniles]')

plot(a$year,a$pw,pch=15,las=1,xlab='Year',ylab='Predicted weight per tow [Adults]')
par(new=TRUE)
plot(b$year,b$pw,pch=15,col='red',xlab='',ylab='',yaxt='n')
plot(a$pw,b$pw)
plot(log10(a$pw+1),log10(b$pw+1))
abline(a=0,b=1)
xx<-merge(a,b,by=c('year'))
ccf(as.ts(log10(xx$pn.x+1)),as.ts(log10(xx$pn.y+1)))
acf(as.ts(xx$pn.x))

xyplot(pw~year|as.factor(jcat),data=ot,pch=15)
xyplot(pn~year|as.factor(jcat),data=ot,pch=15)
a<-subset(ot,jcat=='J')
plot(a$year,a$pw,pch=15)








#############################################################

############ ESIMATE TRENDS IN HERRING BY INFERRED AGE
##################################################################
#LOOK AT LENGTH-AGE DATA TO PREDICT LENGTH OF JUVENILES V ADULTS
setwd(datadir)
adat<-read.csv("herring_assess_2016_len_wt_atage.csv",header=TRUE,na.strings=c('- ',' - '))
names(adat)<-tolower(names(adat))
names(adat)<-gsub('\\.','',names(adat))
adat<- adat %>% gather(age, y, age1:age11)
adat$agen<-as.numeric(gsub('age','',adat$age))
a<-unique(subset(adat,var=='no.x1000',select=c('age','agen','db','y')))
names(a)[4]<-c('no.x1000')
adat<-subset(adat,!(var=='no.x1000'))
adat<-merge(adat,a,by=c('db','age','agen'),all.x=TRUE,all.y=FALSE)


setwd(figsdir)
pdf('length_age_assess2016.pdf',height=7, width=8)
a<-subset(adat,var=='len.cm')
plot(a$agen,a$y,xlim=c(0,11),ylim=c(0,35),las=1,pch=15,col=alpha('black',.3),xaxt='n',xlab='Age',ylab='Length [cm]')
axis(1,seq(0,11,1))
mod <- nls(y ~ b0*(1-exp(-b1 * agen)), a, start=c(b0=1,b1=1),control=nls.control(maxiter=500),algorithm='port')
pdat<-data.frame(agen=seq(0,11,.01))
pdat$p<-predict(mod,newdata=pdat)
lines(pdat$agen,pdat$p)
pdat<-data.frame(agen=seq(0,11,.5))
pdat$p<-predict(mod,newdata=pdat)
aa<-subset(pdat,agen==3.5)
points(aa$agen,aa$p,pch=16,col='firebrick3',cex=2)
lines(c(-1,aa$agen),c(aa$p,aa$p),col='red',lty=2)
lines(c(aa$agen,aa$agen),c(-5,aa$p),col='red',lty=2)
#PREDICTED CUTOFF FOR JUVENILE/ADULT=24.8CM
agedat<-data.frame(agen=seq(0,11,.00001))
agedat$flen<-predict(mod,newdata=agedat)


a<-subset(adat,var=='len.cm')
plot(a$y,a$agen,las=1,pch=15,col=alpha('black',.3),ylab='Age',xlab='Length [cm]',ylim=c(0,11),xlim=c(0,40))
axis(1,seq(0,40,5))
md<-lm(agen~y,data=a)
mds<-summary(md)
k.strt<-mds$coef[2,1]
modage <- nls(agen ~ .1*exp(k * y), a, start=c(k=k.strt),control=nls.control(maxiter=500),algorithm='port')
pdat<-data.frame(y=seq(0,max(a$y,na.rm=TRUE),length.out=1000))
pdat$p<-predict(modage,newdata=pdat)
lines(pdat$y,pdat$p)
dev.off()



setwd(datadir)
#rvl<-read.csv("herring_lengths_RV_survey_spera_spawnar.csv",header=TRUE)
rvl<-read.csv("herring_lengths_RV_survey_spera_allar.csv",header=TRUE)
rvl$flen<-ifelse(rvl$mission=='NED2016016',rvl$flen/10,rvl$flen)
rvl$flen<-ifelse(rvl$mission!='NED2016016',(rvl$flen*1.0866)+0.95632,rvl$flen)
#plot(log10(rvl$flen),log10(rvl$fwt))
rvl<-subset(rvl,log10(rvl$flen)>=.5 & month %in% c(6,7,8)  & lat<46)#REMOVE OUTLIERS
rvl$strat<-as.character(rvl$strat)
rvl<-subset(rvl,!(strat %in% c("5Z1","5Z2","5Z9",'','558','559','440','441','442','445','446','447','443','444','559','447','449','448','450','496','451')))

rvl$lonc<-round(rvl$lon,digits=0)
rvl$lonc<-ifelse(rvl$lon<=rvl$lonc,rvl$lonc-.25,rvl$lonc+.25)
rvl$latc<-round(rvl$lat,digits=0)
rvl$latc<-ifelse(rvl$lat>=rvl$latc,rvl$latc+.25,rvl$latc-.25)
rvl$cell<-gsub(' ','',paste(rvl$lonc,'_',rvl$latc))
rvl$flen<-round(rvl$flen,digits=2)


#plot(log10(rvl$flen),log10(rvl$fwt))
a<-rvl
a$lfwt<-log10(a$fwt)
a$lflen<-log10(a$flen)
mod<-lm(lfwt ~ lflen,data=a)
abline(mod,col='red')
rvl$fwtest<-10^predict(mod,newdata=data.frame(lflen=a$lflen))
rvl$fwtest<-ifelse(is.na(rvl$fwt)==TRUE,rvl$fwtest,rvl$fwt)

agedat$flen<-round(agedat$flen,digits=2)
agedat<-unique(agedat)
agedat<-data.frame(flen=sort(unique(agedat$flen)),
                   agen=tapply(agedat$agen,agedat$flen,mean))
rvl<-merge(rvl,agedat,by=c('flen'),all.x=TRUE,all.y=FALSE)
rvl$agen<-ifelse(is.na(rvl$agen)==TRUE,10,rvl$agen)
rvl$agen<-floor(rvl$agen)+1
rvl$agen<-ifelse(rvl$agen>=8,8,rvl$agen)


rvl2<-rvl
rvl2$id<-gsub(' ','',paste(rvl2$mission,'_',rvl2$setno))
rvl2$id1<-gsub(' ','',paste(rvl2$mission,'_',rvl2$setno,'_',rvl2$fshno))

f<-function(d){
return(data.frame(year=unique(d$year),
                  lon=unique(d$lon),
                  lat=unique(d$lat),
                  tbin5=unique(d$tbin5),
                  strat=unique(d$strat),
                  time=unique(d$time),
                  no=sum(d$clen,na.rm=TRUE),
                  wt=sum(d$fwtest/1000,na.rm=TRUE)))
}
rvll<-ddply(rvl2,.(id,agen),.fun=f,.progress='text')


#ADDS 0'S FOR EACH LENGTH CATEGORY
dt<-unique(subset(rvll,select=c('agen')))
f<-function(d){
d2<-merge(dt,d,by=c('agen'),all.x=TRUE,all.y=FALSE)
#d2<-rbind.fill(dt,d)
d2$id<-unique(d$id)
d2$tbin5<-unique(d$tbin5)
d2$year<-unique(d$year)
d2$lon<-unique(d$lon)
d2$lat<-unique(d$lat)
d2$strat<-unique(d$strat)
d2$time<-unique(d$time)
d2$no<-ifelse(is.na(d2$no)==TRUE,0,d2$no)
d2$wt<-ifelse(is.na(d2$wt)==TRUE,0,d2$wt)
return(d2)
}
rvll2<-ddply(rvll,.(id),.fun=f,.progress='text')

#sdata<-rvl2[,!c('flen','fwt')]


f<-function(d){
gblon<--66.3
gblat<-43.3
if(length(unique(d$year))>=5 & length(unique(d$lat))>10 & length(unique(d$time))>10 & length(unique(d$lat))>10){

#PREDICT ANNUAL AVERAGE VALUES
modw<-gam(wt ~ as.factor(year) + s(lon,lat,k=10) + s(time,bs='cc',k=5),data=d,family='nb'(link='log'))
modn<-gam(no ~ as.factor(year) + s(lon,lat,k=10) + s(time,bs='cc',k=5),data=d,family='nb'(link='log'))
    pdat<-data.frame(year=sort(unique(d$year)),
                     lon=gblon,
                     lat=gblat,
                     time=1200)
    pdat$pn<-predict(modn,newdata=pdat,type='response')
    pdat$pw<-predict(modw,newdata=pdat,type='response')
pdat$lpn<-log10(pdat$pn+1)
pdat$lpw<-log10(pdat$pw+.01)
pdat$pnz<-(pdat$lpn-mean(pdat$lpn))/sd(pdat$lpn)
pdat$pwz<-(pdat$lpw-mean(pdat$lpw))/sd(pdat$lpw)
return(pdat)
} else NULL
}
#ot<-ddply(subset(rvll,lat<=44),.(lcat),.fun=f)#exclude smallest
#ot<-ddply(subset(rvll,lat<=44 & lon< -60),.(lcat),.fun=f)
ot<-ddply(rvll2,.(agen),.fun=f,.progress='text')

f<-function(d){
return(data.frame(totno=sum(d$no),
                  totwgt=sum(d$wt),
                  lon=mean(d$lon),
                  lat=mean(d$lat)))
}
dd<-ddply(rvll2,.(strat,agen),.fun=f)

a<-subset(dd,agen==1)
a<-subset(dd,agen==4)
plot(a$lon,a$lat,pch=16,cex=rescale(a$totno,newrange=c(.5,7)))
map('world',add=TRUE,col='gray',fill=TRUE)

d<-subset(rvll2,year==1970)

f<-function(d){
print(unique(d$year))
j<-subset(d,agen<=4)
a<-subset(d,agen>4)
if(dim(j)[1]>0){
return(data.frame(padult=mean(a$no)/mean(d$no),
                  pjuv=mean(j$no)/mean(d$no),
                  rt=mean(a$no)/mean(j$no)))
} else NULL
}
dd<-ddply(rvll2,.(year),.fun=f)
plot(dd$year,dd$padult)
plot(dd$year,dd$pjuv)
plot(dd$year,dd$rt,log='y',las=1,type='b')
plot(rvll2$lon,rvll2$lat,pch=16)
plot(plg,add=TRUE)

xyplot(pnz~year|as.factor(agen),data=ot,pch=15)
plot(ot$agen,ot$year,pch=16,cex=rescale(ot$lpn,newrange=c(.2,5)),col=alpha('darkred',.5),las=1)

pltfun<-function(ott,ttl){
names(ott)[1]<-'y'
return(ggplot()+
geom_tile(data=ott, aes(x=agen, y=year,fill=y),col='gray80',size=.0001)+
scale_fill_distiller(palette='Spectral')+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.1,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=sort(dt$agen),labels=sort(dt$agen),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(1970,2015,5),labels=as.character(seq(1970,2015,5)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(1969.5,2016.5),xlim=c(.5,8.5))+
    xlab('')+
    ylab('')+
         labs(title = ttl,
       subtitle = "",
       caption = '')
   )
}

pznum<-pltfun(subset(ot,select=c('pnz','agen','year')),'Log numbers')
pzwgt<-pltfun(subset(ot,select=c('pwz','agen','year')),'Log weight')
plnum<-pltfun(subset(ot,select=c('lpn','agen','year')),'Log numbers')
plwgt<-pltfun(subset(ot,select=c('lpw','agen','year')),'Log weight')
pnum<-pltfun(subset(ot,select=c('pn','agen','year')),'Numbers')
pwgt<-pltfun(subset(ot,select=c('pw','agen','year')),'Weight')


setwd(figsdir)
pdf('herring_rv_trends_bylength.pdf',height=8,width=10)
grid.arrange(plnum,plwgt,ncol=2)
grid.arrange(pnum,pwgt,ncol=2)
xyplot(lpn~year|as.factor(lcat),data=ot,type=c('p','spline'),pch=16,main='Log numbers')
xyplot(lpw~year|as.factor(lcat),data=ot,type=c('p','spline'),pch=16,main='Log weight')



f<-function(d){
gblon<--66.3
gblat<-43.3
if(length(unique(d$year))>=5 & length(unique(d$lat))>10 & length(unique(d$time))>10 & length(unique(d$lat))>10){

#PREDICT ANNUAL AVERAGE VALUES
mod<-gam(wt ~ year + s(lon,lat,k=10) + s(time,bs='cc',k=5),data=d,family='nb')
s<-summary(mod)
return(data.frame(b=s$p.table[2,1],
                  se=s$p.table[2,2]))
} else NULL
}
ot2<-ddply(rvll2,.(lcat),.fun=f,.progress='text')
ot2$lcat<-as.numeric(as.character(ot2$lcat))


f<-function(d){
    return(data.frame(lcat=sort(unique(d$lcat)),
              no=sum(d$no),
              wt=sum(d$wt)))
}
o<-ddply(rvll2,.(lcat),.fun=f)
ot2<-merge(ot2,o,by=c('lcat'),all=TRUE)

par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(ot2$lcat,ot2$b,las=1,pch=16,xlab='Length',ylab='Rate of change over time',col='white',ylim=c(-.15,.15),xaxt='n')
axis(1,at=seq(5,45,1))
abline(h=0,lty=2)
f<-function(d){lines(c(d$lcat,d$lcat),c(d$b+(1.96*d$se),d$b-(1.96*d$se)),col=alpha('dodgerblue3',.3),lwd=2) }
zz<-dlply(ot2,.(lcat),.fun=f)
points(ot2$lcat,ot2$b,las=1,pch=16,cex=rescale(ot2$wt,newrange=c(1,4)),col=alpha('darkblue',1))
points(ot2$lcat,ot2$b,las=1,pch=1,cex=rescale(ot2$wt,newrange=c(1,4)),col=alpha('lightgray',1),lwd=.5)
dev.off()


plot(ot2$wt,ot2$b,pch=15,las=1)
abline(h=0,lty=2)
f<-function(d){lines(c(d$wt,d$wt),c(d$b+(1.96*d$se),d$b-(1.96*d$se)))}
zz<-dlply(ot2,.(lcat),.fun=f)

o<-ddply(rvll2,.(lcat,tbin5),.fun=f)
o$lcat<-as.numeric(as.character(o$lcat))
xyplot(wt~lcat |as.factor(tbin5),data=o,type=c('p','spline'),pch=15)
xyplot(no~lcat |as.factor(tbin5),data=o,type=c('p','spline'),pch=15)

plot(o$lcat,o$no,pch=15)
plot(o$lcat,o$wt,pch=15)
plot(o$lcat,log10(o$wt),pch=15)

plot(rvl2$flen,rvl2$fwt)
mod<-gam(fwt ~ s(flen,k=4),data=rvl2)
pdat<-data.frame(flen=seq(min(rvl2$flen),max(rvl2$flen),length.out=100))
pdat$p<-predict(mod,newdata=pdat)
lines(pdat$flen,pdat$p,col='red')
rvll$lwgt<-predict(mod,newdata=data.frame(flen=as.numeric(as.character(rvll$lcat))))

d<-subset(rvll,lcat==7.5)








#CHECKS TO SEE IF GENERIC EARLY WARNINGS INDICATORS ARE RELEVANT
ewfun<-function{
    library(earlywarnings)
mod<-gam(totwgt ~ as.factor(year) + s(lon,lat,k=50),data=rvw)
    pdat<-data.frame(year=sort(unique(rvw$year)),
                     lon=gblon,
                     lat=gblat)
    p<-predict(mod,newdata=pdat,se.fit=TRUE,type='response')
    pdat$p<-10^p$fit
    pdat$se<-10^p$se.fit
ew<-generic_ews(subset(pdat,select=c('year','p')),winsize=10,detrending='gaussian',interpolate=TRUE)
}




##############################################################

##ESTIMATE HERRING TRENDS FROM RV DATA AT STRATUM LEVEL
setwd(datadir)
#rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw<-read.csv("herring_weights_RV_survey_spera_allar.csv",header=TRUE)
rvw<-subset(rvw,month %in% c(6,7,8) & is.na(dmin)==FALSE)


rvw$bank<-ifelse(rvw$strat %in% c(447,448),'banq','no')
rvw$bank<-ifelse(rvw$strat %in% c(443),'mis',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(458),'mid',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(455,456),'sab',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(464),'west',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(463),'em',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(473),'lh',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(474),'rw',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(475),'bac',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(480),'bn',rvw$bank)

rvw$no<-log10(rvw$totno+1)
rvw$wgt<-log10(rvw$totwgt+1)
rvw$sz<-rvw$totwgt/rvw$totno
rvw$sz<-ifelse(rvw$sz==Inf,rvw$totwgt,rvw$sz)

rvw$id<-gsub(' ','',paste(rvw$mission,'_',rvw$setno))
rvw$pres<-ifelse(rvw$totno>0,1,0)
rvw$tbin20<-ifelse(rvw$year<=1990,1980,2010)

rvw$lonc<-round(rvw$lon,digits=0)
rvw$lonc<-ifelse(rvw$lon<=rvw$lonc,rvw$lonc-.25,rvw$lonc+.25)
rvw$latc<-round(rvw$lat,digits=0)
rvw$latc<-ifelse(rvw$lat>=rvw$latc,rvw$latc+.25,rvw$latc-.25)
rvw$cell<-gsub(' ','',paste(rvw$lonc,'_',rvw$latc))

rvw$cell.1<-gsub(' ','',paste(round(rvw$lon,digits=1),'_',round(rvw$lat,digits=1)))
lonc.1<-seq(min(rvw$lon),max(rvw$lon),.1)
latc.1<-seq(min(rvw$lat),max(rvw$lat),.1)
crds<-expand.grid(lonc.1=lonc.1,latc.1=latc.1)
crds$cell.1<-gsub(' ','',paste(round(crds$lonc.1,digits=1),'_',round(crds$latc.1,digits=1)))
rvw<-merge(rvw,crds,by=c('cell.1'),all.x=TRUE,all.y=FALSE)


f<-function(d){
    return(data.frame(mnno=mean(d$totno),
                      mdno=median(d$totno),
                      mnwt=mean(d$totwgt),
                      mdwt=median(d$totwgt),
                      lon=median(d$lonc),
                      lat=median(d$latc)))
}
dd<-ddply(rvw,.(strat),.fun=f)
dd<-ddply(subset(rvw,month==7 & year>=2000),.(cell),.fun=f)
dd<-ddply(subset(rvw,month==8 & year>=2000),.(cell),.fun=f)
plot(dd$lon,dd$lat,pch=16,cex=rescale(log10(dd$mdwt+1),newrange=c(.5,7)),col='dodgerblue')
map('world',add=TRUE,col='gray',fill=TRUE)
plot(plg,add=TRUE,fill=FALSE)

plot(plg,col=alpha('darkblue',.4))
points(rvw$lon,rvw$lat,pch=16,col=alpha('red3',.3),cex=.5)
map('world',add=TRUE,col='gray',fill=TRUE)

dm<-data.frame(year=sort(unique(rvw$year)),
    n=tapply(rvw$totno,rvw$year,sum),
    ntow=tapply(rvw$id,rvw$year,length),
    ncell=tapply(rvw$cell.1,rvw$year,function(x) length(unique(x))),
    ntime=tapply(rvw$time,rvw$year,function(x) length(unique(x))),
    ndepth=tapply(round(rvw$dmax,digits=0),rvw$year,function(x) length(unique(x))),
    nday=tapply(rvw$day,rvw$year,function(x) length(unique(x))),
    dayrng=tapply(rvw$day,rvw$year,function(x) max(x)-min(x)))
a<-subset(dm,year<=1990)
b<-subset(dm,year>1990)
mean(a$ncell)
mean(b$ncell)
mean(a$dayrng)
mean(b$dayrng)

par(mfrow=c(2,2))
dm$n<-log10(dm$n)
plot(dm$year,dm$ncell)
plot(dm$year,dm$ntow)
plot(dm$year,dm$ntow)
plot(dm,pch=15)

setwd(figsdir)
pdf('rv_sample_effort_overtime.pdf',height=10,width=8)
par(mfrow=c(3,2),mar=c(4,4,1,1))
f<-function(d,lg){
 nm<-names(d)[2]
 names(d)<-c('year','y')
 if(lg==TRUE){d$y<-log10(d$y)
          } else NULL
 plot(d$year,d$y,las=1,xlab='Year',ylab=nm,pch=16,col=alpha('gold3',1),cex=2,type='l',lwd=2)
points(d$year,d$y,pch=16,col=alpha('gold3',.3),cex=2)
points(d$year,d$y,pch=1,col='lightgray',lwd=.01,cex=2)
}
f(subset(dm,select=c('year','n')),TRUE)
f(subset(dm,select=c('year','n')),F)
f(subset(dm,select=c('year','ntow')),F)
f(subset(dm,select=c('year','ncell')),F)
f(subset(dm,select=c('year','nday')),F)
f(subset(dm,select=c('year','dayrng')),F)
dev.off()


d<-subset(rvw,bank=='sab')

fsm<-function(d){
d$y<-d$totwgt+.1
mod<-gam(y~as.factor(year) + s(time,bs='cc',k=5) + s(lon,lat,k=4),data=d,gamma=.5,family=Gamma('log'))
pdat<-data.frame(year=sort(unique(d$year)),
                 time=1200,
                 lon=median(d$lon),
                 lat=median(d$lat))
pdat$pwt<-predict(mod,newdata=pdat,type='response')

d$y<-ifelse(d$totwgt>0,1,0)
mod2<-gam(y~as.factor(year) + s(time,bs='cc',k=5) + s(lon,lat,k=4),data=d,gamma=.5,family=binomial)
pdat$pps<-predict(mod2,newdata=pdat,type='response')
pdat<-subset(pdat,select=c('year','pwt','pps'))
pdat$bank<-unique(d$bank)
pdat$lon<-mean(d$lon)
#pdat$p<-(pdat$p-mean(pdat$p))/sd(pdat$p)
#names(pdat)[2]<-unique(as.character(d$cell))
return(pdat)
}
qq<-ddply(rvw,.(bank),.fun=fsm)
xyplot(log(pwt)~year | bank,data=qq,type=c('p','l'),pch=15)
xyplot(pps~year | bank,data=qq,pch=15)

f<-function(d){
    return(data.frame(no=mean(d$totno),
                      ntow=length(unique(d$id))))
}
ot<-ddply(rvw,.(bank),.fun=f)
ot<-ot[order(ot$no,decreasing=TRUE),]


ott<-subset(qq,bank!='no',select=c('pwt','bank','year','lon'))
ttl<-'Presence'
lg<-TRUE

pltfun<-function(ott,ttl,lg){
names(ott)[1]<-'y'
names(ott)[3]<-'year'
if(lg==TRUE){
    lms<-c(max(ott$y)*-1,0)
} else  {   lms<-c(-1,0)}
ott$bank<-as.factor(ott$bank)
dm<-unique(subset(ott,select=c('lon','bank')))
dm<-dm[order(dm$lon),]
dm$id<-seq(1,dim(dm)[1],1)
ott<-merge(ott,dm,by=c('bank'),all.x=TRUE,all.y=FALSE)
ott$y<-ott$y*-1
return(ggplot()+
geom_tile(data=ott, aes(x=id, y=year,fill=y),col='gray',size=.0001) +
scale_fill_distiller(palette='YlOrRd',limits=lms) +
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.1,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(1,10,1),labels=dm$bank,limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(1970,2015,5),labels=as.character(seq(1970,2015,5)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(1969.5,2016.5),xlim=c(.5,10.5))+
    xlab('')+
    ylab('')+
         labs(title = ttl,
       subtitle = "",
       caption = '')
       )
}
p1<-pltfun(subset(qq,bank!='no',select=c('pps','bank','year','lon')),'Presence',FALSE)
p2<-pltfun(subset(qq,bank!='no',select=c('pwt','bank','year','lon')),'Presence',TRUE)
p2
p1





f<-function(d){return(data.frame(nyear=length(unique(d$year)),
                                     myear=min(d$year)))
               }
dm<-ddply(rvw,.(strat),.fun=f)
dmm<-dm[order(dm$nyear),]
dmm<-subset(dmm,nyear>=15 & myear<1990)

fsm<-function(d){
d$y<-d$totwgt+.1
mod<-gam(y~s(year) + s(time,bs='cc',k=5) + s(lon,lat,k=4),data=d,gamma=.5,family=Gamma('log'))
pdat<-data.frame(year=seq(min(d$year),max(d$year),1),
                 time=1200,
                 lon=median(d$lon),
                 lat=median(d$lat))
pdat$pwt<-predict(mod,newdata=pdat,type='response')

d$y<-ifelse(d$totwgt>0,1,0)
mod2<-gam(y~as.factor(year) + s(time,bs='cc',k=5) + s(lon,lat,k=4),data=d,gamma=.5,family=binomial)
pdat$pps<-predict(mod2,newdata=pdat,type='response')
pdat<-subset(pdat,select=c('year','pwt','pps'))
pdat$bank<-unique(d$bank)
pdat$lon<-mean(d$lon)
pdat$lat<-mean(d$lat)
return(pdat)
}
qt<-ddply(subset(rvw,strat %in% dmm$strat),.(strat),.fun=fsm)

#devtools::install_github('dgrtwo/gganimate',force=TRUE)
install.packages('magick')
library(magick)
library(gganimate)
qt2<-subset(qt,select=c('strat','year','pwt','pps','bank'))
qt2$id<-qt2$strat

a<-plg
a$id<-a$stratum
am<-fortify(a,region='id')
am<-subset(am,id%in% qt2$id)
mydat<-merge(am,qt2,by=c('id'))
mydat$lpwt<-log10(mydat$pwt+.1)

p<-ggplot(mydat,aes(x=long,y=lat,group = group, fill=lpwt,frame=year)) +
theme_opts +
  coord_equal() +
  geom_polygon(color = 'grey',size=.0001) +
#  geom_polygon(aes(long,lat,group=group),fill=NA,colour='black',data=a) +
  labs(title = "Herring presence between 1970 and 2015",
   subtitle = "Dynamic Map",
   caption = 'Data source: Statistics Canada') +
  scale_fill_distiller(palette='Spectral')

setwd(figsdir)
  gganimate(p)
  gganimate(p,'output.gif')
  gganimate(p,'output.log.gif')
  gganimate(p,'output.mp4')
  gganimate(p,'output.smooth.mp4')


###################################################################

#ESTIMATES SMOOTH TREND IN DIFFERENT QUANTITIES OVER TIME FOR CLUSTER
f<-function(d){    return(data.frame(nyear=length(unique(d$year)),
                                     myear=min(d$year)))
               }
dm<-ddply(rvw,.(cell),.fun=f)
dmm<-dm[order(dm$nyear),]
dmm<-subset(dmm,nyear>=15 & myear<1990)

dm2<-ddply(subset(rvw,is.na(sz)==FALSE),.(cell),.fun=f)
dmm2<-subset(dm2,nyear>=15 & myear<1990)

fsm<-function(d){
nm<-names(d)[1]
names(d)[1]<-'y'
if(nm=='sz'){
d<-subset(d,is.na(y)==FALSE)
mod<-gam(y~s(year) + s(lon,lat,k=4),data=d,gamma=.5,gamily='nb')
} else {
mod<-gam(y~s(year) + s(time,bs='cc',k=5) + s(lon,lat,k=4),data=d,gamma=.5,gamily='nb')
}
pdat<-data.frame(year=seq(min(d$year),max(d$year),.25),
                 time=1200,
                 lon=median(d$lon),
                 lat=median(d$lat))
pdat$p<-predict(mod,newdata=pdat,type='response')
pdat<-subset(pdat,select=c('year','p'))
pdat$p<-(pdat$p-mean(pdat$p))/sd(pdat$p)
names(pdat)[2]<-unique(as.character(d$cell))
return(pdat)
}
#TOTAL WEIGHT
qq<-dlply(subset(rvw,cell %in% dmm$cell,select=c('totwgt','strat','year','cell','time','lon','lat')),.(cell),.fun=fsm)
qsm<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), qq)#COMBINE
qsm<-qsm[,colSums(is.na(qsm)) != nrow(qsm)]#REMOVES COLUMNS THAT ARE ALL MISSING

#TOTAL NUMBERS
qq<-dlply(subset(rvw,cell %in% dmm$cell,select=c('totno','cell','year','time','lon','lat')),.(cell),.fun=fsm)
qsm2<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), qq)#COMBINE
qsm2<-qsm2[,colSums(is.na(qsm2)) != nrow(qsm2)]#REMOVES COLUMNS THAT ARE ALL MISSING

#AVERAGE SIZE
qq<-dlply(subset(rvw,sz!=Inf & is.na(sz)==FALSE & cell %in% dmm2$cell,select=c('sz','cell','year','time','lon','lat')),.(cell),.fun=fsm)
qsm3<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), qq)#COMBINE
qsm3<-qsm3[,colSums(is.na(qsm3)) != nrow(qsm3)]#REMOVES COLUMNS THAT ARE ALL MISSING


q<-qsm
q2<- q %>% gather(cell, value, -year)
xyplot(value ~ year | strata,data=q2, pch=15, type=c('spline'),col='black')
xyplot(value ~ year | strata,data=q2, pch=15, type=c('p','spline'),col='black')




clfun<-function(df,k,lbl){
rownames(df)<-df$year
df<-df[,-1]
#k<-3#ER OF CLUSTERS
dmat<-1-cor(df,use='pairwise.complete.obs')
dst<-as.dist(dmat)
ff<-fanny(dst,k,maxit=5000,diss=T)

par(mfrow=c(3,1),mar=c(4,12,1,12))
dum<-c('red3','darkblue','gold3')
plot(silhouette(ff),col=dum[1:k],main='')#silhouette plot

dc.pcoa<-cmdscale(dst)
dc.scores<-scores(dc.pcoa,choices=c(1,2))
spefuz.g<-ff$clustering
a<-data.frame(cell=as.character(sort(unique(names(df)))),
              clusters=ff$clustering)
aa<-data.frame(ff$membership)
aa$cell<-rownames(a)

plot(scores(dc.pcoa),asp=1,type='n',xlim=c(-1,1.5),ylim=c(-1,.5),las=1,axes=TRUE,xlab='',ylab='')
stars(ff$membership,location=scores(dc.pcoa),draw.segments=T,add=T,scale=F,len=.1,col.segments=alpha(c(dum[1:k]),.25),byt='n',labels=NULL,xlim=c(-1.1,1.4),ylim=c(-1,.5),lwd=.0001,xpd=TRUE,border=NULL,radius=FALSE,col.radius=alpha('white',.1))
for(i in 1:k){ cl<-dum[i]
    gg<-dc.scores[spefuz.g==i,]
    hpts<-chull(gg)
    hpts<-c(hpts,hpts[1])
    lines(gg[hpts,],col=cl,lwd=3,xlim=c(-1.1,1.4),ylim=c(-1,.5))
}

cx <- data.frame(cell=aa$cell,
                 cx=apply(aa[, 1:k], 1, max))
cx$cxx <- rescale(cx$cx,newrange=c(.05,.5))
#PLOTS THE TIMESERIES OF R FOR EACH CLUSTER
par(mfrow=c(3,1),mar=c(2,12,1,12),oma=c(1,1,1,1))

par(mar=c(4,8,4,8))
mb<-seq(1,k,1)
l<-list()
for(i in 1:length(mb)){
print(mb[i])
one<-subset(a,clusters==mb[i],select=c('cell'))# IN THIS CLUSTER
cx2<-subset(cx,cx>=0.5)#GETS INSTANCES WHERE CLUSTER PROBABILITY>0.5
data5<-subset(df,select=c(as.character(one$cell)))#TS FOR CLUSTER
cx2<-subset(cx2,cell %in% names(data5))
data5<-subset(df,select=c(as.character(cx2$cell)))#TS FOR CLUSTER
x<-rownames(data5)
t<-data.frame(year=as.numeric(x),mn=rowMeans(data5,na.rm=F))
t2<-data.frame(year=as.numeric(x),mn=rowMeans(data5,na.rm=T))

cl<-dum[i]
plot(0,0,pch=16,cex=.01,xlim=c(1970,2020),ylim=c(-2,3),main='',xaxt='n',las=1,axes=FALSE,xlab='',ylab='')
axis(side=2,at=seq(-2,3,1),las=1,lwd=.001,cex.axis=.75)
axis(side=1,at=seq(1970,2020,10),cex.axis=.75)
    for(j in 1:length(data5[1,])){
        try(dat<-data5[,j])
        try(trnsp<-subset(cx,cell==as.character(names(data5[j])))$cxx)
        try(lines(x,dat,cex=.5,ylim=c(0,1),lwd=2,col=alpha(cl,rescale(trnsp,newrange=c(.1,.5)))))
    }
lines(as.numeric(as.character(t$year)),t$mn,col='black',cex=.7,lwd=4)
dm<-subset(t2,mn==max(t2$mn,na.rm=TRUE),select=c('year'))
dm$cluster=mb[i]
names(dm)<-c('clusterday','clusters')
l[[i]]<-dm
}


a2<-a
if(k==3){a2$cx<- apply(aa[,c('X1','X2','X3')], 1, function(x) max(x) )
     } else {a2$cx<- apply(aa[,c('X1','X2')], 1, function(x) max(x) )
     }
a2$id<-a2$cell
a2$cl<-ifelse(a2$clusters==1,dum[1],dum[3])
a2$cl<-ifelse(a2$clusters==2,dum[2],a2$cl)
crds<-unique(subset(rvw,select=c('cell','lonc','latc')))
a2<-merge(a2,crds,by=c('cell'),all.x=TRUE,all.y=FALSE)
dum<-unique(subset(a2,select=c('clusters','cl')))

return(ggplot()+
geom_tile(data=a2, aes(x=lonc, y=latc,fill=as.factor(clusters),alpha=cx),col='gray80',size=.0001)+
scale_fill_manual(breaks=as.character(dum$clusters),values=dum$cl,na.value="transparent",guide=guide_legend(title=''))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.9,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.1, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-56,1),labels=as.character(seq(-68,-56,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,48,1),labels=as.character(seq(41,48,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,48),xlim=c(-68,-56))+
    xlab('')+
    ylab('')+
         labs(title = lbl,
       subtitle = "",
       caption = '')
   )
}


setwd(figsdir)
pdf('herring_rvsurvey_cluster_all_grid_k3.pdf',height=10,width=6.5)
mp1<-clfun(qsm,3,"Total observed herring (weight)")
mp2<-clfun(qsm2,3,"Total observed herring (numb)")
mp3<-clfun(qsm3,3,"Average weight")
dev.off()

pdf('herring_rvsurvey_cluster_map_grid_k3.pdf',height=14,width=7)
grid.arrange(mp1,mp2,mp3,ncol=1)
dev.off()






#MAP OF TRENDS IN ABUNDANCE AND AVERAGE WEIGHT
#ESTIMATES SMOOTH TREND IN DIFFERENT QUANTITIES OVER TIME FOR CLUSTER
f<-function(d){    return(data.frame(nyear=length(unique(d$year)),
                                     myear=min(d$year),
                                     mx=max(d$totwgt)))
               }
dm<-ddply(rvw,.(cell),.fun=f)
dmm<-dm[order(dm$nyear),]
dmm<-subset(dmm,nyear>=15 & myear<1990 & mx>0)


fsm<-function(d){
names(d)[1]<-'y'
#d$y<-(d$y-mean(d$y))/sd(d$y)
mod<-gam(y~s(year) + s(time,bs='cc',k=5) + s(lon,lat,k=4),data=d,gamma=.5,gamily='nb')
mod2<-gam(y~year + s(time,bs='cc',k=5) + s(lon,lat,k=4),data=d,gamma=.5,gamily='nb')
s<-summary(mod2)
pdat<-data.frame(year=seq(min(d$year),max(d$year),.25),
                 time=1200,
                 lon=median(d$lon),
                 lat=median(d$lat))
pdat$p<-predict(mod,newdata=pdat,type='response')
ot<-data.frame(pstart=pdat$p[1],
               pend=pdat$p[dim(pdat)[1]],
               lonc=unique(d$lonc),
               latc=unique(d$latc),
               span=max(pdat$year)-min(pdat$year),
               year1=min(pdat$year),
               beta=s$p.table[2,1],
               pv=s$p.table[2,4])
ot$chng<-ot$pend-ot$pstart
d$y<-(d$y-mean(d$y))/sd(d$y)
mod<-gam(y~s(year,k=4) + s(lon,lat,k=4),data=d,gamma=.5,gamily='nb')
p<-predict(mod,newdata=pdat,type='response',se.fit=FALSE)
ot$chngz<-p[length(p)]-p[1]
return(ot)
}
#TOTAL WEIGHT
mdat<-ddply(subset(rvw,cell %in% dmm$cell,select=c('totwgt','strat','year','cell','time','lon','lat','lonc','latc')),.(cell),.fun=fsm)



###SAME BUT FOR AVERAGE WEIGHT
rvws<-subset(rvw,is.na(sz)==FALSE)
f<-function(d){    return(data.frame(nyear=length(unique(d$year)),
                                     myear=min(d$year),
                                     mx=max(d$totwgt)))
               }
dm<-ddply(rvws,.(cell),.fun=f)
dmm<-dm[order(dm$nyear),]
dmm<-subset(dmm,nyear>=15 & myear<1990 & mx>0)

d<-subset(rvws,cell=="-58.25_45.75",select=c('sz','strat','year','cell','time','lon','lat','lonc','latc'))

fsm<-function(d){
print(unique(d$cell))
names(d)[1]<-'y'
mod<-gam(y~s(year,k=4) + s(lon,lat,k=4),data=d,gamma=.5,gamily='nb')
mod2<-gam(y~year + s(lon,lat,k=4),data=d,gamma=.5,gamily='nb')
s<-summary(mod2)
pdat<-data.frame(year=seq(min(d$year),max(d$year),.25),
                 lon=median(d$lon),
                 lat=median(d$lat))
pdat$p<-predict(mod,newdata=pdat,type='response')
ot<-data.frame(pstart=pdat$p[1],
               pend=pdat$p[dim(pdat)[1]],
               lonc=unique(d$lonc),
               latc=unique(d$latc),
               span=max(pdat$year)-min(pdat$year),
               year1=min(pdat$year),
               beta=s$p.table[2,1],
               pv=s$p.table[2,4])
ot$chng<-ot$pend-ot$pstart

d$y<-(d$y-mean(d$y))/sd(d$y)
mod<-gam(y~s(year,k=4) + s(lon,lat,k=4),data=d,gamma=.5,gamily='nb')
p<-predict(mod,newdata=pdat,type='response',se.fit=FALSE)
ot$chngz<-p[length(p)]-p[1]
return(ot)
}
#AVERAGE WEIGHT
mdat2<-ddply(subset(rvws,cell %in% dmm$cell,select=c('sz','strat','year','cell','time','lon','lat','lonc','latc')),.(cell),.fun=fsm)


pfun<-function(a,mx,dg,lbl){
names(a)[1]<-'y'
a$y<-ifelse(a$y>mx,mx,a$y)
a$y<-ifelse(a$y< -mx,-mx,a$y)
aa<-data.frame(y=seq((-mx)-.001,max(abs(a$y)+.001,na.rm=TRUE),length.out=100))
a<-rbind.fill(a,aa)

n<-21
mxx<-max(abs(a$y))
brks<-seq((-mx)-0.001,mxx+.001,length.out=n)
brks2<-round(seq((-mx)-0.001,mxx+.001,length.out=n),digits=dg)
a$ycat<-cut(a$y,breaks=brks)
lbls<-sort(unique(a$ycat))
lbls2<-sort(unique(cut(a$y,breaks=brks2)))
cls<-matlab.like(length(lbls))

ggplot()+
geom_tile(data=a, aes(x=lonc, y=latc,fill=ycat),col='gray80',size=.0001)+
scale_fill_manual(breaks=as.character(lbls),values=cls,labels=lbls2,na.value="transparent")+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.9,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.1, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-56,1),labels=as.character(seq(-68,-56,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,48,1),labels=as.character(seq(41,48,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,48),xlim=c(-68,-56))+
    xlab('')+
    ylab('')+
         labs(title = lbl,
       subtitle = "",
       caption = '')
}
pfun(subset(mdat2,select=c('chngz','lonc','latc')),4,2,'Average size (Z)')
pfun(subset(mdat2,select=c('chng','lonc','latc')),.4,2,'Average size')
pfun(subset(mdat,select=c('chngz','lonc','latc')),1.5,2,'Biomass trend (Z)')
pfun(subset(mdat,select=c('chng','lonc','latc')),85,2,'Biomass trend')





#######################################################

#GETS CELLS WHERE AT LEAST ONE HERRING WAS CAPTURED DURING ALL YEARS OF SURVEY
zz<-data.frame(cell.1=sort(unique(rvw$cell.1)),
               n=tapply(rvw$totno,rvw$cell.1,sum))
zzz<-subset(zz,n>0)

#GETS SUM OF ALL FISH CAPTURED BY CELL AND YEAR
f<-function(d){
    return(data.frame(n=sum(unique(d$totno)),
                      lon=unique(d$lonc.1),
                      lat=unique(d$latc.1),
                      ntows=length(unique(d$id))))
}
dt<-ddply(rvw,.(cell.1,year),.fun=f,.progress='text')

f<-function(d){
pres<-subset(d,n>0)
return(data.frame(rng=(dim(pres)[1]/dim(d)[1])*100,
                  ncells=length(unique(d$cell.1)),
                  ntows=sum(d$ntows)))
}
oo<-ddply(subset(dt,cell.1 %in% zzz$cell.1),.(year),.fun=f)
plot(oo$year,oo$rng,pch=15)
plot(oo$year,oo$ncells,pch=15)
plot(oo$year,oo$ntows,pch=15)
plot(oo$year,oo$rng/oo$ncells,pch=15)

##########################################################

################       DIURNAL CHANGES

##########################################################
pltfun<-function(ott,ttl){
names(ott)[1]<-'y'
names(ott)[3]<-'year'
ott<-na.omit(ott)
return(ggplot()+
geom_tile(data=ott, aes(x=time, y=year,fill=y),size=.0001) +
scale_fill_distiller(palette='Spectral')       +
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.1,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(0,2400,300),labels=seq(0,2400,300),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(1970,2015,5),labels=as.character(seq(1970,2015,5)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(1969.5,2016.5),xlim=c(0,2350))+
    xlab('')+
    ylab('')+
         labs(title = ttl,
       subtitle = "",
       caption = '')
   )
}

setwd(figsdir)
pdf('herring_rv_diurnal.pdf',height=8,width=10)
par(mfrow=c(2,2),mar=c(4,4,1,1))
f<-function(d){
    if((max(d$time)-min(d$time))>15){
    mod<-gam(totwgt~s(time,bs='cc',k=6)+s(lon,lat,k=20),data=d,gamma=1,family=nb)
s<-summary(mod)
    pdat<-data.frame(time=seq(min(d$time),max(d$time),1),
                     lon=median(d$lon),
                     lat=median(d$lat))
    pdat$p<-predict(mod,newdata=pdat,type='response')
    pdat$pz<-(pdat$p-mean(pdat$p,na.rm=TRUE))/sd(pdat$p,na.rm=TRUE)
    pdat$rng<-max(pdat$pz,na.rm=TRUE)-min(pdat$pz,na.rm=TRUE)
    pdat$pmx<-subset(pdat,p==max(pdat$p)[1])$time[1]
    pdat$pv<-s$s.table[1,4]
    pdat$b<-s$s.table[1,1]
    return(pdat)
} else NULL
}
ot<-ddply(rvw,.(year),.fun=f,.progress='text')
ot2<-acast(ot,year~time,value.var="pz")
image(x=sort(unique(ot$year)),y=sort(unique(ot$time)),ot2,col=palette(rich.colors(500)),xlab='Exploitation rate',ylab='Consumer trophic level',las=1,cex.axis=.8)

plot(ot$year,ot$rng,pch=15,las=1,xlab='Year',ylab='Daily range of catch')
plot(ot$year,ot$pmx,pch=15,las=1,xlab='Year',ylab='Time of max catch')
p1<-pltfun(subset(ot,select=c('pz','time','year')),'Z-score')


#SAME BUT BINNED TO 3 EYAR INTERVALS
f<-function(d){
    mod<-gam(totwgt~s(time,bs='cc',k=6)+s(lon,lat,k=30) + as.factor(year),data=d,gamma=1,family=nb)
    pdat<-data.frame(time=seq(min(d$time),max(d$time),1),
                     lon=median(d$lon),
                     lat=median(d$lat),
                     year=median(d$year))
    pdat$p<-predict(mod,newdata=pdat,type='response')
    pdat$pz<-(pdat$p-mean(pdat$p))/sd(pdat$p)
    pdat$rng<-max(pdat$pz)-min(pdat$pz)
    pdat$pmx<-subset(pdat,p==max(pdat$p)[1])$time
    return(pdat)
}
ot<-ddply(rvw,.(tbin3),.fun=f,.progress='text')
ot2<-acast(ot,tbin3~time,value.var="pz")
image(x=sort(unique(ot$tbin3)),y=sort(unique(ot$time)),ot2,col=palette(rich.colors(500)),xlab='Exploitation rate',ylab='Consumer trophic level',las=1,cex.axis=.8)

plot(ot$tbin3,ot$rng,pch=15,las=1,xlab='Year',ylab='Daily range of catch')
plot(ot$tbin3,ot$pmx,pch=15,las=1,xlab='Year',ylab='Time of max catch')
p2<-pltfun(subset(ot,select=c('pz','time','tbin3')),'Z-score')

cls<-colorRampPalette(c('black','darkmagenta','hotpink','darkblue','lightskyblue','forestgreen','lawngreen','gold','orange','firebrick1','firebrick4'))
n<-length(unique(ot$tbin3))
dum3<-data.frame(tbin3=sort(unique(ot$tbin3)),
                 cls=cls(n+2)[3:(n+2)])
ot<-merge(ot,dum3,by=c('tbin3'),all.x=TRUE,all.y=FALSE)

f<-function(d){
    lines(d$time,d$pz,col=alpha(as.character(unique(d$cl)),.5),lwd=3)
}
plot(0,0,xlim=c(0,2400),ylim=c(-1.5,2.5),las=1,xlab='Time',ylab='Herring',col='white')
zz<-dlply(ot,.(tbin3),.fun=f)
xyplot(pz~time | tbin3,data=ot)

grid.arrange(p1,p2,ncol=2)
dev.off()



dt$lcat<-as.numeric(as.character(dt$lcat))





sfun<-function(d){
    return(d[sample(nrow(d),100,replace=FALSE),])
}
rvw2<-ddply(rvw,.(year),.fun=sfun,.progress='text')

a<-subset(rvw,id==sort(unique(rvw$id))[2])

par(mfrow=c(3,4),mar=c(1,1,1,1))
yrs<-seq(1975,2015,5)
for(i in 1:length(yrs)){
d<-subset(dt,year==yrs[i] & n>0)
plot(d$lon,d$lat,pch=16,col='red',las=1)
    map('world',add=TRUE,fill=TRUE,col='lightgray')
}

oo$rng<-log10(oo$rng)
plot(oo,pch=16)

a<-subset(rvw,year<1990)
b<-subset(rvw,year>=1990)
plot(a$lonc.1,a$latc.1,pch=15)
plot(b$lonc.1,b$latc.1,pch=15)


par(mfrow=c(3,4),mar=c(1,1,1,1))
f<-function(d){
mod<-gam(totwgt~s(time,bs='cc',k=6) + s(lon,lat,k=10),data=d,family='nb',gamma=1)
 pdat<-data.frame(time=seq(min(d$time),max(d$time),length.out=100),
                  lon=-55,
                  lat=42.5)
 p<-predict(mod,newdata=pdat,type='response')
 pdat$p<-p
# plot(pdat$time,pdat$p,type='l',main=unique(d$tbin10))
return(data.frame(time=subset(pdat,p==max(pdat$p))$time[1]))
}
zz<-ddply(rvw,.(year),.fun=f,.progress='text')
plot(zz$year,zz$time,pch=15,type='b')







plot(rvw$lon,rvw$lat,pch='.')
points(rvw$lonc,rvw$latc,pch=15,col='purple')

f<-function(d){ return(data.frame(n=length(unique(d$id)),
                                  lon=unique(d$lonc.1),
                                  lat=unique(d$latc.1),
                                  totno=sum(d$totno)))}
dtt<-ddply(rvw,.(cell.1,year),.fun=f,.progress='text')
dtt$p<-1-dbinom(0,dtt$n,.02)

1-dbinom(0,1,.02)
1-dbinom(0,30,.02)
1-dbinom(0,1,.02)
1-dbinom(0,30,.02)

d<-subset(rvw,cell=="-57.25_44.25")

f<-function(d){
    if(length(unique(d$year))>=5 & length(unique(d$time))>5){
    mod<-gam(pres~as.factor(year),data=d,gamma=1.4,family='binomial')
s<-summary(mod)
    pdat<-data.frame(year=sort(unique(d$year)),
           time=1200,
           ntows=tapply(d$id,d$year,function(x) length(unique(x))))
    pdat$pabs<-1-dbinom(0,pdat$ntows,.02)#probability that 0 is real
length(c(1,s$p.table[,4]))
    pdat$pprs<-s$p.table[,4]
    pdat$pprs[1]<-1
    p<-predict(mod,newdata=pdat,type='response',se.fit=TRUE)
    pdat$p<-p$fit
    return(pdat)
} else NULL
}
ot<-ddply(rvw,.(cell),.fun=f,.progress='text')

a<-subset(ot,p==0)
hist(ot$p,breaks=100,col='black')

plot(log(ot$p+.01),(ot$pb))
plot(log(ot$p+.01),log(ot$pb))


f<-function(d){ return(data.frame(n=length(unique(d$id)),
                                  lon=unique(d$lonc.1),
                                  lat=unique(d$latc.1),
                                  totno=sum(d$totno)))}
dt<-ddply(rvw,.(cell.1),.fun=f,.progress='text')
dt$p<-dbinom(1,dt$n,.02)

#dt<-ddply(rvw,.(cell,tbin20),.fun=f,.progress='text')
plot(log10(dt$n),log10(dt$totno+1),pch=16)
cor(log10(dt$n),log10(dt$totno+1))

hist(dt$n,breaks=100,col='black')


adat<-subset(dt,select=c('p','lon','lat'))
adat<-subset(dt,select=c('n','lon','lat'))
ttl<-'dan'
ct<-30
dg<-2


nm<-names(adat)[1]
names(adat)[1]<-'y'
adat$y<-ifelse(adat$y>ct,ct,adat$y)
adat$y<-ifelse(adat$y< -ct,-ct,adat$y)
adat$y<-round(adat$y,digits=2)
a<-adat
aa<-data.frame(y=seq(0,max(abs(adat$y),na.rm=TRUE),length.out=100))
a<-rbind.fill(a,aa)

n<-21
mxx<-max(abs(adat$y))
brks<-seq(0,mxx+.01,length.out=n)
brks2<-round(seq(0,mxx+.01,length.out=n),digits=dg)
a$ycat<-cut(a$y,breaks=brks)
lbls<-sort(unique(a$ycat))
lbls2<-sort(unique(cut(a$y,breaks=brks2)))
cls<-matlab.like(length(lbls))
#cls<-colorRampPalette(c('magenta4','blue3','green','yellow','red3'))
#cls<-(cls(length(lbls)))
#cls<-colorRampPalette(c('dodgerblue4','white','firebrick4'))
#cls<-cls(length(lbls))

return(
ggplot()+
geom_tile(data=a, aes(x=lon, y=lat,fill=ycat),col='gray80',size=.0001)    +
scale_fill_manual(breaks=as.character(lbls),values=cls,labels=lbls2,na.value="transparent",guide=guide_legend(title=paste(ttl)))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.1, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-70,-56,1),labels=as.character(seq(-70,-56,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,48,1),labels=as.character(seq(41,48,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(min(adat$lat)-.3,max(adat$lat)+.3),xlim=c(min(adat$lon)-.3,max(adat$lon)+.3))+
    xlab('')+
    ylab('')
)



pltfun<-function(adat,ttl,ct,dg){

nm<-names(adat)[1]
names(adat)[1]<-'y'
adat$y<-ifelse(adat$y>ct,ct,adat$y)
adat$y<-ifelse(adat$y< -ct,-ct,adat$y)
adat$y<-round(adat$y,digits=2)
a<-adat
aa<-data.frame(y=seq(-(max(abs(adat$y),na.rm=TRUE)),max(abs(adat$y),na.rm=TRUE),length.out=100))
a<-rbind.fill(a,aa)

n<-21
mxx<-max(abs(adat$y))
brks<-seq(-mxx-.01,mxx+.01,length.out=n)
brks2<-round(seq(-mxx-.01,mxx+.01,length.out=n),digits=dg)
a$ycat<-cut(a$y,breaks=brks)
lbls<-sort(unique(a$ycat))
lbls2<-sort(unique(cut(a$y,breaks=brks2)))
#cls<-matlab.like(length(lbls))
#cls<-colorRampPalette(c('magenta4','blue3','green','yellow','red3'))
#cls<-(cls(length(lbls)))
cls<-colorRampPalette(c('dodgerblue4','white','firebrick4'))
cls<-cls(length(lbls))



return(
ggplot()+
geom_tile(data=a, aes(x=lon, y=lat,fill=ycat),col='gray80',size=.0001)    +
scale_fill_manual(breaks=as.character(lbls),values=cls,labels=lbls2,na.value="transparent",guide=guide_legend(title=paste(ttl)))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.1, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-70,-56,1),labels=as.character(seq(-70,-56,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,48,1),labels=as.character(seq(41,48,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(min(adat$lat)-.3,max(adat$lat)+.3),xlim=c(min(adat$lon)-.3,max(adat$lon)+.3))+
    xlab('')+
    ylab('')
)
}
p1<-pltfun(subset(phen2,select=c('delta.mxfall','lon','lat')),'mxfall',2.8,3)


return(ggplot()+
geom_polygon(aes(long,lat,group=group,fill=as.factor(clusters),alpha=cx),data=mydat,col='black',size=.0001) +
    coord_equal()+
    scale_fill_manual(values=c(dum[1],dum[2],dum[3]))+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='bottomright',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(42,47,1),labels=as.character(seq(42,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(42,47),xlim=c(-68,-57))+
  labs(title = lbl,
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab(''))
}


setwd(figsdir)
pdf('herring_rvsurvey_cluster_all_k3.pdf',height=10,width=8)
mp1<-clfun(qsm,3,"Total observed herring (weight)")





summary(rvw$time)
sort(unique(rvw$time))
plot(rvw$time,log10(rvw$no+.1))
mod<-gam(totwgt ~ s(time,bs='cc'),data=rvw,family='nb')
pdat<-data.frame(time=seq(min(rvw$time),max(rvw$time),length.out=100))
p<-predict(mod,newdata=pdat,type='response')
pdat$p<-p
plot(pdat$time,pdat$p)

f<-function(d){    return(data.frame(nyear=length(unique(d$year)),
                                     myear=min(d$year)))
               }
dm<-ddply(rvw,.(strat),.fun=f)
dmm<-dm[order(dm$nyear),]
dmm<-subset(dmm,nyear>=20 & myear<1990)




f<-function(d){
names(d)[1]<-'y'
mod<-gam(y~as.factor(year),data=d,gamma=.5,gamily=negbin)
pdat<-data.frame(year=sort(unique(d$year)),
               ntows=tapply(d$id,d$year,function(x) length(unique(x))))
pdat$p<-predict(mod,newdata=pdat,type='response')
return(pdat)
}
#TOTAL NUMBERS
qq<-ddply(subset(rvw,strat %in% dmm$strat,select=c('totwgt','strat','year','id')),.(strat),.fun=f)
par(mfrow=c(2,2))
plot(qq$ntows,qq$p,col=alpha('darkred',.4),pch=16,cex=2)
plot(log10(qq$ntows),log10(qq$p+1),col=alpha('darkred',.4),pch=16,cex=2)
cor(log10(qq$ntows),log10(qq$p+1))#.019



f<-function(d){
names(d)[1]<-'y'
mod<-gam(y~as.factor(year),data=d,gamma=.5,gamily=binomial)
pdat<-data.frame(year=sort(unique(d$year)),
               ntows=tapply(d$id,d$year,function(x) length(unique(x))))
pdat$p<-predict(mod,newdata=pdat,type='response')
return(pdat)
}
#TOTAL NUMBERS
qq<-ddply(subset(rvw,strat %in% dmm$strat,select=c('pres','strat','year','id')),.(strat),.fun=f)
par(mfrow=c(2,2))
plot(qq$ntows,qq$p,col=alpha('darkred',.4),pch=16,cex=2)
plot(log10(qq$ntows),log10(qq$p+.1),col=alpha('darkred',.4),pch=16,cex=2)
cor(log10(qq$ntows),log10(qq$p+.1))#-.1
plot(qq$year,log10(qq$p+.1),col=alpha('darkred',.4),pch=16,cex=2)
xyplot(log10(p+1)~year|strat,data=qq,pch=15,col=alpha('darkred',.5))


#ESTIMATES SMOOTH TREND IN DIFFERENT QUANTITIES OVER TIME FOR CLUSTER
f<-function(d){
names(d)[1]<-'y'
#mod<-gam(y~as.factor(year),data=d,gamma=.5,gamily=Gamma('log'))
mod<-gam(y~as.factor(year),data=d,gamma=.5,family=nb)
pdat<-data.frame(year=sort(unique(d$year)),
               ntows=tapply(d$id,d$year,function(x) length(unique(x))),
               n=tapply(d$y,d$year,sum),
               sz=tapply(d$sz,d$year,function(x) mean(x,na.rm=TRUE)))
pdat$p<-predict(mod,newdata=pdat,type='response')
#pdat$p<-(pdat$p-mean(pdat$p))/sd(pdat$p)
return(pdat)
}
#TOTAL NUMBERS
qq<-ddply(subset(rvw,strat %in% dmm$strat,select=c('totno','strat','year','id','sz')),.(strat),.fun=f)
plot(log10(qq$sz+.1),log10(qq$p+.01),pch=15)
cor(log10(qq$sz+.1),log10(qq$p+.01),use='pairwise.complete.obs')

plot(log10(qq$ntows),log10(qq$n+1),col=alpha('darkred',.3),pch=16,cex=2)
plot(log10(qq$ntows),log10(qq$p+.01),col=alpha('darkred',.3),pch=16,cex=2)
cor(log10(qq$ntows),log10(qq$p+.01))#.10

plot((qq$ntows),log10(qq$n+.01),col=alpha('darkred',.4),pch=16,cex=2)
plot(log10(qq$ntows),log10(qq$n+.01),col=alpha('darkred',.4),pch=16,cex=2)
cor(log10(qq$ntows),log10(qq$n+.01))#.40


f<-function(d){
return(data.frame(r=cor(log10(d$ntows),log10(d$p+.01),use='pairwise.complete.obs'),
                  r2=cor(log10(d$ntows),log10(d$n+.01),use='pairwise.complete.obs')))
    }
dt<-ddply(qq,.(strat),.fun=f)
dt$dir<-ifelse(dt$r>0,1,-1)
dt$dir2<-ifelse(dt$r2>0,1,-1)

cor(log10(qq$ntows),log10(qq$p+1))
xyplot(log10(p+.1)~log10(ntows) |strat,data=qq,col=alpha('dodgerblue3',.7),pch=15, type=c('p','r'))
xyplot(log10(p+.01)~log10(ntows) |strat,data=qq,col=alpha('dodgerblue3',.7),pch=15, type=c('p','r'))
xyplot(log10(p+1)~log10(ntows) |strat,data=qq,col=alpha('dodgerblue3',.7),pch=15, type=c('p','r'))




#ESTIMATES SMOOTH TREND IN DIFFERENT QUANTITIES OVER TIME FOR CLUSTER
fsm<-function(d){
names(d)[1]<-'y'
mod<-gam(y~s(year),data=d,gamma=.5,gamily=Gamma('log'))
pdat<-data.frame(year=seq(min(d$year),max(d$year),.25))
pdat$p<-predict(mod,newdata=pdat,type='response')

pdat$p<-(pdat$p-mean(pdat$p))/sd(pdat$p)
names(pdat)[2]<-unique(as.character(d$strat))
return(pdat)
}

#TOTAL WEIGHT
qq<-dlply(subset(rvw,strat %in% dmm$strat,select=c('totwgt','strat','year')),.(strat),.fun=fsm)
qsm<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), qq)#COMBINE
qsm<-qsm[,colSums(is.na(qsm)) != nrow(qsm)]#REMOVES COLUMNS THAT ARE ALL MISSING

#TOTAL NUMBERS
qq<-dlply(subset(rvw,strat %in% dmm$strat,select=c('totno','strat','year')),.(strat),.fun=fsm)
qsm2<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), qq)#COMBINE
qsm2<-qsm2[,colSums(is.na(qsm2)) != nrow(qsm2)]#REMOVES COLUMNS THAT ARE ALL MISSING

#AVERAGE SIZE
qq<-dlply(subset(rvw,strat %in% dmm$strat,select=c('sz','strat','year')),.(strat),.fun=fsm)
qsm3<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), qq)#COMBINE
qsm3<-qsm3[,colSums(is.na(qsm3)) != nrow(qsm3)]#REMOVES COLUMNS THAT ARE ALL MISSING


q<-qsm3
q2<- q %>% gather(strata, value, -year)
xyplot(value ~ year | strata,data=q2, pch=15, type=c('spline'),col='black')
xyplot(value ~ year | strata,data=q2, pch=15, type=c('p','spline'),col='black')






clfun<-function(df,k,lbl){
rownames(df)<-df$year
df<-df[,-1]
#k<-3#ER OF CLUSTERS
dmat<-1-cor(df,use='pairwise.complete.obs')
dst<-as.dist(dmat)
ff<-fanny(dst,k,maxit=5000,diss=T)

par(mfrow=c(2,1),mar=c(4,4,1,1))
dum<-c('red3','forestgreen','darkblue','cornflowerblue','darkblue')
plot(silhouette(ff),col=dum[1:k],main='')#silhouette plot

dc.pcoa<-cmdscale(dst)
dc.scores<-scores(dc.pcoa,choices=c(1,2))
spefuz.g<-ff$clustering
a<-data.frame(strat=as.character(sort(unique(names(df)))),
              clusters=ff$clustering)
aa<-data.frame(ff$membership)
aa$strat<-rownames(a)

#par(mar=c(1,1,1,8),oma=c(1,1,1,1))
plot(scores(dc.pcoa),asp=1,type='n',xlim=c(-1.5,1),ylim=c(-1,1.2),las=1,axes=TRUE,xlab='',ylab='')
stars(ff$membership,location=scores(dc.pcoa),draw.segments=T,add=T,scale=F,len=.1,col.segments=alpha(c(dum[1:k]),.25),byt='n',labels=NULL,xlim=c(-1.1,1.4),ylim=c(-1,.5),lwd=.0001,xpd=TRUE,border=NULL,radius=FALSE,col.radius=alpha('white',.1))
for(i in 1:k){ cl<-dum[i]
    gg<-dc.scores[spefuz.g==i,]
    hpts<-chull(gg)
    hpts<-c(hpts,hpts[1])
    lines(gg[hpts,],col=cl,lwd=3,xlim=c(-1.1,1.4),ylim=c(-1,.5))
}

cx <- data.frame(strat=aa$strat,
                 cx=apply(aa[, 1:k], 1, max))
cx$cxx <- rescale(cx$cx,newrange=c(.05,.5))
#PLOTS THE TIMESERIES OF R FOR EACH CLUSTER
par(mfrow=c(3,1),mar=c(2,12,1,12),oma=c(1,1,1,1))

mb<-seq(1,k,1)
l<-list()
for(i in 1:length(mb)){
print(mb[i])
one<-subset(a,clusters==mb[i],select=c('strat'))# IN THIS CLUSTER
cx2<-subset(cx,cx>=0.5)#GETS INSTANCES WHERE CLUSTER PROBABILITY>0.5
data5<-subset(df,select=c(as.character(one$strat)))#TS FOR CLUSTER
cx2<-subset(cx2,strat %in% names(data5))
data5<-subset(df,select=c(as.character(cx2$strat)))#TS FOR CLUSTER
x<-rownames(data5)
t<-data.frame(year=as.numeric(x),mn=rowMeans(data5,na.rm=F))
t2<-data.frame(year=as.numeric(x),mn=rowMeans(data5,na.rm=T))

cl<-dum[i]
plot(0,0,pch=16,cex=.01,xlim=c(1970,2015),ylim=c(-2,3),main='',xaxt='n',las=1,axes=FALSE,xlab='',ylab='')
axis(side=2,at=seq(-2,3,1),las=1,lwd=.001,cex.axis=.75)
axis(side=1,at=seq(1970,2020,10),cex.axis=.75)
    for(j in 1:length(data5[1,])){
        try(dat<-data5[,j])
        try(trnsp<-subset(cx,strat==as.character(names(data5[j])))$cxx)
        try(lines(x,dat,cex=.5,ylim=c(0,1),lwd=2,col=alpha(cl,rescale(trnsp,newrange=c(.1,.5)))))
    }
lines(as.numeric(as.character(t$year)),t$mn,col='gold3',cex=.7,lwd=3)
dm<-subset(t2,mn==max(t2$mn,na.rm=TRUE),select=c('year'))
dm$cluster=mb[i]
names(dm)<-c('clusterday','clusters')
l[[i]]<-dm
}


a2<-a
if(k==3){a2$cx<- apply(aa[,c('X1','X2','X3')], 1, function(x) max(x) )
     } else {a2$cx<- apply(aa[,c('X1','X2')], 1, function(x) max(x) )
     }
a2$id<-a2$strat
am<-fortify(plg,region='stratum')
am<-subset(am,id%in% a2$id)
mydat<-merge(am,a2,by=c('id'))
mydat$cx<-rescale(mydat$cx,newrange=c(.4,1))

par(mfrow=c(1,1),mar=c(3,3,3,3))
return(ggplot()+
geom_polygon(aes(long,lat,group=group,fill=as.factor(clusters),alpha=cx),data=mydat,col='black',size=.0001) +
    coord_equal()+
    scale_fill_manual(values=c(dum[1],dum[2],dum[3]))+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='bottomright',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(42,47,1),labels=as.character(seq(42,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(42,47),xlim=c(-68,-57))+
  labs(title = lbl,
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab(''))
}


setwd(figsdir)
pdf('herring_rvsurvey_cluster_all_k3.pdf',height=10,width=8)
mp1<-clfun(qsm,3,"Total observed herring (weight)")
mp2<-clfun(qsm2,3,"Total observed herring (numb)")
mp3<-clfun(qsm3,3,"Average weight")
dev.off()

pdf('herring_rvsurvey_cluster_map_k3.pdf',height=14,width=8)
grid.arrange(mp1,mp2,mp3,ncol=1)
dev.off()


setwd(figsdir)
pdf('herring_rvsurvey_cluster_all_k2.pdf',height=10,width=8)
mp1<-clfun(qsm,2,"Total observed herring (weight)")
mp2<-clfun(qsm2,2,"Total observed herring (numb)")
mp3<-clfun(qsm3,2,"Average weight")
dev.off()

pdf('herring_rvsurvey_cluster_map_k2.pdf',height=14,width=8)
grid.arrange(mp1,mp2,mp3,ncol=1)
dev.off()

















###################################################

#TOTAL NUMBER OF HERRING RECORDED FROM ALL TRAWLS ALL YEARS
f<-function(d){
    return(data.frame(wt=mean(d$totwgt),
                 no=mean(d$totno),
                 dep=mean(d$depth,na.rm=TRUE),
                 sz=mean(d$sz)))
}
d1<-ddply(rvw,.(strat),.fun=f)
d1<-d1[order(d1$sz,decreasing=TRUE),]

cor(subset(d1,select=c('wt','no','dep','sz')),use='pairwise.complete.obs')
plot(subset(d1,select=c('wt','no','dep','sz')),pch=15)

a2<-d1
a2$id<-a2$strat
am<-fortify(plg,region='stratum')
am<-subset(am,id%in% a2$id)
mydat<-merge(am,a2,by=c('id'))

pwt<-ggplot()+
geom_polygon(aes(long,lat,group=group,fill=wt),data=mydat,col='black',size=.0001) +
    coord_equal()+
scale_fill_distiller(palette='Spectral')+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='bottomright',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,47,1),labels=as.character(seq(41,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,47),xlim=c(-68,-57))+
  labs(title = "Total observed herring (weight)",
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab('')

pno<-ggplot()+
geom_polygon(aes(long,lat,group=group,fill=no),data=mydat,col='black',size=.0001) +
    coord_equal()+
scale_fill_distiller(palette='Spectral')+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='bottomright',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,47,1),labels=as.character(seq(41,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,47),xlim=c(-68,-57))+
  labs(title = "Total observed herring (numbers)",
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab('')


psz<-ggplot()+
geom_polygon(aes(long,lat,group=group,fill=sz),data=mydat,col='black',size=.0001) +
    coord_equal()+
scale_fill_distiller(palette='Spectral')+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='bottomright',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,47,1),labels=as.character(seq(41,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,47),xlim=c(-68,-57))+
  labs(title = "Average weight of all observed herring",
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab('')

setwd(figsdir)
pdf('herring_rvsurvey_allobs_map2.pdf',height=14,width=8)
grid.arrange(pwt,pno,psz,ncol=1)
dev.off()



rvw$dayc<-cut(rvw$day,breaks=seq(170,230,10),labels=seq(175,225,10))

f<-function(d){
    if(length(unique(d$dayc))>1){
    mod<-gam(totwgt~as.factor(dayc),data=d,gamma=1)
    pdat<-data.frame(dayc=sort(unique(d$dayc)))
    pdat$p<-predict(mod,newdata=pdat)
    return(pdat)
} else NULL
}
phdat<-ddply(rvw,.(strat),.fun=f)


a2<-phdat
a2$id<-a2$strat
am<-fortify(plg,region='stratum')
am<-subset(am,id%in% a2$id)
mydat<-merge(am,a2,by=c('id'))
mydat$p<-(mydat$p-mean(mydat$p))/sd(mydat$p)

dys<-sort(unique(phdat$dayc))
l<-list()
for(i in 1:length(dys)){

d<-subset(mydat,dayc==dys[i])
d$p<-log10(d$p+1)
print(dim(d))
l[[i]]<-ggplot()+
geom_polygon(aes(long,lat,group=group,fill=p),data=d,col='black',size=.0001) +
    coord_equal()+
scale_fill_distiller(palette='Spectral')+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.8,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,47,1),labels=as.character(seq(41,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,47),xlim=c(-68,-57))+
  labs(title = paste("Total observed herring (day=",dys[i]),
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab('')
}
setwd(figsdir)
pdf('herring_rvsurvey_seasonal_daysbin.pdf',height=14,width=10)
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],ncol=2)
dev.off()


a<-subset(mydat,dayc==205)
hist(a$p,breaks=50)
hist(log(a$p),breaks=50)



d<-subset(rvw,strat==462)
a2<-a2[order(a2$wt,decreasing=TRUE),]

f<-function(d){
    if(length(unique(d$month))>1){
    mod<-gam(totwgt~as.factor(month),data=d,gamma=1)
    pdat<-data.frame(month=sort(unique(d$month)))
    pdat$p<-predict(mod,newdata=pdat)
    return(pdat)
} else NULL
}
phdat<-ddply(rvw,.(strat),.fun=f)
plot(rvw$day,rvw$totno)


a2<-phdat
a2$id<-a2$strat
am<-fortify(plg,region='stratum')
am<-subset(am,id%in% a2$id)
mydat<-merge(am,a2,by=c('id'))

mydat1<-subset(mydat,month==6)
p1<-ggplot()+
geom_polygon(aes(long,lat,group=group,fill=p),data=mydat1,col='black',size=.0001) +
    coord_equal()+
scale_fill_distiller(palette='Spectral')+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='bottomright',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,47,1),labels=as.character(seq(41,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,47),xlim=c(-68,-57))+
  labs(title = "Total observed herring (June)",
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab('')

mydat2<-subset(mydat,month==7)
p2<-ggplot()+
geom_polygon(aes(long,lat,group=group,fill=p),data=mydat2,col='black',size=.0001) +
    coord_equal()+
scale_fill_distiller(palette='Spectral')+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='bottomright',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,47,1),labels=as.character(seq(41,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,47),xlim=c(-68,-57))+
  labs(title = "Total observed herring (July)",
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab('')

mydat3<-subset(mydat,month==8)
p3<-ggplot()+
geom_polygon(aes(long,lat,group=group,fill=p),data=mydat3,col='black',size=.0001) +
    coord_equal()+
scale_fill_distiller(palette='Spectral')+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='bottomright',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,47,1),labels=as.character(seq(41,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,47),xlim=c(-68,-57))+
  labs(title = "Total observed herring (August)",
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab('')


setwd(figsdir)
pdf('herring_rvsurvey_bymonth_map.pdf',height=14,width=8)
grid.arrange(p1,p2,p3,ncol=1)
dev.off()












rvw$yearc<-cut(rvw$year,breaks=seq(1970,2020,5),labels=seq(1972.5,2017.5,5),include.lowest=TRUE)

f<-function(d){
    return(data.frame(no=median(d$totno),
                      wt=median(d$totwgt),
                      sz=median(d$sz)))
}
d1<-ddply(rvw,.(strat,yearc),.fun=f)

a2<-d1
a2$id<-a2$strat
am<-fortify(plg,region='stratum')
am<-subset(am,id%in% a2$id)
mydat<-merge(am,a2,by=c('id'))
#mydat$p<-(mydat$p-mean(mydat$p))/sd(mydat$p)

yrs<-sort(unique(mydat$yearc))
l<-list()
for(i in 1:length(yrs)){

d<-subset(mydat,yearc==yrs[i])
d$wt<-log10(d$wt)
print(dim(d))
l[[i]]<-ggplot()+
geom_polygon(aes(long,lat,group=group,fill=wt),data=d,col='black',size=.0001) +
    coord_equal()+
scale_fill_distiller(palette='Spectral')+
    geom_polygon(aes(long,lat,group=group),fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.8,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.11, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,47,1),labels=as.character(seq(41,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,47),xlim=c(-68,-57))+
  labs(title = paste("Total observed herring (year=",yrs[i]),
       subtitle = "",
       caption = 'Data source: Summer RV survey') +
    xlab('')+
    ylab('')
}
setwd(figsdir)
pdf('herring_rvsurvey_wt_byyear.pdf',height=16,width=16)
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]],l[[8]],l[[9]],l[[10]],ncol=3)
dev.off()




f<-function(d){
if(length(unique(d$year))>=20){
print(unique(d$strat))
print(length(unique(d$year)))

mod<-gam(no~year,data=d,gamma=1)
s<-summary(mod)
no.pv<-s$p.table[2,4]
no.b<-(10^s$p.table[2,1])-1

mod<-gam(wgt~year,data=d,gamma=1)
s<-summary(mod)
wt.pv<-s$p.table[2,4]
wt.b<-(10^s$p.table[2,1])-1

d$depth<-(d$dmin+d$dmax)/2
d$totno<-ifelse(d$totno==0,1,d$totno)
d$M<-d$totwgt/d$totno
d$temperature<-(d$surface_temperature+d$bottom_temperature)/2
d<-subset(d,is.na(temperature)==FALSE)
d$kelvins<-d$temperature+273.15
k<-0.00008617#BOLTZMANN CONSTANT
d$metai<-(d$M^.25)*(exp(-1*(1/(k*d$kelvins))))
d$metai<-d$metai*10^18
mod<-gam(metai~as.factor(year),data=d,gamma=1)
s<-summary(mod)
met.pv<-s$p.table[2,4]
met.b<-(10^s$p.table[2,1])-1

return(data.frame(lon=mean(d$lon),
                  lat=mean(d$lat),
                  no.b=no.b,
                  no.pv=no.pv,
                  wt.b=wt.b,
                  wt.pv=wt.pv,
                  met.b=met.b,
                  met.pv=met.pv))
} else NULL
}
sdat1<-ddply(rvw,.(strat),.fun=f)


map('world',xlim=c(-70,-60),ylim=c(42,46),fill=TRUE,col='gray')
map.axes()
points(sdat1$lon,sdat1$lat,pch=16,cex=rescale(abs(sdat1$no.b),newrange=c(3,14)),col=ifelse(sdat1$no.b>0,'firebrick3','dodgerblue3'))

map('world',xlim=c(-70,-60),ylim=c(42,46),fill=TRUE,col='gray')
map.axes()
points(sdat1$lon,sdat1$lat,pch=16,cex=rescale(abs(sdat1$wt.b),newrange=c(3,14)),col=ifelse(sdat1$wt.b>0,'firebrick3','dodgerblue3'))

map('world',xlim=c(-70,-60),ylim=c(42,46),fill=TRUE,col='gray')
map.axes()
points(sdat1$lon,sdat1$lat,pch=16,cex=rescale(abs(sdat1$met.b),newrange=c(3,14)),col=ifelse(sdat1$met.b>0,'firebrick3','dodgerblue3'))

cor(sdat1$no.b,sdat1$wt.b)
cor(sdat1$no.b,sdat1$met.b)
plot(sdat1$met.b,sdat1$no.b,pch=15)


setwd(datadir)
rvl<-read.csv("herring_lengths_RV_survey_spera_spawnar.csv",header=TRUE)
rvl$flen<-ifelse(rvl$mission=='NED2016016',rvl$flen/10,rvl$flen)
rvl$flen<-ifelse(rvl$mission!='NED2016016',(rvl$flen*1.0866)+0.95632,rvl$flen)
#plot(log10(rvl$flen),log10(rvl$fwt))
rvl<-subset(rvl,log10(rvl$flen)>=.5 & month %in% c(6,7,8))#REMOVE OUTLIERS

rvl2<-rvl
rvl2$lcat<-cut(rvl$flen,breaks=seq(5,45,10),labels=seq(10,40,10))
rvl2$id<-gsub(' ','',paste(rvl2$mission,'_',rvl2$setno))

f<-function(d){
return(data.frame(year=unique(d$year),
                  lon=unique(d$lon),
                  lat=unique(d$lat),
                  strat=unique(d$strat),
                  no=sum(d$clen,na.rm=TRUE),
                  wt=sum(d$fwt/1000,na.rm=TRUE)))
}
rvll<-ddply(rvl2,.(id,lcat),.fun=f,.progress='text')

d<-subset(rvll,strat==sort(unique(rvll$strat))[1] & lcat==27.5)

f<-function(d){
if(length(unique(d$year))>=25){
print(unique(d$strat))
print(unique(d$lcat))

mod<-gam(no~year,data=d,gamma=1.4)
s<-summary(mod)
no.pv<-s$p.table[2,4]
no.b<-s$p.table[2,1]

return(data.frame(lon=mean(d$lon),
                  lat=mean(d$lat),
                  no.b=no.b,
                  no.pv=no.pv))
} else NULL
}
sdat2<-ddply(rvll,.(strat,lcat),.fun=f,.progress='text')

a<-subset(sdat2,lcat==20)
a<-subset(sdat2,lcat==30)
a<-subset(sdat2,lcat==40)
map('world',xlim=c(-70,-60),ylim=c(42,46),fill=TRUE,col='gray')
map.axes()
points(a$lon,a$lat,pch=16,cex=rescale(abs(a$no.b),newrange=c(3,14)),col=alpha(ifelse(a$no.b>0,'firebrick3','dodgerblue3'),.75))
















###################################################################

#################   PHENOLOGY

###################################################################



##################################################################
setwd(datadir)
rvw<-read.csv("herring_weights_RV_survey_spera_allar.csv",header=TRUE)
rvw$strat<-as.character(rvw$strat)
#rvw<-subset(rvw,!(strat %in% c('493','494','')))

rvw$bank<-ifelse(rvw$strat %in% c(447,448),'banq','no')
rvw$bank<-ifelse(rvw$strat %in% c(443),'mis',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(458),'mid',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(455,456),'sab',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(464),'west',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(463),'em',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(473),'lh',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(474),'rw',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(475),'bac',rvw$bank)
rvw$bank<-ifelse(rvw$strat %in% c(480),'bn',rvw$bank)

rvw$no<-log10(rvw$totno+1)
rvw$wgt<-log10(rvw$totwgt+1)
rvw$sz<-rvw$totwgt/rvw$totno
rvw$sz<-ifelse(rvw$sz==Inf,rvw$totwgt,rvw$sz)

rvw$id<-gsub(' ','',paste(rvw$mission,'_',rvw$setno))
rvw$pres<-ifelse(rvw$totno>0,1,0)
rvw$tbin20<-ifelse(rvw$year<=1990,1980,2010)

rvw$lonc<-round(rvw$lon,digits=0)
rvw$lonc<-ifelse(rvw$lon<=rvw$lonc,rvw$lonc-.25,rvw$lonc+.25)
rvw$latc<-round(rvw$lat,digits=0)
rvw$latc<-ifelse(rvw$lat>=rvw$latc,rvw$latc+.25,rvw$latc-.25)
rvw$cell<-gsub(' ','',paste(rvw$lonc,'_',rvw$latc))

rvw$cell.1<-gsub(' ','',paste(round(rvw$lon,digits=1),'_',round(rvw$lat,digits=1)))
lonc.1<-seq(min(rvw$lon),max(rvw$lon),.1)
latc.1<-seq(min(rvw$lat),max(rvw$lat),.1)
crds<-expand.grid(lonc.1=lonc.1,latc.1=latc.1)
crds$cell.1<-gsub(' ','',paste(round(crds$lonc.1,digits=1),'_',round(crds$latc.1,digits=1)))
rvw<-merge(rvw,crds,by=c('cell.1'),all.x=TRUE,all.y=FALSE)


#ESTIMATES SMOOTH TREND IN PHENOLOGY BY STRATA - SHORT FORMAT
rvw2<-subset(rvw,!(geardesc %in% c('Newston net','Campelen 1800 survey trawl')))
f<-function(d){
if(dim(subset(d,month %in% c(1,2,3)))[1]>10){wint<-1
                                         } else { wint<-0 }
if(dim(subset(d,month %in% c(6,7,8)))[1]>10){summ<-1
                                         } else { summ<-0 }
if(dim(subset(d,month %in% c(9,10,11)))[1]>10){fall<-1
                                         } else { fall<-0 }
    return(data.frame(nmo=length(unique(d$month)),
                      ndy=length(unique(d$day)),
                      nyr=length(unique(d$year)),
                      wint=wint,
                      summ=summ,
                      fall=fall,
                      tot=wint+summ+fall,
                      lon=mean(d$lon),
                      lat=mean(d$lat)))
}
dt<-ddply(rvw2,.(strat),.fun=f)
dt<-dt[order(dt$nmo,decreasing=TRUE),]
dmm<-subset(dt,nmo>= 4 & ndy>=30 & nyr>=3 & tot>=3)




f<-function(d){
if(dim(subset(d,month %in% c(1,2,3)))[1]>10){wint<-1
                                         } else { wint<-0 }
if(dim(subset(d,month %in% c(6,7,8)))[1]>10){summ<-1
                                         } else { summ<-0 }
if(dim(subset(d,month %in% c(9,10,11)))[1]>10){fall<-1
                                         } else { fall<-0 }
    return(data.frame(nmo=length(unique(d$month)),
                      ndy=length(unique(d$day)),
                      nyr=length(unique(d$year)),
                      wint=wint,
                      summ=summ,
                      fall=fall,
                      tot=wint+summ+fall,
                      lon=mean(d$lon),
                      lat=mean(d$lat)))
}
dt<-ddply(rvw2,.(cell),.fun=f)
dt<-dt[order(dt$nmo,decreasing=TRUE),]
dmm<-subset(dt,nmo>= 4 & ndy>=30 & nyr>=3 & tot>=3)


f1<-function(d0){
f<-function(d){
day<-unique(d$day)
d2<-subset(d0,select=c('day'))
d2$day2<-d2$day-(day-1)
d2$day2<-ifelse(d2$day2<=0,d2$day2+(365),d2$day2)
    return(data.frame(day=day,
                   span=diff(range(min(d2$day2),max(d2$day2)))))
}
z<-ddply(d0,.(day),.fun=f)
z<-subset(z,span==max(z$span))
return(subset(z,day==min(z$day)))
}
#dayshift<-ddply(unique(subset(rvw2,strat %in% dmm$strat,select=c('strat','day'))),.(strat),.fun=f1,.progress='text')
dayshift<-ddply(unique(subset(rvw2,cell %in% dmm$cell,select=c('strat','day','cell'))),.(cell),.fun=f1,.progress='text')

##############################################
#ADJUSTS DAY VALUES TO MAXIMIZE PHENOLOGY SPAN
dayshiftfun<-function(d,shf){
d$day2<-d$day-(shf)
d$day2<-ifelse(d$day2<=0,d$day2+(365),d$day2)
return(d)
}

#BACK-CALCULATE ORIGINAL DAY OF THE YEAR
revdayshiftfun<-function(d,shf){
d$day3<-d$day2+shf
d$day3<-ifelse(d$day3>365,d$day3-365,d$day3)
return(d)
}

ff<-function(d){
shf<-subset(dayshift,cell==unique(d$cell))$day
print(shf)
names(d)[1]<-'y'
#SHIFT DAYS TO CENTER ON JULY AND MAKE CONTINUOUS
d<-dayshiftfun(d,shf)
mod<-gam(y~s(year) + s(time,k=4,bs='cc') + s(day2,k=4,bs='cc'),data=d,gamma=1.4,family='nb'(link='log'))
    pdat<-data.frame(year=max(d$year),
                     time=1200,
                     day2=seq(min(d$day2),max(d$day2),1))
    p<-predict(mod,newdata=pdat,type='response')
    pdat$p<-p
    pdat$pz<-(pdat$p-mean(pdat$p))/sd(pdat$p)
    pdat<-revdayshiftfun(pdat,shf)
    pdat<-subset(pdat,select=c('day3','p','pz'))
pdat$lon<-median(d$lon)
pdat$lat<-median(d$lat)
pdat$did<-ifelse(pdat$day %in% d$day,1,0)
names(pdat)[1]<-'day'
return(data.frame(pdat))
}
#TOTAL WEIGHT
#out<-ddply(subset(rvw2, strat %in% dmm$strat,select=c('totwgt','strat','year','day','time','lon','lat')),.(strat),.fun=ff,.progress='text')
out<-ddply(subset(rvw2, cell %in% dmm$cell,select=c('totwgt','strat','year','day','time','lon','lat','cell')),.(cell),.fun=ff,.progress='text')
xyplot(pz~day|strat,data=out,col=ifelse(out$did==1,'red','blue'))
xyplot(pcx~day|strat,data=out,col=ifelse(out$did==1,'red','blue'))







###########################################################
#ESTIMATE PHENOLOGY FOR EACH 1/4 GRID CELL AND ANIMATE
rvw2<-subset(rvw,!(geardesc %in% c('Newston net','Campelen 1800 survey trawl')))
f<-function(d){
if(dim(subset(d,month %in% c(1,2,3)))[1]>10){wint<-1
                                         } else { wint<-0 }
if(dim(subset(d,month %in% c(6,7,8)))[1]>10){summ<-1
                                         } else { summ<-0 }
if(dim(subset(d,month %in% c(9,10,11)))[1]>10){fall<-1
                                         } else { fall<-0 }
    return(data.frame(nmo=length(unique(d$month)),
                      ndy=length(unique(d$day)),
                      nyr=length(unique(d$year)),
                      wint=wint,
                      summ=summ,
                      fall=fall,
                      tot=wint+summ+fall,
                      lon=mean(d$lon),
                      lat=mean(d$lat)))
}
dt<-ddply(rvw2,.(cell),.fun=f)
dt<-dt[order(dt$nmo,decreasing=TRUE),]
dmm<-subset(dt,nmo>= 4 & ndy>=30 & nyr>=3 & tot>=3)


f1<-function(d0){
f<-function(d){
day<-unique(d$day)
d2<-subset(d0,select=c('day'))
d2$day2<-d2$day-(day-1)
d2$day2<-ifelse(d2$day2<=0,d2$day2+(365),d2$day2)
    return(data.frame(day=day,
                   span=diff(range(min(d2$day2),max(d2$day2)))))
}
z<-ddply(d0,.(day),.fun=f)
z<-subset(z,span==max(z$span))
return(subset(z,day==min(z$day)))
}
#dayshift<-ddply(unique(subset(rvw2,strat %in% dmm$strat,select=c('strat','day'))),.(strat),.fun=f1,.progress='text')
dayshift<-ddply(unique(subset(rvw2,cell %in% dmm$cell,select=c('strat','day','cell'))),.(cell),.fun=f1,.progress='text')

##############################################
#ADJUSTS DAY VALUES TO MAXIMIZE PHENOLOGY SPAN
dayshiftfun<-function(d,shf){
d$day2<-d$day-(shf)
d$day2<-ifelse(d$day2<=0,d$day2+(365),d$day2)
return(d)
}

#BACK-CALCULATE ORIGINAL DAY OF THE YEAR
revdayshiftfun<-function(d,shf){
d$day3<-d$day2+shf
d$day3<-ifelse(d$day3>365,d$day3-365,d$day3)
return(d)
}

library(akima)
ff<-function(d){
shf<-subset(dayshift,cell==unique(d$cell))$day
print(shf)
names(d)[1]<-'y'
#SHIFT DAYS TO CENTER ON JULY AND MAKE CONTINUOUS
d<-dayshiftfun(d,shf)
mod<-gam(y~s(year) + s(time,k=4,bs='cc') + s(day2,k=4,bs='cc'),data=d,gamma=1.4,family='nb'(link='log'))
    pdat0<-data.frame(year=2000,
                     time=1200,
                     day2=seq(min(d$day2),max(d$day2),1))
    p<-predict(mod,newdata=pdat0,type='response')
    pdat0$p<-p
pdat<-data.frame(day2=seq(1,365,1))
pdat$p<-pasp<-aspline(pdat0$day2,pdat0$p,xout=pdat$day2)$y
    pdat$pz<-(pdat$p-mean(pdat$p))/sd(pdat$p)
    pdat<-revdayshiftfun(pdat,shf)
    pdat<-subset(pdat,select=c('day3','p','pz'))
pdat$lon<-unique(d$lonc)
pdat$lat<-unique(d$latc)
pdat$depth<-mean(d$dmax,na.rm=TRUE)
pdat$did<-ifelse(pdat$day %in% d$day,1,0)
names(pdat)[1]<-'day'
return(data.frame(pdat))
}
#TOTAL WEIGHT
out<-ddply(subset(rvw2, cell %in% dmm$cell,select=c('totwgt','strat','year','day','time','lon','lat','cell','lonc','latc','dmax')),.(cell),.fun=ff,.progress='text')
xyplot(pz~day|cell,data=out,col=ifelse(out$did==1,'red','blue'))
xyplot(log10(p+1)~day|cell,data=out,col=ifelse(out$did==1,'red','blue'))
xyplot(p~day|cell,data=out,col=ifelse(out$did==1,'red','blue'))


bnk<-subset(plg,stratum %in% c(443,458,455,456,464,463,473,474,475,480))
x1<--68
x2<--57
y<-41.5
out$dy<-rescale(out$day,newrange=c(x1,x2))
frames <- length(unique(out$day))

rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',i,'plot.jpg',sep=''))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste('00',i,'plot.jpg', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', i,'plot.jpg', sep=''))
  }
}


out$pcx<-rescale(out$pz,newrange=c(.2,7))
setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/Figures/prac')
#loop through plots
for(i in 1:frames){
  name <- rename(i)
  d<-subset(out,day==i)

#saves the plot as a file in the working directory
jpeg(name,width=6,height=5,units='in',quality=100,res=300)
map('worldHires',fill=TRUE,col='gray',border=NA,xlim=c(x1,x2),ylim=c(y,46))
#plot(plg,add=TRUE,lwd=.01,border=alpha('lightgray',.5))
plot(bnk,add=TRUE,lwd=.01,border=NA,col=alpha('firebrick3',.2))
points(d$lon,d$lat,pch=16,cex=d$pcx,col=alpha('dodgerblue3',.3),xlim=c(x1,x2),ylim=c(y,46))
points(d$lon,d$lat,pch=1,lwd=.01,cex=d$pcx,col='dodgerblue3',xlim=c(-x1,x2),ylim=c(y,46))
legend('bottomright',paste('Day=',i),bty='n')
points(unique(d$dy),y,pch=17,col='red3',cex=2)
axis(1,at=seq(x1,x2,length.out=14),labels=c('','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC',''),cex.axis=.7)
dev.off()
}



#################################################################
# SAME BUT FOR NON STANDARDIZED BIOMSSS
out$dy<-rescale(out$day,newrange=c(x1,x2))
out$pcx<-rescale(log10(out$p+1),newrange=c(.2,7))
setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/Figures/prac')
#loop through plots
for(i in 1:frames){
  name <- rename(i)
  d<-subset(out,day==i)

#saves the plot as a file in the working directory
jpeg(name,width=6,height=5,units='in',quality=100,res=300)
map('worldHires',fill=TRUE,col='gray',border=NA,xlim=c(x1,x2),ylim=c(y,46))
#plot(plg,add=TRUE,lwd=.01,border=alpha('lightgray',.5))
plot(bnk,add=TRUE,lwd=.01,border=NA,col=alpha('firebrick3',.2))
points(d$lon,d$lat,pch=16,cex=d$pcx,col=alpha('dodgerblue3',.3),xlim=c(x1,x2),ylim=c(y,46))
points(d$lon,d$lat,pch=1,lwd=.01,cex=d$pcx,col='dodgerblue3',xlim=c(-x1,x2),ylim=c(y,46))
legend('bottomright',paste('Day=',i),bty='n')
points(unique(d$dy),y,pch=17,col='red3',cex=2)
#axis(1,at=seq(x1,x2,length.out=10),labels=round(seq(0,365,length.out=10),digits=0))
axis(1,at=seq(x1,x2,length.out=14),labels=c('','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC',''),cex.axis=.7)
dev.off()
}

xyplot(p~day|cell,data=out,pch=16)
xyplot(log10(p+1)~day|cell,data=out,pch=16)
xyplot(pcx~day|cell,data=out,pch=16)

#FROM COMMAND LINE WITH IMAGEMAGIK INSTALLED RUN: CONVERT -DELAY 2 -QUALITY 100 *.JPG MOVIE.MP4
#run ImageMagick
#library(magick)
#my_command <- 'convert *.png -delay 2 -loop 0 animation.gif'
#system(my_command)

f<-function(d){
return(data.frame(lon=unique(d$lon),
       lat=unique(d$lat),
       depth=unique(d$depth),
       day=subset(d,pz==max(d$pz))$day[1],
       amp=max(d$pz)[1]-min(d$pz)[1]))
}
tm<-ddply(out,.(cell),.fun=f)
plot(subset(tm,select=c('lon','lat','depth','day','amp')),pch=15)

plot(tm$lon,tm$lat,pch=16,cex=rescale(tm$day,newrange=c(.2,6)),col='purple')
plot(tm$lon,tm$lat,pch=16,cex=rescale(tm$amp,newrange=c(.2,6)),col='purple')




###########################################################
#ESTIMATE PHENOLOGY OF HERRING SIZE AND ANIMATE
b<-subset(rvw2,(totwgt==0 & totno>0) | (totwgt>0 & totno==0))
rvw3<-subset(rvw2,!(id %in% b$id))

f<-function(d){
if(dim(subset(d,month %in% c(1,2,3)))[1]>10){wint<-1
                                         } else { wint<-0 }
if(dim(subset(d,month %in% c(6,7,8)))[1]>10){summ<-1
                                         } else { summ<-0 }
if(dim(subset(d,month %in% c(9,10,11)))[1]>10){fall<-1
                                         } else { fall<-0 }
    return(data.frame(nmo=length(unique(d$month)),
                      ndy=length(unique(d$day)),
                      nyr=length(unique(d$year)),
                      wint=wint,
                      summ=summ,
                      fall=fall,
                      tot=wint+summ+fall,
                      lon=mean(d$lon),
                      lat=mean(d$lat)))
}
dt<-ddply(rvw3,.(cell),.fun=f)
dt<-dt[order(dt$nmo,decreasing=TRUE),]
dmm<-subset(dt,nmo>= 4 & ndy>=30 & nyr>=3 & tot>=3)


f1<-function(d0){
f<-function(d){
day<-unique(d$day)
d2<-subset(d0,select=c('day'))
d2$day2<-d2$day-(day-1)
d2$day2<-ifelse(d2$day2<=0,d2$day2+(365),d2$day2)
    return(data.frame(day=day,
                   span=diff(range(min(d2$day2),max(d2$day2)))))
}
z<-ddply(d0,.(day),.fun=f)
z<-subset(z,span==max(z$span))
return(subset(z,day==min(z$day)))
}
dayshift<-ddply(unique(subset(rvw3,cell %in% dmm$cell,select=c('strat','day','cell'))),.(cell),.fun=f1,.progress='text')

##############################################
#ADJUSTS DAY VALUES TO MAXIMIZE PHENOLOGY SPAN
dayshiftfun<-function(d,shf){
d$day2<-d$day-(shf)
d$day2<-ifelse(d$day2<=0,d$day2+(365),d$day2)
return(d)
}

#BACK-CALCULATE ORIGINAL DAY OF THE YEAR
revdayshiftfun<-function(d,shf){
d$day3<-d$day2+shf
d$day3<-ifelse(d$day3>365,d$day3-365,d$day3)
return(d)
}

library(akima)

ff<-function(d){
shf<-subset(dayshift,cell==unique(d$cell))$day
print(unique(d$cell))
print(length(unique(d$day)))
#SHIFT DAYS TO CENTER ON JULY AND MAKE CONTINUOUS
d<-dayshiftfun(d,shf)
modw<-gam(totwgt~s(year) + s(day2,k=4,bs='cc') + s(time,bs='cc',k=4),data=d,gamma=1.4,family='nb'(link='log'))
    pdat0<-data.frame(year=sort(unique(d$year),decreasing=TRUE)[1],
                      time=1200,
                     day2=seq(min(d$day2),max(d$day2),1))
    p<-predict(modw,newdata=pdat0,type='response')
    pdat0$p<-p
pdat<-data.frame(day2=seq(1,365,1))
pdat$p<-pasp<-aspline(pdat0$day2,pdat0$p,xout=pdat$day2)$y
    pdat<-revdayshiftfun(pdat,shf)
    pdatw<-subset(pdat,select=c('day3','p'))
    names(pdatw)<-c('day','pwt')


modn<-gam(totno~as.factor(year)+s(day2,k=4,bs='cc') + s(time,bs='cc',k=4),data=d,gamma=1.4,family='nb'(link='log'))
    pdat0<-data.frame(year=sort(unique(d$year),decreasing=TRUE)[1],
                      time=1200,
                     day2=seq(min(d$day2),max(d$day2),1))
    p<-predict(modn,newdata=pdat0,type='response')
    pdat0$p<-p
pdat<-data.frame(day2=seq(1,365,1))
pdat$p<-pasp<-aspline(pdat0$day2,pdat0$p,xout=pdat$day2)$y
    pdat<-revdayshiftfun(pdat,shf)
    pdatn<-subset(pdat,select=c('day3','p'))
    names(pdatn)<-c('day','pno')

pdat<-merge(pdatw,pdatn,by=c('day'),all=FALSE)
pdat$sz<-pdat$pwt/pdat$pno
pdat$szz<-(pdat$sz-mean(pdat$sz))/sd(pdat$sz)
pdat$lon<-unique(d$lonc)
pdat$lat<-unique(d$latc)
pdat$did<-ifelse(pdat$day %in% d$day,1,0)
return(data.frame(pdat))
}
#LENGTH
out<-ddply(subset(rvw3, cell %in% dmm$cell),.(cell),.fun=ff,.progress='text')
xyplot(log10(sz+1)~day|cell,data=out,col=ifelse(out$did==1,'red','blue'))
xyplot(sz~day|cell,data=out,col=ifelse(out$did==1,'red','blue'))
xyplot(szz~day|cell,data=out,col=ifelse(out$did==1,'red','blue'))


bnk<-subset(plg,stratum %in% c(443,458,455,456,464,463,473,474,475,480))
x1<--68
x2<--57
y<-41.5
out$dy<-rescale(out$day,newrange=c(x1,x2))
frames <- length(unique(out$day))

rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',i,'plot.jpg',sep=''))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste('00',i,'plot.jpg', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', i,'plot.jpg', sep=''))
  }
}


out$pcx<-rescale(log10(out$sz+1),newrange=c(.2,7))
setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/Figures/prac')
#loop through plots
for(i in 1:frames){
  name <- rename(i)
  d<-subset(out,day==i)

#saves the plot as a file in the working directory
jpeg(name,width=6,height=5,units='in',quality=100,res=300)
map('worldHires',fill=TRUE,col='gray',border=NA,xlim=c(x1,x2),ylim=c(y,46))
#plot(plg,add=TRUE,lwd=.01,border=alpha('lightgray',.5))
plot(bnk,add=TRUE,lwd=.01,border=NA,col=alpha('firebrick3',.2))
points(d$lon,d$lat,pch=16,cex=d$pcx,col=alpha('dodgerblue3',.3),xlim=c(x1,x2),ylim=c(y,46))
points(d$lon,d$lat,pch=1,lwd=.01,cex=d$pcx,col='dodgerblue3',xlim=c(-x1,x2),ylim=c(y,46))
legend('bottomright',paste('Day=',i),bty='n')
points(unique(d$dy),y,pch=17,col='red3',cex=2)
axis(1,at=seq(x1,x2,length.out=14),labels=c('','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC',''),cex.axis=.7)
dev.off()
}


out$pcx<-rescale(out$szz,newrange=c(.2,7))
setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/Figures/prac')
#loop through plots
for(i in 1:frames){
  name <- rename(i)
  d<-subset(out,day==i)

#saves the plot as a file in the working directory
jpeg(name,width=6,height=5,units='in',quality=100,res=300)
map('worldHires',fill=TRUE,col='gray',border=NA,xlim=c(x1,x2),ylim=c(y,46))
#plot(plg,add=TRUE,lwd=.01,border=alpha('lightgray',.5))
plot(bnk,add=TRUE,lwd=.01,border=NA,col=alpha('firebrick3',.2))
points(d$lon,d$lat,pch=16,cex=d$pcx,col=alpha('dodgerblue3',.3),xlim=c(x1,x2),ylim=c(y,46))
points(d$lon,d$lat,pch=1,lwd=.01,cex=d$pcx,col='dodgerblue3',xlim=c(-x1,x2),ylim=c(y,46))
legend('bottomright',paste('Day=',i),bty='n')
points(unique(d$dy),y,pch=17,col='red3',cex=2)
axis(1,at=seq(x1,x2,length.out=14),labels=c('','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC',''),cex.axis=.7)
dev.off()
}









a<-subset(rvw2,geardesc=='Campelen 1800 survey trawl')#2002 and 2005 during october
plot(a$lon,a$lat,pch=16,col='red')
setwd(figsdir)
pdf('campelen.survey.scotianshelf.pdf')
map('world',col='gray',fill=TRUE,xlim=c(-67,-58),ylim=c(43,45))
map.axes()
points(a$lon,a$lat,pch=16,col=alpha('red',.5))
dev.off()






#####################################################################

#####################################################################
#SMOOTH TREND IN PHENOLOGY BY 1/4 CELL FOR CLUSTER (LONG FORM)
rvw2<-subset(rvw,!(geardesc %in% c('Newston net','Campelen 1800 survey trawl')))

f<-function(d){
if(dim(subset(d,month %in% c(1,2,3)))[1]>10){wint<-1
                                         } else { wint<-0 }
if(dim(subset(d,month %in% c(6,7,8)))[1]>10){summ<-1
                                         } else { summ<-0 }
if(dim(subset(d,month %in% c(9,10,11)))[1]>10){fall<-1
                                         } else { fall<-0 }
    return(data.frame(nmo=length(unique(d$month)),
                      ndy=length(unique(d$day)),
                      nyr=length(unique(d$year)),
                      wint=wint,
                      summ=summ,
                      fall=fall,
                      tot=wint+summ+fall,
                      lon=mean(d$lon),
                      lat=mean(d$lat)))
}
dt<-ddply(rvw2,.(cell),.fun=f)
dt<-dt[order(dt$nmo,decreasing=TRUE),]
dmm<-subset(dt,nmo>= 4 & ndy>=30 & nyr>=3 & tot>=3)


f1<-function(d0){
f<-function(d){
day<-unique(d$day)
d2<-subset(d0,select=c('day'))
d2$day2<-d2$day-(day-1)
d2$day2<-ifelse(d2$day2<=0,d2$day2+(365),d2$day2)
    return(data.frame(day=day,
                   span=diff(range(min(d2$day2),max(d2$day2)))))
}
z<-ddply(d0,.(day),.fun=f)
z<-subset(z,span==max(z$span))
return(subset(z,day==min(z$day)))
}
dayshift<-ddply(unique(subset(rvw2,cell %in% dmm$cell,select=c('strat','day','cell'))),.(cell),.fun=f1,.progress='text')

##############################################
#ADJUSTS DAY VALUES TO MAXIMIZE PHENOLOGY SPAN
dayshiftfun<-function(d,shf){
d$day2<-d$day-(shf)
d$day2<-ifelse(d$day2<=0,d$day2+(365),d$day2)
return(d)
}

#BACK-CALCULATE ORIGINAL DAY OF THE YEAR
revdayshiftfun<-function(d,shf){
d$day3<-d$day2+shf
d$day3<-ifelse(d$day3>365,d$day3-365,d$day3)
return(d)
}

d<-subset(rvw2,cell=="-57.75_44.25",select=c('totwgt','cell','year','time','lon','lat','day'))


fsm<-function(d){
shf<-subset(dayshift,cell==unique(d$cell))$day
names(d)[1]<-'y'
#SHIFT DAYS TO CENTER ON JULY AND MAKE CONTINUOUS
d<-dayshiftfun(d,shf)
mod<-gam(y~s(year) + s(time,k=4,bs='cc') + s(day2,k=4,bs='cc'),data=d,gamma=1.4,family='nb'(link='log'))
    pdat0<-data.frame(year=2000,
                     time=1200,
                     day2=seq(min(d$day2),max(d$day2),1))
    p<-predict(mod,newdata=pdat0,type='response')
    pdat0$p<-p
    pdat<-data.frame(day2=seq(1,365,1))
    pdat$p<-aspline(pdat0$day2,pdat0$p,xout=pdat$day2)$y
    pdat$pz<-(pdat$p-mean(pdat$p))/sd(pdat$p)
    pdat<-revdayshiftfun(pdat,shf)
    pdat<-subset(pdat,select=c('day3','pz'))
names(pdat)<-c('day',gsub(' ','',paste('X_',unique(d$cell))))
pdat<-pdat[order(pdat$day),]
return(data.frame(pdat))
}
#TOTAL WEIGHT
qq<-dlply(subset(rvw2,cell %in% dmm$cell,select=c('totwgt','cell','year','time','lon','lat','day')),.(cell),.fun=fsm)
qsm<-Reduce(function(x, y) merge(x, y, by=c('day'),all=TRUE), qq)
qsm<-qsm[,colSums(is.na(qsm)) != nrow(qsm)]


q<-qsm
q2<- q %>% gather(cell, value, -day)
xyplot(value ~ day | cell,data=q2, pch=15, type=c('spline'),col='black')
xyplot(value ~ year | strata,data=q2, pch=15, type=c('p','spline'),col='black')




clfun<-function(df,k,lbl){
rownames(df)<-df$year
df<-df[,-1]
#k<-3#ER OF CLUSTERS
dmat<-1-cor(df,use='pairwise.complete.obs')
dst<-as.dist(dmat)
ff<-fanny(dst,k,maxit=5000,diss=T)

par(mfrow=c(3,1),mar=c(4,12,1,12))
dum<-c('red3','darkblue','gold3')
plot(silhouette(ff),col=dum[1:k],main='')#silhouette plot

dc.pcoa<-cmdscale(dst)
dc.scores<-scores(dc.pcoa,choices=c(1,2))
spefuz.g<-ff$clustering
a<-data.frame(cell=as.character(sort(unique(names(df)))),
              clusters=ff$clustering)
aa<-data.frame(ff$membership)
aa$cell<-rownames(a)

plot(scores(dc.pcoa),asp=1,type='n',xlim=c(-1,1.5),ylim=c(-1,1),las=1,axes=TRUE,xlab='',ylab='')
stars(ff$membership,location=scores(dc.pcoa),draw.segments=T,add=T,scale=F,len=.1,col.segments=alpha(c(dum[1:k]),.25),byt='n',labels=NULL,xlim=c(-1.1,1.4),ylim=c(-1,.5),lwd=.0001,xpd=TRUE,border=NULL,radius=FALSE,col.radius=alpha('white',.1))
for(i in 1:k){ cl<-dum[i]
    gg<-dc.scores[spefuz.g==i,]
    hpts<-chull(gg)
    hpts<-c(hpts,hpts[1])
    lines(gg[hpts,],col=cl,lwd=3,xlim=c(-1.1,1.4),ylim=c(-1,.5))
}

cx <- data.frame(cell=aa$cell,
                 cx=apply(aa[, 1:k], 1, max))
cx$cxx <- rescale(cx$cx,newrange=c(.05,.5))
#PLOTS THE TIMESERIES OF R FOR EACH CLUSTER
par(mfrow=c(3,1),mar=c(2,12,1,12),oma=c(1,1,1,1))
#par(mfrow=c(2,2),mar=c(2,12,1,12),oma=c(1,1,1,1))

mb<-seq(1,k,1)
l<-list()
for(i in 1:length(mb)){
print(mb[i])
one<-subset(a,clusters==mb[i],select=c('cell'))# IN THIS CLUSTER
cx2<-subset(cx,cx>=0.5)#GETS INSTANCES WHERE CLUSTER PROBABILITY>0.5
data5<-subset(df,select=c(as.character(one$cell)))#TS FOR CLUSTER
cx2<-subset(cx2,cell %in% names(data5))
data5<-subset(df,select=c(as.character(cx2$cell)))#TS FOR CLUSTER
x<-rownames(data5)
t<-data.frame(day=as.numeric(x),mn=rowMeans(data5,na.rm=F))
t2<-data.frame(day=as.numeric(x),mn=rowMeans(data5,na.rm=T))

cl<-dum[i]
plot(0,0,pch=16,cex=.01,xlim=c(1,365),ylim=c(-2,3),main='',xaxt='n',las=1,axes=FALSE,xlab='',ylab='')
axis(side=2,at=seq(-2,3,1),las=1,lwd=.001,cex.axis=.75)
axis(side=1,at=seq(1,365,10),cex.axis=.75)
    for(j in 1:length(data5[1,])){
        try(dat<-data5[,j])
        try(trnsp<-subset(cx,cell==as.character(names(data5[j])))$cxx)
        try(lines(x,dat,cex=.5,ylim=c(0,1),lwd=2,col=alpha(cl,rescale(trnsp,newrange=c(.1,.5)))))
    }
lines(as.numeric(as.character(t$day)),t$mn,col='gold3',cex=.7,lwd=3)
dm<-subset(t2,mn==max(t2$mn,na.rm=TRUE),select=c('day'))
dm$cluster=mb[i]
names(dm)<-c('clusterday','clusters')
l[[i]]<-dm
}


a2<-a
a2$cell<-gsub('X_\\.','',a2$cell)
if(k==3){a2$cx<- apply(aa[,c('X1','X2','X3')], 1, function(x) max(x) )
     } else {a2$cx<- apply(aa[,c('X1','X2')], 1, function(x) max(x) )
     }
a2$id<-a2$cell
a2$cl<-ifelse(a2$clusters==1,dum[1],dum[3])
a2$cl<-ifelse(a2$clusters==2,dum[2],a2$cl)
crds<-unique(subset(rvw,select=c('cell','lonc','latc')))
crds$cell<-gsub('-','',crds$cell)
a2<-merge(a2,crds,by=c('cell'),all.x=TRUE,all.y=FALSE)
dum<-unique(subset(a2,select=c('clusters','cl')))

return(
    ggplot()+
geom_polygon(aes(long,lat, group=group),fill='black', data=bnk,size=1)+
geom_tile(data=a2, aes(x=lonc, y=latc,fill=as.factor(clusters),alpha=cx),col='gray80',size=.0001)+
scale_fill_manual(breaks=as.character(dum$clusters),values=dum$cl,na.value="transparent",guide=guide_legend(title=''))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.9,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.1, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-68,-57,1),labels=as.character(seq(-68,-57,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(41,47,1),labels=as.character(seq(41,47,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,47),xlim=c(-68,-57))+
    xlab('')+
    ylab('')+
         labs(title = lbl,
       subtitle = "",
       caption = '')
   )
}


setwd(figsdir)
pdf('herring_rvsurvey_cluster_phenology_grid_k3.pdf',height=10,width=6.5)
mp1<-clfun(qsm,3,"Average weight")
dev.off()

pdf('herring_rvsurvey_cluster_phenology_map_grid_k3.pdf',height=14,width=7)
grid.arrange(mp1,mp1,mp1,ncol=1)
dev.off()




#MAKES BARPLOT OF CLUSTER DISTRIBUTION BY LONGITUDE
setwd(figsdir)
pdf('herring_rvsurvey_cluster_phenology_extras.pdf',height=14,width=7)
f<-function(d){
f2<-function(d2){    return(data.frame(prp=dim(d2)[1]/dim(d)[1]))}
return(ddply(d,.(clusters),.fun=f2))
}
dd<-ddply(a2,.(lonc),.fun=f)

p2<-ggplot(dd, aes(x=lonc, fill=as.factor(clusters),y=prp))+
    geom_bar(data=dd,stat='identity',aes(width=.7,order=clusters),size=.0001,col='gray')+
    theme(legend.position=c(.8,.8),text = element_text(size=8),axis.text.x = element_text(angle=0, vjust=1), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),panel.background = element_blank(), legend.key.size =  unit(0.2, "in"),legend.text=element_text(size=6),axis.line = element_line(color="black", size = .1))+
    scale_fill_manual(values=as.character(dum$cl),breaks=as.factor(dum$clusters),labels=as.character(dum$clusters),name='Cluster')+
    scale_x_continuous(expand=c(0,0),breaks=seq(min(dd$lonc),max(dd$lonc),10),labels=round(seq(min(dd$lonc),max(dd$lonc),10),digits=0))+
    expand_limits(x=c(min(dd$lonc),max(dd$lonc)))+
    xlab('Longitude') +
    ylab('Proportion of cells')


#MODEL PROBABILITY OF CLUSTER ASSIGNMENT BY LONGITUDE
f<-function(d,cl){
d$ps<-ifelse(d$clusters==cl,1,0)
mod<-gam(ps~s(lonc,k=6),data=d,family='binomial',weights=d$cx,gamma=1.5)
pdat<-data.frame(lonc=seq(min(d$lonc),max(d$lonc),length.out=100))
p<-predict(mod,newdata=pdat,se.fit=TRUE,type='response')
pdat$p<-p$fit
return(pdat)
}
c1<-f(a2,1)
c2<-f(a2,2)
c3<-f(a2,3)
par(mfrow=c(3,1),mar=c(4,4,1,1))
plot(c1$lonc,c1$p,ylim=c(0,1),col=dum$cl[1],las=1,pch=16,type='l',lwd=2,ylab='Proportion of all cells')
points(c2$lonc,c2$p,col=dum$cl[2],pch=16,type='l',lwd=2)
points(c3$lonc,c3$p,col=dum$cl[3],pch=16,type='l',lwd=2)

cx<-.3
plot(0,0,col='white',ylim=c(0,1),las=1,pch=16,type='l',lwd=2,ylab='Proportion of all cells',xlim=c(min(c1$lonc),max(c1$lonc)))
polygon(c(c2$lonc,c2$lonc[length(c2$lonc):1]),c(c2$p,rep(-5,dim(c2)[1])[length(c2$p):1]),col=alpha(dum$cl[2],cx),border=alpha(dum$cl[2],.8))
polygon(c(c1$lonc,c1$lonc[length(c1$lonc):1]),c(c1$p,rep(-5,dim(c1)[1])[length(c1$p):1]),col=alpha(dum$cl[1],cx),border=alpha(dum$cl[1],.8))
polygon(c(c3$lonc,c3$lonc[length(c3$lonc):1]),c(c3$p,rep(-5,dim(c3)[1])[length(c3$p):1]),col=alpha(dum$cl[3],cx),border=alpha(dum$cl[3],.8))

grid.arrange(p2,p2,p2,ncol=1)
dev.off()








#devtools::install_github('dgrtwo/gganimate',force=TRUE)
library(magick)
library(gganimate)
plot(coast.mc)
p<-bmap+
    geom_point(aes(x=lon, y=lat,size=pz,frame=day),data=out,colour='purple',alpha=.5)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.1,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'))

out$strat<-as.character(out$strat)
o<-data.frame(strat=sort(unique(out$strat)),
n=tapply(out$day,out$strat,function(x) length(unique(x))))

library(ggthemes)
library(animation)
setwd(figsdir)
ani.options(interval = 0.2)
  gganimate(p)
  gganimate(p,'output.phen.gif')
  gganimate(p,'output.phen.mp4')
  gganimate(p,'output.smooth.mp4')


p



a<-subset(rvw2,geardesc=='Campelen 1800 survey trawl')#2002 and 2005 during october

plot(a$lon,a$lat,pch=16,col='red')
setwd(figsdir)
pdf('campelen.survey.scotianshelf.pdf')
map('world',col='gray',fill=TRUE,xlim=c(-67,-58),ylim=c(43,45))
map.axes()
points(a$lon,a$lat,pch=16,col=alpha('red',.5))
dev.off()









####################################################################

####################################################################
#ANIMATE: ESTIMATE TIME TREND IN AVERAGE BIOMASS AND ANIMATE
###########################################################
#ESTIMATE PHENOLOGY FOR EACH GRID CELL AND ANIMATE
rvw2<-subset(rvw,!(geardesc %in% c('Newston net','Campelen 1800 survey trawl')) & month %in% c(6,7,8))
f<-function(d){
    return(data.frame(nmo=length(unique(d$month)),
                      ndy=length(unique(d$day)),
                      nyr=length(unique(d$year)),
                      myr=min(d$year),
                      nlon=length(unique(d$lon)),
                      nlat=length(unique(d$lat)),
                      ntm=length(unique(d$time)),
                      lon=unique(d$lonc),
                      lat=unique(d$latc)))
}
dt<-ddply(rvw2,.(cell),.fun=f)
dt<-dt[order(dt$nmo,decreasing=TRUE),]
dmm<-subset(dt,nyr>=10 & myr<1990)



library(akima)
ff<-function(d){
names(d)[1]<-'y'
#SHIFT DAYS TO CENTER ON JULY AND MAKE CONTINUOUS
mod<-gam(y~s(year,k=6) + s(time,k=4,bs='cc'),data=d,gamma=1.4,family='nb'(link='log'))
    pdat<-data.frame(year=seq(min(d$year),max(d$year),.25),
                      time=1200,
                      lon=unique(d$lonc),
                      lat=unique(d$latc))
    p<-predict(mod,newdata=pdat,type='response')
    pdat$p<-p
    pdat$pz<-(pdat$p-mean(pdat$p))/sd(pdat$p)
    pdat<-subset(pdat,select=c('year','p','pz'))
pdat$lon<-unique(d$lonc)
pdat$lat<-unique(d$latc)
pdat$depth<-mean(d$dmax,na.rm=TRUE)
pdat$did<-ifelse(pdat$year %in% d$year,1,0)
return(data.frame(pdat))
}
#TOTAL WEIGHT
out<-ddply(subset(rvw2, cell %in% dmm$cell,select=c('totwgt','strat','year','day','time','lon','lat','cell','lonc','latc','dmax')),.(cell),.fun=ff,.progress='text')
xyplot(pz~year|cell,data=out,col=ifelse(out$did==1,'red','blue'))
xyplot(log10(p+1)~year|cell,data=out,col=ifelse(out$did==1,'red','blue'))


bnk<-subset(plg,stratum %in% c(443,458,455,456,464,463,473,474,475,480))
x1<--68
x2<--57
y<-41.5
out$dy<-rescale(out$year,newrange=c(x1,x2))
frames <- length(unique(out$year))

rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',i,'plot.jpg',sep=''))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste('00',i,'plot.jpg', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', i,'plot.jpg', sep=''))
  }
}


out$pcx<-rescale(log10(out$p+1),newrange=c(.2,7))
setwd('C:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/Figures/prac')
#loop through plots
yr<-sort(unique(out$year))
for(i in 1:frames){
  name <- rename(i)
  d<-subset(out,year==yr[i])

#saves the plot as a file in the working directory
jpeg(name,width=6,height=5,units='in',quality=100,res=300)
map('worldHires',fill=TRUE,col='gray',border=NA,xlim=c(x1,x2),ylim=c(y,46))
#plot(plg,add=TRUE,lwd=.01,border=alpha('lightgray',.5))
plot(bnk,add=TRUE,lwd=.01,border=NA,col=alpha('firebrick3',.2))
points(d$lon,d$lat,pch=16,cex=d$pcx,col=alpha('dodgerblue3',.3),xlim=c(x1,x2),ylim=c(y,46))
points(d$lon,d$lat,pch=1,lwd=.01,cex=d$pcx,col='dodgerblue3',xlim=c(-x1,x2),ylim=c(y,46))
legend('bottomright',paste('Year=',round(yr[i],digits=0)),bty='n')
points(unique(d$dy),y,pch=17,col='red3',cex=2)
axis(1,at=seq(x1,x2,length.out=10),labels=seq(1970,2015,5),cex.axis=.7)
dev.off()
}


#FROM COMMAND LINE WITH IMAGEMAGIK INSTALLED RUN: CONVERT -DELAY 2 -QUALITY 100 *.JPG MOVIE.MP4
#run ImageMagick
#library(magick)
#my_command <- 'convert *.png -delay 2 -loop 0 animation.gif'
#system(my_command)


