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


#datadir1<-'/scratch/dboyce/spera/data/stagingdat'
#datadir<-'/scratch/dboyce/spera/data/finaldat'
#figsdir<-'/scratch/dboyce/spera/figs'

datadir1<-'N://cluster_2017//scratch//spera//data//stagingdat'
datadir<-'N://cluster_2017//scratch//spera//data//finaldat_v2'
figsdir<-'C://Users//copepod//Documents//aalldocuments//literature//research//active//SPERA//Figures'
#datadir1<-'/scratch/dboyce/spera/data/stagingdat'
#datadir<-'/scratch/dboyce/spera/data/finaldat_v2'
#figsdir<-'/scratch/dboyce/spera/figs_v2'



setwd('N://cluster_2017//scratch//spera//data//finaldat')
setwd(datadir)
load("pathfinder_sst_daily_1981_2012_spawnar.RData")
a<-unique(subset(dat,select=c('lon','lat')))
map('world',xlim=c(-70,-65),ylim=c(42,46))
points(a$lon,a$lat,pch=16,col='red')
dm<-data.frame(year=sort(unique(dat$year)),
               sst=tapply(dat$sst,dat$year,mean))
plot(dm$year,dm$sst,pch=15)
dat$cell<-gsub(' ','',paste(dat$lon,'_',dat$lat))

d<-subset(dat,cell==sort(unique(dat$cell))[1])

#ESTIMATE AVERAGE SST DURING FALL FOR EACH YEAR
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
pdat0<-subset(pdat,p>=14)
dur20<-dim(pdat0)[1]
tim20<-min(pdat0$day)

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
                    dur20=dur20,
                    tim20=tim20,
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

#######################################
#ABUNDANCE AT AGE FROM LANDINGS
hera<-read.csv('numbersAtAge_1965-2016_SWNS-BoF_stock.csv',header=TRUE,skip=1)
names(hera)[1]<-'year'

#CHANGE DATA FRAME FORMAT
l<-list()
for(i in 2:dim(hera)[2]){
print(i)
d<-data.frame(cbind(hera$year,hera[,i]))
d$age<-i-1
names(d)[1:2]<-c('year','catch')
l[[i-1]]<-d
}
hera<-data.frame(do.call('rbind',l))

#LINEAR TIME TREND IN LANDINGS OF ALL AGES
f<-function(d){
print(unique(d$age))
d$catch<-(d$catch-mean(d$catch))/sd(d$catch)
mod<-lm(catch~year,data=d)
    s<-summary(mod)
return(data.frame(age=unique(d$age),
                  beta=s$coef[2,1]))
}
tdat<-ddply(hera,.(age),.fun=f)
plot(tdat$age,tdat$beta)

#TREND IN PROPORTION OF OLD HERRING IN CATCH OVER TIME
f<-function(d){
    a<-subset(d,age>=5)
    return(data.frame(prp=(sum(a$catch)/sum(d$catch))*100))
}
tdat2<-ddply(hera,.(year),.fun=f)
plot(tdat2$year,tdat2$prp,pch=15,type='b')

#EVENNESS OF AGES IN CATCH OVER TIME
f<-function(d){
d2<-subset(d,catch>0)
tcatch<-sum(d$catch)

#SHANNON I
f2<-function(d3){ return(data.frame(p=sum(d3$catch)/tcatch)) }
sdat<-ddply(d,.(age),.fun=f2)
sdat$lp<-log(sdat$p)
sdat$pp<--1*(sdat$p*sdat$lp)
sdat<-subset(sdat,is.na(pp)==FALSE)

dout<-data.frame(shannon=sum(sdat$pp,na.rm=TRUE),
                 pe=sum(sdat$pp)/(log(length(unique(d2$age)))))
return(dout)
}
tdat3<-ddply(hera,.(year),.fun=f)
plot(tdat3$year,tdat3$pe,pch=15,type='b')

#ABUNDANCE WEIGHTED AGE OVER TIME
mod<-lm(age~as.factor(year),weights=hera$catch,data=hera)
pdat<-data.frame(year=sort(unique(hera$year)))
p<-predict(mod,newdata=data.frame(year=pdat$year))
pdat$p<-p
plot(pdat$year,pdat$p,pch=15,type='b')
par(new=TRUE)
plot(her$year,her$ssb,pch=15,col='red',type='b',yaxt='n')

#PLOTS VPA AND ACOUSTIC SSB ESTIMATES; CALIBRATES ACOUSTIC
plot(her$year,her$ssb,pch=16,xlim=c(1965,2016),type='l')
points(her$year,her$acoustic,pch=16,col='red',type='l')
points(her$year,her$ssb,pch=16,xlim=c(1965,2016),col=alpha('black',.3))
points(her$year,her$acoustic,pch=16,col=alpha('red',.3))
mod<-lm(ssb~acoustic,data=her)
her$acoustic2<-predict(mod,newdata=data.frame(acoustic=her$acoustic))
plot(her$year,her$ssb,pch=16,xlim=c(1965,2016),type='l',xlab='Year',ylab='SSB')
points(her$year,her$acoustic2,pch=16,col='red',type='l')
points(her$year,her$ssb,pch=16,xlim=c(1965,2016),col=alpha('black',.3))
points(her$year,her$acoustic2,pch=16,col=alpha('red',.3))

mod<-lm(acoustic~ssb,data=her)
her$ssb2<-predict(mod,newdata=data.frame(ssb=her$ssb))
plot(her$year,her$ssb2,pch=16,xlim=c(1965,2016),type='l',las=1,xlab='Year',ylab='SSB')
points(her$year,her$acoustic,pch=16,col='red',type='l')
points(her$year,her$ssb2,pch=16,xlim=c(1965,2016),col=alpha('black',.3))
points(her$year,her$acoustic,pch=16,col=alpha('red',.3))

cor(her$ssb,her$acoustic,use='pairwise.complete.obs')
cor(her$ssb,her$acoustic,use='pairwise.complete.obs',method='kendall')



#HERRING DATA FROM ASSESSMENTS: BIOMASS, WEIGHT, F, ETC...
setwd(datadir1)
her<-read.csv('herring_assess_ssb_r_f_w_spera.csv',header=TRUE)
#AVERAGE F PER YEAR
f<-function(d){
    return(data.frame(f=mean(d$fage2,d$fage3,d$fage4,d$fage5,d$fage6,d$fage7,d$fage8,d$fage9,d$fage10,d$fage11,na.rm=TRUE)))
}
d2<-ddply(her,.(year),.fun=f)
her<-merge(her,d2,by=c('year'),all.x=TRUE,all.y=FALSE)

#AVERAGE WEIGHT PER YEAR
f<-function(d){
 return(data.frame(wgt=mean(d$wage7,d$wage8,d$wage9,d$wage10,d$wage11,na.rm=TRUE)))
}
d2<-ddply(her,.(year),.fun=f)
her<-merge(her,d2,by=c('year'),all.x=TRUE,all.y=FALSE)


setwd(datadir)
nao<-read.csv('nao.csv',header=TRUE)
load('phytoplankton_biomass_all_spera_spawnar.RData');phyt<-dat
path<-load('pathfinder_sst_daily_1981_2012_spawnar.RData');path<-dat
nit<-read.csv("Atlas1999nuts_nit_spawnar.csv",header=TRUE)
phos<-read.csv("Atlas1999nuts_phos_spawnar.csv",header=TRUE)
sil<-read.csv("Atlas1999nuts_sil_spawnar.csv",header=TRUE)
sz<-read.csv("phyto_size_monthly_1997_2007_spera_spawnar.csv",header=TRUE)
fg<-read.csv("phyto_functional_monthly_1997_2010_spera_spawnar.csv",header=TRUE)
strat<-read.csv("physical_stratification_spera_spawnar.csv",header=TRUE)
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
zto<-read.csv("zoop_turnover_spera_spawnar.csv",header=TRUE)
pto<-read.csv("phyto_turnover_spera_spawnar.csv",header=TRUE)
plank<-read.csv('plank_cpr_richness_spera_spawnar.csv',header=TRUE)
vtemp<-read.csv('vtemp.csv',header=TRUE)
vwind<-read.csv('vwind.csv',header=TRUE)
#ctdat<-read.csv('ctdat.csv',header=TRUE)
ctdat<-read.csv('ct_data.csv',header=TRUE)
names(ctdat)[1]<-'year'
load("BoF_larval_herring_lengths_spera_spawnar.RData")
larvs<-dat

###########################################################
#EXAMINE THERMAL TOLERANCE RANGE OF HERRING FROM RV SURVERY
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw2<-subset(rvw,is.na(surface_temperature)==FALSE)#ONLY 358 MISSING
rvw2$temp<-round((rvw2$bottom_temperature+rvw2$surface_temperature)/2,digits=1)
rvw2$btemp<-round(rvw2$bottom_temperature,digits=0)
rvw2$stemp<-round(rvw2$surface_temperature,digits=0)
rvw2$occ<-ifelse(rvw2$totno>0,TRUE,FALSE)
save(rvw2,file='rvw2.RData')

stdat<-data.frame(temp=sort(unique(rvw2$temp)),
stemp=tapply(rvw2$surface_temperature,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
btemp=tapply(rvw2$bottom_temperature,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
swgt=tapply(rvw2$sampwgt,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
twgt=tapply(rvw2$totwgt,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
tno=tapply(rvw2$totno,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
n=tapply(rvw2$totno,rvw2$temp,function(x) length(x)))
save(stdat,file='stdat.RData')




#CORRECTION AS GIVEN BY MIKE MCMAHON: 2016 SWITCH FROM CM TO MM AND FROM FORK LENGTH TO TOTAL LENGTH
rvl<-read.csv("herring_lengths_RV_survey_spera_spawnar.csv",header=TRUE)
rvl$flen<-ifelse(rvl$mission=='NED2016016',rvl$flen/10,rvl$flen)
rvl$flen<-ifelse(rvl$mission!='NED2016016',(rvl$flen*1.0866)+0.95632,rvl$flen)
plot(log10(rvl$flen),log10(rvl$fwt))
rvl<-subset(rvl,log10(rvl$flen)>=.5)#REMOVE OUTLIERS

#HERRING BIOMASS FROM SUMMER RV
setwd(datadir1)
rvww<-read.csv('herring_weights_RV_survey_spera.csv',header=TRUE)
names(rvww)<-tolower(names(rvww))

crds<-subset(rvww,select=c('lon','lat'))
map('world',xlim=c(-70,-64),ylim=c(42,46))
points(crds$lon,crds$lat,pch=16,col='red')

#POSITION OF SLOPE WATERS DURING FALL
setwd(datadir)
cur<-read.csv('GS_SS_distances.csv',header=TRUE)
f<-function(d){
    return(data.frame(GS.Dist=mean(d$GS.Dist),
               SS.Dist=mean(d$SS.Dist),
               GS.Dist.sd=sd(d$GS.Dist),
               SS.Dist.sd=sd(d$SS.Dist)))
}
cur<-ddply(subset(cur,month>=8 & month<12),.(year),.fun=f,.progress='text')


#HERRING LARVAE
#larv<-read.csv("BoF_larval_herring_counts_spera_spawnar.csv",header=TRUE)
load("BoF_larval_herring_counts_spera_spawnar.RData")
larv<-subset(dat,bathy<0)
#STANDARDIZE PER VOLUME STRAINED
larv$stdno<-(larv$totno/larv$volm3)*abs(larv$bathy)#NUMBER PER M2 - INTEGRATED OVER DEPTH
#larv$stdno<-(larv$totno/larv$volm3)
larv<-subset(larv,is.na(stdno)==FALSE)
larv<-subset(larv,gearprocdesc=='Sawtooth Oblique')
#SPAWNING ON GERMAN BANK AUG-SEPT, EARLIER ON LURCHER
larv<-subset(larv,month>=8)
larv$phylum<-tolower(larv$phylum)
larv$order<-tolower(larv$order)
larv$kingdom<-tolower(larv$kingdom)

#PREDICT ANNUAL AVERAGES FOR HERRING LARVAE LENGTH
load("BoF_larval_herring_lengths_spera_spawnar.RData");larvs<-dat
modf<-gam(length~as.factor(year) + s(lon,lat,k=20) + s(day,bs='cc',k=5),weights=larvs$clen,data=larvs,gamma=1.4)
hlen<-data.frame(year=sort(unique(larvs$year)),
                 lon=median(larvs$lon),
                 lat=median(larvs$lat),
                 day=230)
p<-predict(modf,newdata=pdat,se.fit=TRUE)
hlen$p<-p$fit
hlen$se<-p$se.fit

plot(hlen$year,hlen$p,pch=15,xlim=c(1970,2010),las=1,ylab='Herring larvae length',xlab='Year')
par(new=TRUE)
plot(her$year,her$larv,pch=16,col='red',ylim=c(0,100),xlim=c(1970,2010),axes=FALSE,xlab='',ylab='')
axis(4,seq(0,100,25),labels=TRUE,las=1)
a<-merge(hlen,her,by=c('year'),all=FALSE)
a<-na.omit(subset(a,select=c('year','p','larv')))
r<-cor(log(a$p),log(a$larv),use='pairwise.complete.obs')
par(fig = c(0,1,0,1))
plot(log(a$p),log(a$larv),pch=16,las=1,xlab='Herring larvae length',ylab='Herring larvae abundance')
legend('topright',legend=paste('r =',round(r,digits=2)),bty='n')

a<-merge(hlen,her,by=c('year'),all=FALSE)
a<-na.omit(subset(a,select=c('p','f')))
cor(a$p,a$f,use='pairwise.complete.obs')
ccf(a$p,a$f)
ccf(log(a$p),log(a$f))
plot(a$p,a$f)
plot(log(a$p),log(a$f))



#SUM OF COUNTS (NUMBER OF INDIVIDUALS) PER SAMPLE FOR ZOOPLANKTON
avfun<-function(d){
dout<-unique(subset(d,select=c('tday','tbin3','tbin5','tbin10','surftemp','bottemp','airtemp')))
dout$stdno<-sum(d$stdno)

dout$var<-'T zoop'
return(dout)
}
zp.juv<-ddply(subset(larv,order %in% c('calanoida','cyclopoida','decapoda') | phylum %in% c('mollusca')),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')

zp.ad<-ddply(subset(larv,order %in% c('euphausiacea','calanoida','cyclopoida') | phylum %in% c('chaetognatha')),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')

zp<-ddply(larv,.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')




#CALCULATE SUM OF COUNTS (NUMBER OF INDIVIDUALS) PER SAMPLE
avfun<-function(d){
 dout<-unique(subset(d,select=c('tday','tbin3','tbin5','tbin10','surftemp','bottemp','airtemp')))
 dout$stdno<-sum(d$stdno)

cname<-unique(na.omit(d$comname))
pname<-unique(na.omit(d$phylum))
oname<-unique(na.omit(d$order))
kname<-unique(na.omit(d$kingdom))

 if(length(cname)==1){
     if(cname=='herring(atlantic)'){ dout$var<-'Herring'}
 } else NULL

if(length(pname)==1){
if(pname=='arthropoda'){dout$var<-'Arthropod'}
if(pname=='foraminifera'){dout$var<-'Foraminifera'}
if(pname=='chaetognatha'){dout$var<-'Chaetognath'}
if(pname=='annelida'){dout$var<-'Annelids'}
if(pname=='mollusca'){dout$var<-'Mollusc'}
if(pname=='cnidaria'){dout$var<-'Cnidarian'}
     if(pname=='ctenophora'){dout$var<-'Ctenophore'}
 } else NULL

if(length(oname)==1){
if(oname=='gadiformes'){dout$var<-'Gadids'}
if(oname=='stomiiformes'){dout$var<-'Dragonfish'}
if(oname=='myctophiformes'){dout$var<-'Myctophids'}
if(oname=='pleuronectiformes'){dout$var<-'Flatfish'}
 } else NULL

if(length(pname)>1){
if(kname=='animalia'){dout$var<-'T L Fish'}
 } else NULL

 return(dout)
}
l<-list()
l[[1]]<-ddply(subset(larv,comname=='herring(atlantic)'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[2]]<-ddply(subset(larv,phylum=='arthropoda'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[3]]<-ddply(subset(larv,phylum=='foraminifera'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[4]]<-ddply(subset(larv,phylum=='chaetognatha'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[5]]<-ddply(subset(larv,phylum=='annelida'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[6]]<-ddply(subset(larv,phylum=='mollusca'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[7]]<-ddply(subset(larv,phylum=='cnidaria'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[8]]<-ddply(subset(larv,phylum=='ctenophora'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[9]]<-ddply(subset(larv,order=='gadiformes'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[10]]<-ddply(subset(larv,order=='stomiiformes'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[11]]<-ddply(subset(larv,order=='myctophiformes'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[12]]<-ddply(subset(larv,order=='pleuronectiformes'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
l[[13]]<-ddply(subset(larv,kingdom=='animalia'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
lrv<-rbind.fill(l)

#OUTPUT HERRING LARVAE DATA FOR THERMAL ENVELOPE FIGURES
ldat<-ddply(subset(larv,comname=='herring(atlantic)'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
ldat$occ<-factor(ifelse(ldat$stdno>0,TRUE,FALSE))
ldat$stdno<-ldat$stdno+.1
ldat<-subset(ldat,is.na(surftemp)==FALSE)
save(ldat,file='ldat.RData')











#CALCULATE RICHNESS, EVENNESS, FOR LARVAE
#NEED TO REMOVE 0'S OR CAN'T CALCULATE SHANNON INDEX
larv2<-subset(larv,stdno>0)
richfun<-function(d){
d<-na.omit(d)
#RICHNESS
dout<-data.frame(richness=length(unique(d$valid_name)))

#RICHNESS WITH 1% CUTOFF LOCALLY
tcells<-sum(d$stdno)
d$prp<-(d$stdno/tcells)*100
d2<-subset(d,prp>=1)
dout$richness1<-length(unique(d2$valid_name))

#SHANNON I
f2<-function(d3){ return(data.frame(p=sum(d3$stdno)/tcells)) }
sdat<-ddply(subset(d,select=c('valid_name','stdno')),.(valid_name),.fun=f2)
sdat$lp<-log(sdat$p)
sdat$pp<--1*(sdat$p*sdat$lp)
sdat<-subset(sdat,is.na(pp)==FALSE)

dout$shannon<-sum(sdat$pp,na.rm=TRUE)
dout$pe<-sum(sdat$pp)/(log(length(unique(d$valid_name))))
dout$totalabundance<-sum(d$stdno)
dout$nspecies<-length(unique(d$valid_name))

return(dout)
}
lrich<-ddply(subset(larv2,kingdom=='animalia',select=c('year','month','day','tday','lon','lat','tbin3','tbin5','tbin10','valid_name','bathy','stdno')),.(year,month,day,tday,tbin3,tbin5,tbin10,bathy,lon,lat),.fun=richfun,.progress='text')
names(lrich)[11:16]<-gsub(' ','',paste(names(lrich)[11:16],'.lrv'))



#ZOOPLANKTON: RETAIN LARGEST SPECIES
zoop<-read.csv("zoop_cpr_counts_spera_spawnar.csv",header=TRUE)
zoop<-subset(zoop,scientificname %in% c('Calanoida','Oithona','Calanus', 'Calanus finmarchicus','Calanus glacialis','Calanus hyperboreus','Candacia','Candacia armata','Centropages','Centropages hamatus','Centropages typicus','Pseudocalanus elongatus','Temora longicornis','Copepoda'))

#SUM COUNTS (INDIVIDUALS) PER SAMPLE
f<-function(d){  return(data.frame(zoop=sum(d$counts)))}
zoop<-ddply(zoop,.(lon,lat,year, month,djul,day,tbin3,tbin5,tbin10),.fun=f,.progress='text')










d<-subset(strat,select=c('s25','year','day','lon','lat','month'))
nm<-'Strat'
bp<-1982
segtrue<-FALSE


###########################################################################
##MODELING FUNCTIONS: FIRST FUNCTION TAKES DATA AND FIRST ESTIMATES ANNUAL AVERAGES, SECOND JUST TAKES ANNUAL AVERAGES
#TAKE ANNUAL TIMESERIES AS INPUT AND:
#1) DETERMINES WHICH MODEL FITS BEST BY AIC (LINEAR, NONLINEAR, SEGMENTED, OR BREAKPOINT
#2) ESTIMATES TOTAL CHANGE AS DIFFERENCE BETWEEN 1965-2015 AND 1985-2005
#3) OUTPUTS TOTAL CHANGE, BREAKPOINT (IF PRESENT), R2, PV, NYEARS
#4) MAKES TIMESERIES PLOT WITH BEST FITTING MODEL OVERLAID
modf<-function(d,nm,bp,segtrue){
    options(warn=-1)
    nmm<-names(d)[1]
    d<-subset(d,year>1960 & month>=8)
    names(d)[1]<-'y'
    d<-na.omit(d)
    #ESTIMATE ANNUAL AVERAGES
    if(nm %in% c('Lrv Flatf','Foram','Lrv Mollusc','Lrv Mycto')){
        print('lognormal')
        d$y<-d$y+.1
    mod<-gam(y~as.factor(year)+s(day,bs='cs',k=5)+s(lon,lat,k=30),data=d,gamma=1.4,family='gaussian'(link='log'))
    } else {
        mod<-gam(y~as.factor(year) + s(day,bs='cs',k=5) + s(lon,lat,k=30),data=d,gamma=1.4)
    }

    #PREDICT ANNUAL AVERAGES AND CIS
    lon<-mean(d$lon)
    lat<-mean(d$lat)
    day<-300
    pdat<-data.frame(year=sort(unique(d$year)))
    pdat<-subset(pdat,year>=1965)
    p<-predict(mod,newdata=data.frame(year=pdat$year,lon=lon,lat=lat,day=day),se.fit=TRUE,type='response')
    pdat$p<-p$fit
    pdat$se<-p$se.fit
    pdat$upr<-pdat$p+(1.96*pdat$se)
    pdat$lwr<-pdat$p-(1.96*pdat$se)
    pdat<-subset(pdat,p!=Inf)

    dum<-data.frame(year=sort(unique(d$year)),
                    n=tapply(d$y,d$year,length),
                    nmonths=tapply(d$month,d$year,function(x) length(unique(x))))
    pdat<-merge(pdat,dum,by=c('year'),all.x=TRUE,all.y=FALSE)
    if(nm %in% c('Nit','Sil','Phos','Lrv Fish','Strat','Temp 50m','Zoop Cells')){
            pdat<-subset(pdat,nmonths>0)
    } else { pdat<-subset(pdat,nmonths>1)
    }

    ylm<-c(min(pdat$lwr),max(pdat$upr))
    plot(pdat$year,pdat$p,pch=15,main=as.character(unique(d$db)),las=1,xlim=c(1965,2015),xaxt='n',ylim=ylm,xlab='Year',ylab=nm)
    mtext(nm,side=3,line=0,outer=FALSE)
    axis(1,at=seq(1965,2015,by=10))
    pdat$w<-1/pdat$se
    if(length(unique(d$year))<=10){dff<-4
    } else {dff<-5
    }

    #CHARACTERIZE TRENDS
    modl<-lm(p~year,data=pdat,weights=w)
    m<-breakpoints(p~year,data=pdat,h=3,breaks=1)
    modst<-lm(p~breakfactor(m,breaks=length(unique(m$breakpoints))),data=pdat,weights=w)
    modnl<-lm(p~bs(year,degree=3,df=dff),data=pdat,weights=w)
    if(segtrue==TRUE){
        modseg<-segmented(modl,seg.Z = ~year,weights=w,psi=bp,control=seg.control(it.max=200))
        dt<-data.frame(AIC(modl,modnl,modseg,modst))
    } else {rm(modseg)
            dt<-data.frame(AIC(modl,modnl,modst))
    }
    names(dt)<-ifelse(names(dt) %in% c('BIC'),'AIC',names(dt))
    dt$md<-rownames(dt)
    dt$AIC<-ifelse(dt$md=='modnl',dt$AIC+4,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modseg',dt$AIC-2,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modl',dt$AIC-2,dt$AIC)
    dt<-subset(dt,AIC==min(dt$AIC))

if(dt$md=='modnl'){
       modl<-modnl
    pdat2<-data.frame(year=seq(min(pdat$year),max(pdat$year),length.out=100))
} else if (dt$md=='modseg'){modl<-modseg
    pdat2<-data.frame(year=seq(min(pdat$year),max(pdat$year),length.out=100))
} else if (dt$md=='modst'){modl<-modst
    pdat2<-data.frame(year=sort(unique(pdat$year)))
} else { modl<-modl
    pdat2<-data.frame(year=seq(min(pdat$year),max(pdat$year),length.out=100))
}
    s<-summary(modl)
    r2<-round(s$r.squared,digits=2)
    p<-predict(modl,newdata=data.frame(year=pdat2$year,lon=lon,lat=lat,day=day),se.fit=TRUE,type='response')
    pdat2$p<-p$fit
    pdat2$se<-p$se.fit
    pdat2$upr<-pdat2$p+(1.96*pdat2$se)
    pdat2$lwr<-pdat2$p-(1.96*pdat2$se)
    polygon(c(pdat2$year,pdat2$year[length(pdat2$year):1]),c(pdat2$upr,pdat2$lwr[length(pdat2$lwr):1]),col=alpha('lightgray',.4),border=NA)
    points(pdat$year,pdat$p,pch=15)
    f<-function(dd){ lines(c(dd$year,dd$year),c(dd$lwr,dd$upr),col=alpha('black',.4)) }
    zz<-dlply(pdat,.(year),.fun=f)
    lines(pdat2$year,pdat2$p,col='red3')
    legend('topleft',paste('r2 =',r2),pch=16,col='white',bty='n')

if(dt$md!='modst'){
    pout<-data.frame(year=seq(min(pdat$year),max(pdat$year),1))
    pout$p<-predict(modl,newdata=data.frame(year=pout$year),type='response')
} else {
    p<-predict(modl,newdata=data.frame(year=sort(unique(pdat$year))),type='response')
    pout<-data.frame(year=seq(min(pdat$year),max(pdat$year),1))

    #GET BREAKPOINT
    dm<-data.frame(m$X)
    dm$id<-seq(1,dim(dm)[1],1)
    dm$bp<-breakfactor(m,breaks=length(unique(m$breakpoints)))
    bpt<-max(subset(dm,bp=='segment1')$year)
    pout$p<-ifelse(pout$year<=bpt,p[1],p[length(p)])
}

    #GET BREAKPOINT
    if(dt$md=='modst'){
    dm<-data.frame(m$X)
    dm$id<-seq(1,dim(dm)[1],1)
    dm$bp<-breakfactor(m,breaks=length(unique(m$breakpoints)))
    pout$bpt<-max(subset(dm,bp=='segment1')$year)
    } else if (dt$md=='modseg'){
    pout$bpt<-round(data.frame(modseg$psi)$Est.,digits=0)
    } else {
    pout$bpt<-NA
    }

    if(dt$md=='modst'){
    p<-predict(modl,type='response')
    t1.65<-p[1]
    t2.65<-p[length(p)]
    t1.05<-p[1]
    t2.05<-p[length(p)]
    } else {
    pdat3<-data.frame(year=c(1965,2015))
    p<-predict(modl,newdata=data.frame(year=pdat3$year),type='response',se.fit=TRUE)
    pdat3$p<-p$fit
    t1.65<-pdat3$p[1]
    t2.65<-pdat3$p[2]
    pdat3<-data.frame(year=c(1985,2005))
    p<-predict(modl,newdata=data.frame(year=pdat3$year),type='response',se.fit=TRUE)
    pdat3$p<-p$fit
    t1.05<-pdat3$p[1]
    t2.05<-pdat3$p[2]
    }

#################################################
#NOW ESTIMATE CHANGE IN ZUNITS
#CHARACTERIZE TRENDS
if(regexpr('CT', nm)>0 | regexpr('timing', nm)>0){
    pdat$pz<-(pdat$p-182)/sd(pdat$p,na.rm=TRUE)
} else {
    pdat$pz<-(pdat$p-mean(pdat$p,na.rm=TRUE))/sd(pdat$p,na.rm=TRUE)
}
    modl<-lm(pz~year,data=pdat,weights=w)
    m<-breakpoints(pz~year,data=pdat,h=3,breaks=1)
    modst<-lm(pz~breakfactor(m,breaks=length(unique(m$breakpoints))),data=pdat,weights=w)
    modnl<-lm(pz~bs(year,degree=3,df=dff),data=pdat,weights=w)
    if(segtrue==TRUE){
    modseg<-segmented(modl,seg.Z = ~year,weights=w,psi=bp,control=seg.control(it.max=200))
    dt<-data.frame(AIC(modl,modnl,modseg,modst))
    } else {rm(modseg)
    dt<-data.frame(AIC(modl,modnl,modst))
    }
    names(dt)<-ifelse(names(dt) %in% c('BIC'),'AIC',names(dt))
    dt$md<-rownames(dt)
    dt$AIC<-ifelse(dt$md=='modnl',dt$AIC+4,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modseg',dt$AIC-2,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modl',dt$AIC-2,dt$AIC)
    dt<-subset(dt,AIC==min(dt$AIC))

if(dt$md=='modnl'){
       modl<-modnl
} else if (dt$md=='modseg'){modl<-modseg
} else if (dt$md=='modst'){modl<-modst
} else { modl<-modl
}
    s<-summary(modl)
    pv<-round(s$coef[2,4],digits=4)
    r2<-round(s$r.squared,digits=2)
    if(dt$md=='modst'){
    p<-predict(modl,type='response')
    t1<-p[1]
    t2<-p[length(p)]
    } else {
    pdat3<-data.frame(year=c(1965,2015))
    p<-predict(modl,newdata=data.frame(year=pdat3$year),type='response',se.fit=TRUE)
    pdat3$p<-p$fit
    t1<-pdat3$p[1]
    t2<-pdat3$p[2]
    }

    pout<-subset(pout,select=c('year','p','bpt'))
    pout$t1<-t1
    pout$t2<-t2
    pout$t1.65<-t1.65
    pout$t2.65<-t2.65
    pout$t1.05<-t1.05
    pout$t2.05<-t2.05
    pout$zr2<-r2
    pout$zpv<-pv
    pout$n<-length(unique(pdat$year))
    names(pout)[2]<-gsub(' ','.',paste(names(pout)[2],'.',nm))
    return(pout)
    options(warn=0)
}



modf2<-function(d,nm,bp,ylbb,segtrue){
    options(warn=-1)
    nmm<-names(d)[1]
    d<-subset(d,year>1960)
    names(d)[1]<-'y'
    pdat<-na.omit(d)
    if(nmm=='sbb'){d$y<-d$y/10000
    } else NULL
    if(dim(pdat)[2]>=3 & nm %in% c('larv')){
        pdat$w<-1/pdat[,3]
    } else if(dim(pdat)[2]>=3){pdat$w<-pdat[,3]
    } else {pdat$w<-1
    }
    if(length(unique(pdat$year))<=10){dff<-4
    } else {dff<-5
    }

    #CHARACTERIZE TRENDS
    modl<-lm(y~year,data=pdat,weights=w)
    m<-breakpoints(y~year,data=pdat,h=3,breaks=1)
    modst<-lm(y~breakfactor(m,breaks=length(unique(m$breakpoints))),data=pdat,weights=w)
    modnl<-lm(y~bs(year,degree=3,df=dff),data=pdat,weights=w)
    if(segtrue==TRUE){
        modseg<-segmented(modl,seg.Z = ~year,psi=bp,control=seg.control(it.max=200),weights=w)
        dt<-data.frame(AIC(modl,modnl,modseg,modst))
    } else {rm(modseg)
            dt<-data.frame(AIC(modl,modnl,modst))
    }
    names(dt)<-ifelse(names(dt) %in% c('BIC'),'AIC',names(dt))
    dt$md<-rownames(dt)
    dt$AIC<-ifelse(dt$md=='modnl',dt$AIC+4,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modseg',dt$AIC-2,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modl',dt$AIC-2,dt$AIC)
    dt<-subset(dt,AIC==min(dt$AIC))

if(dt$md=='modnl'){
       modl<-modnl
    pdat2<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
} else if (dt$md=='modseg'){modl<-modseg
    pdat2<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
} else if (dt$md=='modst'){modl<-modst
    pdat2<-data.frame(year=sort(unique(pdat$year)))
} else { modl<-modl
    pdat2<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
}

    s<-summary(modl)
    r2<-round(s$r.squared,digits=2)
    p<-predict(modl,newdata=data.frame(year=pdat2$year),se.fit=TRUE,type='response')
    pdat2$p<-p$fit
    pdat2$se<-p$se.fit
    pdat2$upr<-pdat2$p+(1.96*pdat2$se)
    pdat2$lwr<-pdat2$p-(1.96*pdat2$se)
    ylm<-c(min(pdat$y),max(pdat$y))
    plot(pdat$year,pdat$y,pch=15,las=1,xlim=c(1965,2015),xaxt='n',ylim=ylm,xlab='Year',ylab=ylbb)
    mtext(nm,side=3,line=0,outer=FALSE)
    axis(1,at=seq(1965,2015,by=10))
    polygon(c(pdat2$year,pdat2$year[length(pdat2$year):1]),c(pdat2$upr,pdat2$lwr[length(pdat2$lwr):1]),col=alpha('lightgray',.4),border=NA)
    points(pdat$year,pdat$y,pch=15)
    f<-function(dd){ lines(c(dd$year,dd$year),c(dd$lwr,dd$upr),col=alpha('black',.4)) }
    zz<-dlply(pdat,.(year),.fun=f)
    lines(pdat2$year,pdat2$p,col='red3')
    legend('topright',paste('r2 =',r2),pch=16,col='white',bty='n')

    pd<-data.frame(year=c(1960,2010))
    p<-predict(modl,newdata=data.frame(year=pd$year),se.fit=TRUE,type='response')
    print(p$fit[1]-p$fit[2])


if(dt$md!='modst'){
    pout<-data.frame(year=seq(min(pdat$year),max(pdat$year),1))
    pout$p<-predict(modl,newdata=data.frame(year=pout$year),type='response')
} else {
    p<-predict(modl,newdata=data.frame(year=sort(unique(pdat$year))),type='response')
    pout<-data.frame(year=seq(min(pdat$year),max(pdat$year),1))

    #GET BREAKPOINT
    dm<-data.frame(m$X)
    dm$id<-seq(1,dim(dm)[1],1)
    dm$bp<-breakfactor(m,breaks=length(unique(m$breakpoints)))
    bpt<-max(subset(dm,bp=='segment1')$year)
    pout$p<-ifelse(pout$year<=bpt,p[1],p[length(p)])
}

    #GET BREAKPOINT
    if(dt$md=='modst'){
    dm<-data.frame(m$X)
    dm$id<-seq(1,dim(dm)[1],1)
    dm$bp<-breakfactor(m,breaks=length(unique(m$breakpoints)))
    pout$bpt<-max(subset(dm,bp=='segment1')$year)
    } else if (dt$md=='modseg'){
    pout$bpt<-round(data.frame(modseg$psi)$Est.,digits=0)
    } else {
    pout$bpt<-NA
    }

    if(dt$md=='modst'){
    p<-predict(modl,type='response')
    t1.65<-p[1]
    t2.65<-p[length(p)]
    t1.05<-p[1]
    t2.05<-p[length(p)]
    } else {
    pdat3<-data.frame(year=c(1965,2015))
    p<-predict(modl,newdata=data.frame(year=pdat3$year),type='response',se.fit=TRUE)
    pdat3$p<-p$fit
    t1.65<-pdat3$p[1]
    t2.65<-pdat3$p[2]
    pdat3<-data.frame(year=c(1985,2005))
    p<-predict(modl,newdata=data.frame(year=pdat3$year),type='response',se.fit=TRUE)
    pdat3$p<-p$fit
    t1.05<-pdat3$p[1]
    t2.05<-pdat3$p[2]
    }

#################################################
#NOW ESTIMATE CHANGE IN ZUNITS
 #CHARACTERIZE TRENDS
if(regexpr('CT', nm)>0 | regexpr('timing', nm)>0){
    pdat$pz<-(pdat$y-182)/sd(pdat$y,na.rm=TRUE)
} else {
    pdat$pz<-(pdat$y-mean(pdat$y,na.rm=TRUE))/sd(pdat$y,na.rm=TRUE)
}
    modl<-lm(pz~year,data=pdat,weights=w)
    m<-breakpoints(pz~year,data=pdat,h=3,breaks=1)
    modst<-lm(pz~breakfactor(m,breaks=length(unique(m$breakpoints))),data=pdat,weights=w)
    modnl<-lm(pz~bs(year,degree=3,df=dff),data=pdat,weights=w)
    if(segtrue==TRUE){
    modseg<-segmented(modl,seg.Z = ~year,weights=w,psi=bp,control=seg.control(it.max=200))
    dt<-data.frame(AIC(modl,modnl,modseg,modst))
    } else {rm(modseg)
    dt<-data.frame(AIC(modl,modnl,modst))
    }
    names(dt)<-ifelse(names(dt) %in% c('BIC'),'AIC',names(dt))
    dt$md<-rownames(dt)
    dt$AIC<-ifelse(dt$md=='modnl',dt$AIC+4,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modseg',dt$AIC-2,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modl',dt$AIC-2,dt$AIC)
    dt<-subset(dt,AIC==min(dt$AIC))

if(dt$md=='modnl'){
       modl<-modnl
} else if (dt$md=='modseg'){modl<-modseg
} else if (dt$md=='modst'){modl<-modst
} else { modl<-modl
}
    s<-summary(modl)
    pv<-round(s$coef[2,4],digits=4)
    r2<-round(s$r.squared,digits=2)
    if(dt$md=='modst'){
    p<-predict(modl,type='response')
    t1<-p[1]
    t2<-p[length(p)]
    } else {
    pdat3<-data.frame(year=c(1965,2015))
    p<-predict(modl,newdata=data.frame(year=pdat3$year),type='response',se.fit=TRUE)
    pdat3$p<-p$fit
    t1<-pdat3$p[1]
    t2<-pdat3$p[2]
    }

    pout<-subset(pout,select=c('year','p','bpt'))
    pout$t1<-t1
    pout$t2<-t2
    pout$t1.65<-t1.65
    pout$t2.65<-t2.65
    pout$t1.05<-t1.05
    pout$t2.05<-t2.05
    pout$zr2<-r2
    pout$zpv<-pv
    pout$n<-length(unique(pdat$year))
    names(pout)[2]<-gsub(' ','.',paste(names(pout)[2],'.',nm))
    return(pout)
    options(warn=0)
}



setwd(figsdir)
pdf('yearly_tsplots2.pdf',height=12,width=12)
par(mfrow=c(4,3),mar=c(4,4,1,1))
l<-list()
l[[1]]<-modf(subset(path,select=c('sst','year','day','lon','lat','month')),'SST',1995,TRUE)
l[[2]]<-modf(subset(path,select=c('wind','year','day','lon','lat','month')),'Wind',1984,FALSE)
l[[3]]<-modf(subset(strat,select=c('s25','year','day','lon','lat','month')),'Strat',1982,FALSE)
l[[4]]<-modf(subset(strat,select=c('temp50','year','day','lon','lat','month')),'Temp 50m',2005,TRUE)#R WITH SST: 76
l[[5]]<-modf(subset(lrv,var=='T L Fish',select=c('stdno','year','day','lon','lat','month')),'Lrv Fish',1980,FALSE)
l[[6]]<-modf(subset(plank,select=c('copepods.z','year','day','lon','lat','month')),'Copepods',1995,FALSE)
l[[7]]<-modf(subset(nit,select=c('nit','year','day','lon','lat','month')),'Nit',1982,FALSE)
l[[8]]<-modf(subset(sil,select=c('sil','year','day','lon','lat','month')),'Sil',1992,FALSE)#
l[[9]]<-modf(subset(phos,select=c('phos','year','day','lon','lat','month')),'Phos',1985,FALSE)
l[[10]]<-modf(subset(plank,select=c('shannon.p','year','day','lon','lat','month')),'Phyt Div',1995,TRUE)
l[[11]]<-modf(subset(plank,select=c('totalabundance.p','year','day','lon','lat','month')),'Phyt Cells',1995,TRUE)
l[[12]]<-modf(subset(plank,select=c('richness.p','year','day','lon','lat','month')),'Phyt Rich',1996,TRUE)
l[[13]]<-modf(subset(plank,select=c('pe.z','year','day','lon','lat','month')),'Zoop Even',1995,FALSE)
l[[14]]<-modf(subset(plank,select=c('pe.p','year','day','lon','lat','month')),'Phyt Even',1986,FALSE)
l[[15]]<-modf2(subset(vtemp,select=c('daypeak','year')),'SST Timing',1995,'Day',TRUE)
l[[16]]<-modf2(subset(vtemp,select=c('amp','year')),'SST Amp',1995,'SST',TRUE)
l[[17]]<-modf2(subset(vtemp,select=c('dur','year')),'SST Dur',1995,'Days',TRUE)
l[[18]]<-modf2(subset(vtemp,select=c('mx','year')),'SST Seas max',1985,'SST',FALSE)
l[[19]]<-modf2(subset(vwind,select=c('daypeak','year')),'Wind Timing',1995,'Day',TRUE)
l[[20]]<-modf2(subset(vwind,select=c('amp','year')),'Wind Amp',1990,'Wind',FALSE)
l[[21]]<-modf2(subset(vwind,select=c('dur','year')),'Wind Dur',1995,'Days',FALSE)
l[[22]]<-modf2(subset(vwind,select=c('mx','year')),'Wind Seas max',1995,'Wind',FALSE)
l[[23]]<-modf2(subset(pto,select=c('turnover','year')),'Phyt TO',1990,'Phyt Turn',TRUE)
l[[24]]<-modf2(subset(zto,select=c('turnover','year')),'Zoop TO',1995,'Zoop Turn',FALSE)
l[[25]]<-modf2(subset(ctdat,var=='Chl.mn',select=c('ct','year','nmonth')),'Chl CT',1995,'CT Chl',TRUE)
l[[26]]<-modf2(subset(ctdat,var=='F Larvae',select=c('ct','year','nmonth')),'Lrv Fish CT',1980,'CT Larvae',FALSE)
l[[27]]<-modf2(subset(ctdat,var=='P diversity',select=c('ct','year','nmonth')),'Phyt Div CT',1995,'CT P diversity',TRUE)
l[[28]]<-modf2(subset(ctdat,var=='P evenness',select=c('ct','year','nmonth')),'Phyt Even CT',1995,'CT P evenness',TRUE)
l[[29]]<-modf2(subset(ctdat,var=='P richness',select=c('ct','year','nmonth')),'Phyt Rich CT',1995,'CT P richness',TRUE)
l[[30]]<-modf2(subset(ctdat,var=='P total',select=c('ct','year','nmonth')),'Phyt Cells CT',1995,'CT P P total',TRUE)
l[[31]]<-modf2(subset(ctdat,var=='Stratification',select=c('ct','year','nmonth')),'Strat CT',1995,'CT Stratification',TRUE)
l[[32]]<-modf2(subset(ctdat,var=='Temperature',select=c('ct','year','nmonth')),'Temp 50 CT',1998,'CT Temp 50',FALSE)
l[[33]]<-modf2(subset(ctdat,var=='Zooplankton',select=c('ct','year','nmonth')),'Zoop CT',1970,'CT Zoop',FALSE)

l[[34]]<-modf2(subset(ctdat,var=='L Diversity',select=c('ct','year','nmonth')),'Lrv Fish Div CT',1970,'CT L Diversity',FALSE)
l[[35]]<-modf2(subset(ctdat,var=='L Richness',select=c('ct','year','nmonth')),'Lrv Fish Rich CT',1970,'CT L Richness',FALSE)
l[[36]]<-modf2(subset(ctdat,var=='L Evenness',select=c('ct','year','nmonth')),'Lrv Fish Even CT',1970,'CT L Evenness',FALSE)
l[[37]]<-modf(subset(plank,select=c('totalabundance.z','year','day','lon','lat','month')),'Zoop Cells',1995,FALSE)
l[[38]]<-modf2(subset(ctdat,var=='GS Dist',select=c('ct','year','nmonth')),'GS Dist CT',1970,'CT GS Dist',FALSE)
l[[39]]<-modf2(subset(ctdat,var=='SS Dist',select=c('ct','year','nmonth')),'SS Dist CT',1970,'CT SS Dist',FALSE)

l[[40]]<-modf(subset(lrv,var=='Herring',select=c('stdno','year','day','lon','lat','month')),'Herr Lrv',1995,TRUE)
l[[41]]<-modf(subset(lrv,var=='Arthropod',select=c('stdno','year','day','lon','lat','month')),'Lrv Arthro',1995,TRUE)
l[[42]]<-modf(subset(lrv,var=='Foraminifera',select=c('stdno','year','day','lon','lat','month')),'Foram',1995,TRUE)
l[[43]]<-modf(subset(lrv,var=='Chaetognath',select=c('stdno','year','day','lon','lat','month')),'Chaeto',1995,TRUE)
l[[44]]<-modf(subset(lrv,var=='Annelids',select=c('stdno','year','day','lon','lat','month')),'Lrv Annelid',1995,TRUE)
l[[45]]<-modf(subset(lrv,var=='Mollusc',select=c('stdno','year','day','lon','lat','month')),'Lrv Moll',1995,FALSE)
l[[46]]<-modf(subset(lrv,var=='Cnidarian',select=c('stdno','year','day','lon','lat','month')),'Lrv Cnid',1995,TRUE)
l[[47]]<-modf(subset(lrv,var=='Ctenophore',select=c('stdno','year','day','lon','lat','month')),'Lrv Cteno',1995,TRUE)
l[[48]]<-modf(subset(lrv,var=='Gadids',select=c('stdno','year','day','lon','lat','month')),'Lrv Gadids',1990,FALSE)
l[[49]]<-modf(subset(lrv,var=='Dragonfish',select=c('stdno','year','day','lon','lat','month')),'Lrv Dragon',1985,FALSE)
l[[50]]<-modf(subset(lrv,var=='Myctophids',select=c('stdno','year','day','lon','lat','month')),'Lrv Mycto',1995,FALSE)
l[[51]]<-modf(subset(lrv,var=='Flatfish',select=c('stdno','year','day','lon','lat','month')),'Lrv Flatf',1985,FALSE)
l[[52]]<-modf2(subset(cur,select=c('GS.Dist','year','GS.Dist.sd')),'GS Dist',1985,'GS Distance',TRUE)
l[[53]]<-modf2(subset(cur,select=c('SS.Dist','year','SS.Dist.sd')),'SS Dist',1985,'SS Distance',TRUE)
l[[54]]<-modf2(subset(her,select=c('r','year')),'Her Recruitment',1985,'Her Recruitment',TRUE)
l[[55]]<-modf2(subset(her,select=c('s','year')),'Her Survivorship',1985,'Her Survivorship',TRUE)
l[[56]]<-modf2(subset(her,select=c('f','year')),'Her F',1985,'Her F',TRUE)
l[[57]]<-modf(subset(plank,select=c('richness.z','year','day','lon','lat','month')),'Zoop Rich',1965,FALSE)
l[[58]]<-modf2(subset(her,select=c('obs.w','year','n.obs.inrange')),'Her CF obs',1985,'Her CF obs',TRUE)
l[[59]]<-modf2(subset(her,select=c('ssb','year')),'Her SSB',1985,'Her SSB',TRUE)
l[[60]]<-modf2(subset(her,select=c('wgt','year')),'Her Weight',1985,'Her Weight',TRUE)
l[[61]]<-modf2(subset(nao,select=c('nao','year','nao.sd')),'NAO',1985,'NAO',TRUE)
l[[62]]<-modf(subset(lrich,select=c('richness.lrv','year','day','lon','lat','month')),'Lrv Fish Rich',1985,TRUE)
l[[63]]<-modf(subset(lrich,select=c('shannon.lrv','year','day','lon','lat','month')),'Lrv Fish Div',1985,TRUE)
l[[64]]<-modf(subset(lrich,select=c('pe.lrv','year','day','lon','lat','month')),'Lrv Fish Even',1985,TRUE)
l[[65]]<-modf(subset(plank,select=c('shannon.z','year','day','lon','lat','month')),'Zoop Div',1995,FALSE)

dev.off()



#GETS RATES OF CHANGE FOR METAPLOT
f<-function(d){return(unique(data.frame(t1=d[,4],
                                        t2=d[,5],
                                        t1.65=d[,6],
                                        t2.65=d[,7],
                                        t1.05=d[,8],
                                        t2.05=d[,9],
                                        r2=d[,10],
                                        pv=d[,11],
                                        n=d[,12],
                                 var=names(d)[2])))}
mpdat<-ldply(l,.fun=f)
mpdat$var<-gsub('\\.',' ',mpdat$var)
mpdat$var<-gsub('p  ',' ',mpdat$var)

mpdat$pct.65<-ifelse(mpdat$t1.65<mpdat$t2.65,round(((mpdat$t2.65-mpdat$t1.65)/mpdat$t1.65)*100,digits=2),
                  round(-1*((mpdat$t1.65-mpdat$t2.65)/mpdat$t1.65)*100,digits=2))
mpdat$pct.05<-ifelse(mpdat$t1.05<mpdat$t2.05,round(((mpdat$t2.05-mpdat$t1.05)/mpdat$t1.05)*100,digits=2),
                  round(-1*((mpdat$t1.05-mpdat$t2.05)/mpdat$t1.05)*100,digits=2))
mpdat0<-subset(mpdat,is.na(r2)==FALSE)


#GETS BREAKPOINTS IF PRESENT
f<-function(d){return(unique(data.frame(bpt=d[,3],
                                 var=names(d)[2])))}
bpdat<-ldply(l,.fun=f)
bpdat$var<-gsub('\\.',' ',bpdat$var)
bpdat$var<-gsub('p  ',' ',bpdat$var)
bpdat<-subset(bpdat,is.na(bpt)==FALSE)

#COMBINE YEARLY PREDICITONS IN TO SINGLE DB
f<-function(d){return(d[,1:2])}
z<-llply(l,.fun=f)

#FORMAT AS DATA.FRAME
pdat.lin<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), z)#COMBINE
rownames(pdat.lin)<-pdat.lin$year
pdat.lin<-pdat.lin[,-1]
names(pdat.lin)<-gsub('\\.',' ',names(pdat.lin))
names(pdat.lin)<-gsub('p  ',' ',names(pdat.lin))


#GET ORDER FOR IMAGE PLOT - SORT FROM RAPID DECLINE TO RAPID INCREASE
p<-pdat.lin
p$year<-as.numeric(rownames(pdat.lin))
vbls<-names(pdat.lin)
ll<-list()
for(i in 1:length(vbls)){
    print(vbls[i])
    d<-na.omit(subset(p,select=c('year',paste(vbls[i]))))
    names(d)[2]<-'y'
    d$mn<-mean(d$y,na.rm=TRUE)
    d$sdd<-sd(d$y,na.rm=TRUE)
    d$z<-(d$y-d$mn)/d$sdd

    d1<-subset(d,year==min(d$year))
    d2<-subset(d,year==max(d$year))
    a<-data.frame(var=vbls[i],
                  chng=d1$z-d2$z)
    ll[[i]]<-a
}
ord<-data.frame(do.call('rbind',ll))
ord<-ord[order(ord$chng),]


#FORMAT AS MATRIX
vbls<-as.character(ord$var)
ll<-list()
for(i in 1:length(vbls)){
    d<-subset(pdat.lin,select=c(paste(vbls[i])))
    names(d)[1]<-'y'
    d$mn<-mean(d$y,na.rm=TRUE)
    d$sdd<-sd(d$y,na.rm=TRUE)
    d$z<-(d$y-d$mn)/d$sdd
    dd<-data.frame(year=rownames(pdat.lin),
                   lbl=vbls[i],
                   vbl=i,
                   z=d$z)
    ll[[i]]<-dd
}
mat.lin<-rbind.fill(ll)
mat.lin$year<-as.numeric(as.character(mat.lin$year))


#FORMAT AS MATRIX
bpdat<-bpdat[order(bpdat$bpt,decreasing=TRUE),]
vbls<-as.character(bpdat$var)
ll<-list()
for(i in 1:length(vbls)){
    print(vbls[i])
    d<-subset(pdat.lin,select=c(paste(vbls[i])))
    names(d)[1]<-'y'
    d$mn<-mean(d$y,na.rm=TRUE)
    d$sdd<-sd(d$y,na.rm=TRUE)
    d$z<-(d$y-d$mn)/d$sdd
    dd<-data.frame(year=rownames(pdat.lin),
                       vbl=i,
                       z=d$z,
                   lbl=vbls[i])
    ll[[i]]<-dd
}
mat.lin.bp<-rbind.fill(ll)
mat.lin.bp$year<-as.numeric(as.character(mat.lin.bp$year))


#ADD PLOTTING COORDINATES TO BREAKPOINTS DATA
dm<-unique(subset(mat.lin.bp,select=c('vbl','lbl')))
bpdat<-merge(bpdat,dm,by.x='var',by.y='lbl',all.x=TRUE,all.y=FALSE)
bpdat$bpt<-bpdat$bpt+.5

ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }


pltfun<-function(a){
aa<-data.frame(z=seq((-max(abs(a$z),na.rm=TRUE)),max(abs(a$z),na.rm=TRUE),length.out=100))
a<-rbind.fill(a,aa)
n<-21
brks<-seq((min(aa$z,na.rm=TRUE)-.01),max(aa$z,na.rm=TRUE)+.01,length.out=n)
brks2<-round(seq((min(aa$z,na.rm=TRUE)-.01),max(aa$z,na.rm=TRUE)+.01,length.out=n),digits=2)
a$zcat<-cut(a$z,breaks=brks)
a$zcat2<-cut(a$z,breaks=brks2)
lbls<-sort(unique(a$zcat2))
#cls<-colorRampPalette(brewer.pal(name='RdYlBu',8))
aa<-subset(a,is.na(z)==FALSE)
ylm<-seq(min(a$vbl,na.rm=TRUE),max(a$vbl,na.rm=TRUE),1)
lb<-unique(subset(a,select=c('vbl','lbl')))
lb<-na.omit(lb[order(lb$vbl),])

return(
    ggplot()+
    geom_tile(data=a, aes(x=year, y=vbl,fill=zcat),colour='white',size=.001)+
    geom_tile(data=aa, aes(x=year, y=vbl,fill=NA),colour='black',size=.1)+
    scale_fill_manual(breaks=as.character(lbls),values=matlab.like(length(lbls)),labels=lbls,na.value="transparent",guide=guide_legend(title=paste('Z-score'))) +
    scale_alpha(guide = 'none')+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(.2, "in"),legend.text=element_text(size=6))+
    scale_y_continuous(expand=c(0,0),breaks=ylm,labels=lb$lbl,limits=c(0,max(ylm)+1))
)
}

setwd(figsdir)
pdf('linear_trends_image.pdf',height=7,width=9)
pltfun(subset(mat.lin.bp,year>=1965))+ geom_point(aes(x=bpt,y=vbl),data=bpdat,col='black',alpha=1,size=3,shape=23)

dev.off()

setwd(figsdir)
pdf('linear_trends_image2.pdf',height=9,width=8)
pltfun(subset(mat.lin,year>=1965))

pltfun(subset(mat.lin,year>=1965))+ geom_point(aes(x=bpt,y=vbl),data=bpdat,col='black',alpha=1,size=3,shape=23)
dev.off()










#########################################################

#########################################################
mpdat<-mpdat0
dum<-data.frame(t1=mpdat$t1+150,
                t2=mpdat$t2+150)
mpdat$chng<-dum$t2-dum$t1


brks<-seq(-0.001,1,.1)
mpdat$r2f<-cut(mpdat$r2,breaks=brks)
dm<-data.frame(r2f=sort(unique(mpdat$r2f)),
               cl=matlab.like(length(unique(mpdat$r2f))),
               cl2=rev(gray.colors(length(unique(mpdat$r2f)))),
               alph=seq(0,1,length.out=10))
mpdat<-merge(mpdat,dm,by=c('r2f'),all.x=TRUE,all.y=FALSE)


mpdat$grp<-ifelse(regexpr('Her', mpdat$var)>0,'e_her','a_env')
mpdat$grp<-ifelse(regexpr('Lrv', mpdat$var)>0,'d_lrv',mpdat$grp)
mpdat$grp<-ifelse(regexpr('Zoop', mpdat$var)>0 |
                  regexpr('Foram', mpdat$var)>0 |
                  regexpr('Copepods', mpdat$var)>0 |
                  regexpr('Chaeto', mpdat$var)>0,'c_zoop',mpdat$grp)
mpdat$grp<-ifelse(regexpr('Phyt', mpdat$var)>0 |
                 regexpr('Chl', mpdat$var)>0,'b_phyt',mpdat$grp)

f<-function(d){    d<-d[order(d$chng),]}
mpdat<-ddply(mpdat,.(grp),.fun=f)
mpdat$mvar<-seq(1,dim(mpdat)[1],1)
cx<-1.25

setwd(figsdir)
pdf('pedicted_model_changes.pdf',height=10,width=9)
par(mar=c(4,8,1,1))
plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-5,15),col=as.character(mpdat$cl),pch=16,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1965-2015 Change [Standard Deviations]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(10,5,zr,col=as.character(dm$cl),horizontal=TRUE,strip.width=.02,strip.length=.2)
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl),pch=16,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=1,cex=cx,lwd=.01,col='gray')

plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-5,15),col=as.character(mpdat$cl2),pch=16,yaxt='n',ylab='',cex=cx)
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(12,10,zr,col=rev(gray.colors(10)),horizontal=TRUE,strip.width=.04,strip.length=.3)
abline(h=lns+.5,col='gray',lwd=2)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl2),pch=16,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=1,cex=cx,lwd=.1)



#1965 TO 2015
mpdat<-mpdat0
mpdat$chng<-mpdat$pct.65

brks<-seq(-0.001,1,.1)
mpdat$r2f<-cut(mpdat$r2,breaks=brks)
dm<-data.frame(r2f=sort(unique(mpdat$r2f)),
               cl=matlab.like(length(unique(mpdat$r2f))),
               cl2=rev(gray.colors(length(unique(mpdat$r2f)))),
               alph=seq(0,1,length.out=10))
mpdat<-merge(mpdat,dm,by=c('r2f'),all.x=TRUE,all.y=FALSE)
mpdat$pchh<-ifelse(mpdat$pv<=.05,16,17)
mpdat$pchh2<-ifelse(mpdat$pv<=.05,21,2)

mpdat$grp<-ifelse(regexpr('Her', mpdat$var)>0,'e_her','a_env')
mpdat$grp<-ifelse(regexpr('Lrv', mpdat$var)>0,'d_lrv',mpdat$grp)
mpdat$grp<-ifelse(regexpr('Zoop', mpdat$var)>0 |
                  regexpr('Foram', mpdat$var)>0 |
                  regexpr('Copepods', mpdat$var)>0 |
                  regexpr('Chaeto', mpdat$var)>0,'c_zoop',mpdat$grp)
mpdat$grp<-ifelse(regexpr('Phyt', mpdat$var)>0 |
                 regexpr('Chl', mpdat$var)>0,'b_phyt',mpdat$grp)

f<-function(d){    d<-d[order(d$chng),]}
mpdat<-ddply(mpdat,.(grp),.fun=f)
mpdat$mvar<-seq(1,dim(mpdat)[1],1)
cx<-1.25

par(mar=c(4,8,1,1))
plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-400,300),col=as.character(mpdat$cl),pch=mpdat$pchh,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1965-2015 Change [Percent]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(200,5,zr,col=as.character(dm$cl),horizontal=TRUE,strip.width=.02,strip.length=.2)
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl),pch=mpdat$pchh,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=mpdat$pchh2,cex=cx,lwd=.01,col='gray')
legend(200,10,legend=c('P-value<=0.05','P-value>0.05'),pch=c(16,17),bty='n')

par(mar=c(4,8,1,1))
plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-400,300),col=as.character(mpdat$cl2),pch=mpdat$pchh,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1965-2015 Change [Percent]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(200,5,zr,col=as.character(dm$cl2),horizontal=TRUE,strip.width=.02,strip.length=.2)
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl2),pch=mpdat$pchh,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=mpdat$pchh2,cex=cx,lwd=.01,col='gray')
legend(150,10,legend=c('P-value<=0.05','P-value>0.05'),pch=c(16,17),bty='n')




###1980-2005
mpdat<-mpdat0
mpdat$chng<-mpdat$pct.05

brks<-seq(-0.001,1,.1)
mpdat$r2f<-cut(mpdat$r2,breaks=brks)
dm<-data.frame(r2f=sort(unique(mpdat$r2f)),
               cl=matlab.like(length(unique(mpdat$r2f))),
               cl2=rev(gray.colors(length(unique(mpdat$r2f)))),
               alph=seq(0,1,length.out=10))
mpdat<-merge(mpdat,dm,by=c('r2f'),all.x=TRUE,all.y=FALSE)
mpdat$pchh<-ifelse(mpdat$pv<=.05,16,17)
mpdat$pchh2<-ifelse(mpdat$pv<=.05,21,2)

mpdat$grp<-ifelse(regexpr('Her', mpdat$var)>0,'e_her','a_env')
mpdat$grp<-ifelse(regexpr('Lrv', mpdat$var)>0,'d_lrv',mpdat$grp)
mpdat$grp<-ifelse(regexpr('Zoop', mpdat$var)>0 |
                  regexpr('Foram', mpdat$var)>0 |
                  regexpr('Copepods', mpdat$var)>0 |
                  regexpr('Chaeto', mpdat$var)>0,'c_zoop',mpdat$grp)
mpdat$grp<-ifelse(regexpr('Phyt', mpdat$var)>0 |
                 regexpr('Chl', mpdat$var)>0,'b_phyt',mpdat$grp)

f<-function(d){    d<-d[order(d$chng),]}
mpdat<-ddply(mpdat,.(grp),.fun=f)
mpdat$mvar<-seq(1,dim(mpdat)[1],1)
cx<-1.25

par(mar=c(4,8,1,1))
plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-100,100),col=as.character(mpdat$cl),pch=mpdat$pchh,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1985-2005 Change [Percent]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl),pch=mpdat$pchh,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=mpdat$pchh2,cex=cx,lwd=.01,col='gray')
colorbar.plot(50,5,zr,col=as.character(dm$cl),horizontal=TRUE,strip.width=.02,strip.length=.2)
legend(50,10,legend=c('P-value<=0.05','P-value>0.05'),pch=c(16,17),bty='n')
md<-subset(mpdat,abs(chng)>150)
md$x1<-ifelse(md$chng<150,-110,110)
f<-function(d){
    if(d$x1<0){
        arrows(d$x1+6,d$mvar,d$x1+2,d$mvar,length=.1)
    } else {
        arrows(d$x1-6,d$mvar,d$x1-2,d$mvar,length=.1)
    }
}
zz<-dlply(md,.(mvar),.fun=f)


par(mar=c(4,8,1,1))
plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-100,100),col=as.character(mpdat$cl2),pch=mpdat$pchh,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1985-2005 Change [Percent]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(50,5,zr,col=as.character(dm$cl2),horizontal=TRUE,strip.width=.02,strip.length=.2)
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl2),pch=mpdat$pchh,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=mpdat$pchh2,cex=cx,lwd=.01,col='gray')
legend(50,10,legend=c('P-value<=0.05','P-value>0.05'),pch=c(16,17),bty='n')
md<-subset(mpdat,abs(chng)>150)
md$x1<-ifelse(md$chng<150,-110,110)
f<-function(d){
    if(d$x1<0){
        arrows(d$x1+6,d$mvar,d$x1+2,d$mvar,length=.1)
    } else {
        arrows(d$x1-6,d$mvar,d$x1-2,d$mvar,length=.1)
    }
}
zz<-dlply(md,.(mvar),.fun=f)

dev.off()










####################################################################
#PREDICTS ANNUAL AVERAGES
modfac<-function(d,nm,bp,nmcols){
    nmm<-names(d)[1]
    d<-subset(d,year>1960 & month>=8)
    names(d)[1]<-'y'
    d<-na.omit(d)
    if(nm %in% c('dia','pro','slc')){
        print('quasi')
        mod<-gam(y~as.factor(year) + s(day,bs='cs',k=5) + s(lon,lat,k=30),data=d,gamma=1.4,family='quasibinomial')
    } else {
        mod<-gam(y~as.factor(year) + s(day,bs='cs',k=5) + s(lon,lat,k=30),data=d,gamma=1.4)
    }

    lon<-mean(d$lon)
    lat<-mean(d$lat)
    day<-300
    pdat<-data.frame(year=sort(unique(d$year)))
    pdat<-subset(pdat,year>=1965)
    p<-predict(mod,newdata=data.frame(year=pdat$year,lon=lon,lat=lat,day=day),se.fit=TRUE,type='response')
    pdat$p<-p$fit
    pdat$w<-1/p$se.fit

    dum<-data.frame(year=sort(unique(d$year)),
                    n=tapply(d$y,d$year,length),
                    nmonths=tapply(d$month,d$year,function(x) length(unique(x))))
    pdat<-merge(pdat,dum,by=c('year'),all.x=TRUE,all.y=FALSE)
    if(nm %in% c('Stratification','Temperature 50m','Wind','SST')){
            pdat<-subset(pdat,nmonths>1)
    } else { pdat<-subset(pdat,nmonths>0)
    }

    pdat<-subset(pdat,select=c('year','p','w'))


    if(nmcols==TRUE){
        pdat<-subset(pdat,select=c('year','p'))
       names(pdat)[2]<-gsub(' ','.',paste(names(pdat)[2],'.',nm))
#        names(pdat)[2]<-gsub(' ','',paste(names(pdat)[2],'.',nm))
    } else {pdat$var<-nm
    }

return(pdat)
}



modfac2<-function(d,nm,x1,x2,nmcols){
    nmm<-names(d)[1]
    d<-subset(d,year>1960)
    names(d)[1]<-'p'
    pdat<-na.omit(d)
    if(nmm=='sbb'){d$p<-d$y/10000
    } else NULL

    if(dim(pdat)[2]>=3 & nmm=='larv'){
        pdat$w<-1/pdat[,3]
    } else if(dim(pdat)[2]>=3){
        pdat$w<-pdat[,3]
    } else {
        pdat$w<-1
    }

    if(nm %in% c('CT Nutrients')){
        pdat<-data.frame(year=sort(unique(pdat$year)),
                         p=tapply(pdat$p,pdat$year,mean),
                         w=tapply(pdat$w,pdat$year,mean))
    } else NULL

    pdat<-subset(pdat,select=c('year','p','w'))

    if(nmcols==TRUE){
        pdat<-subset(pdat,select=c('year','p'))
       names(pdat)[2]<-gsub(' ','.',paste(names(pdat)[2],'.',nm))
    } else {pdat$var<-nm
    }
return(pdat)
}


ll<-list()
ll[[1]]<-modfac(subset(path,select=c('sst','year','day','lon','lat','month')),'SST',1995,TRUE)
ll[[2]]<-modfac(subset(path,select=c('wind','year','day','lon','lat','month')),'Wind',1984,TRUE)
ll[[3]]<-modfac(subset(strat,select=c('s25','year','day','lon','lat','month')),'Strat',1982,TRUE)
ll[[4]]<-modfac(subset(strat,select=c('temp50','year','day','lon','lat','month')),'Temp 50m',2005,TRUE)#R WITH SST: 76
ll[[5]]<-modfac(subset(lrv,var=='T L Fish',select=c('stdno','year','day','lon','lat','month')),'Lrv Fish',1980,TRUE)
ll[[6]]<-modfac(subset(plank,select=c('copepods.z','year','day','lon','lat','month')),'Copepods',1995,TRUE)
ll[[7]]<-modfac(subset(nit,select=c('nit','year','day','lon','lat','month')),'Nit',1982,TRUE)
ll[[8]]<-modfac(subset(sil,select=c('sil','year','day','lon','lat','month')),'Sil',1992,TRUE)#
ll[[9]]<-modfac(subset(phos,select=c('phos','year','day','lon','lat','month')),'Phos',1985,TRUE)
ll[[10]]<-modfac(subset(plank,select=c('shannon.p','year','day','lon','lat','month')),'Phyt Div',1995,TRUE)
ll[[11]]<-modfac(subset(plank,select=c('totalabundance.p','year','day','lon','lat','month')),'Phyt Cells',1995,TRUE)
ll[[12]]<-modfac(subset(plank,select=c('richness.p','year','day','lon','lat','month')),'Phyt Rich',1996,TRUE)
ll[[13]]<-modfac(subset(plank,select=c('pe.z','year','day','lon','lat','month')),'Zoop Even',1995,TRUE)
ll[[14]]<-modfac(subset(plank,select=c('pe.p','year','day','lon','lat','month')),'Phyt Even',1986,TRUE)
ll[[15]]<-modfac2(subset(vtemp,select=c('daypeak','year')),'SST Timing',1995,'Day',TRUE)
ll[[16]]<-modfac2(subset(vtemp,select=c('amp','year')),'SST Amp',1995,'SST',TRUE)
ll[[17]]<-modfac2(subset(vtemp,select=c('dur','year')),'SST Dur',1995,'Days',TRUE)
ll[[18]]<-modfac2(subset(vtemp,select=c('mx','year')),'SST Seas max',1985,'SST',TRUE)
ll[[19]]<-modfac2(subset(vwind,select=c('daypeak','year')),'Wind Timing',1995,'Day',TRUE)
ll[[20]]<-modfac2(subset(vwind,select=c('amp','year')),'Wind Amp',1990,'Wind',TRUE)
ll[[21]]<-modfac2(subset(vwind,select=c('dur','year')),'Wind Dur',1995,'Days',TRUE)
ll[[22]]<-modfac2(subset(vwind,select=c('mx','year')),'Wind Seas max',1995,'Wind',TRUE)
ll[[23]]<-modfac2(subset(pto,select=c('turnover','year')),'Phyt TO',1990,'Phyt Turn',TRUE)
ll[[24]]<-modfac2(subset(zto,select=c('turnover','year')),'Zoop TO',1995,'Zoop Turn',TRUE)
ll[[25]]<-modfac2(subset(ctdat,var=='Chl.mn',select=c('ct','year','nmonth')),'Chl CT',1995,'CT Chl',TRUE)
ll[[26]]<-modfac2(subset(ctdat,var=='F Larvae',select=c('ct','year','nmonth')),'Lrv Fish CT',1980,'CT Larvae',TRUE)
ll[[27]]<-modfac2(subset(ctdat,var=='P diversity',select=c('ct','year','nmonth')),'Phyt Div CT',1995,'CT P diversity',TRUE)
ll[[28]]<-modfac2(subset(ctdat,var=='P evenness',select=c('ct','year','nmonth')),'Phyt Even CT',1995,'CT P evenness',TRUE)
ll[[29]]<-modfac2(subset(ctdat,var=='P richness',select=c('ct','year','nmonth')),'Phyt Rich CT',1995,'CT P richness',TRUE)
ll[[30]]<-modfac2(subset(ctdat,var=='P total',select=c('ct','year','nmonth')),'Phyt Cells CT',1995,'CT P P total',TRUE)
ll[[31]]<-modfac2(subset(ctdat,var=='Stratification',select=c('ct','year','nmonth')),'Strat CT',1995,'CT Stratification',TRUE)
ll[[32]]<-modfac2(subset(ctdat,var=='Temperature',select=c('ct','year','nmonth')),'Temp 50 CT',1998,'CT Temp 50',TRUE)
ll[[33]]<-modfac2(subset(ctdat,var=='Zooplankton',select=c('ct','year','nmonth')),'Zoop CT',1970,'CT Zoop',TRUE)

ll[[34]]<-modfac2(subset(ctdat,var=='L Diversity',select=c('ct','year','nmonth')),'Lrv Fish Div CT',1970,'CT L Diversity',TRUE)
ll[[35]]<-modfac2(subset(ctdat,var=='L Richness',select=c('ct','year','nmonth')),'Lrv Fish Rich CT',1970,'CT L Richness',TRUE)
ll[[36]]<-modfac2(subset(ctdat,var=='L Evenness',select=c('ct','year','nmonth')),'Lrv Fish Even CT',1970,'CT L Evenness',TRUE)
ll[[37]]<-modfac(subset(plank,select=c('totalabundance.z','year','day','lon','lat','month')),'Zoop Cells',1995,TRUE)
ll[[38]]<-modfac2(subset(ctdat,var=='GS Dist',select=c('ct','year','nmonth')),'GS Dist CT',1970,'CT GS Dist',TRUE)
ll[[39]]<-modfac2(subset(ctdat,var=='SS Dist',select=c('ct','year','nmonth')),'SS Dist CT',1970,'CT SS Dist',TRUE)

ll[[40]]<-modfac(subset(lrv,var=='Herring',select=c('stdno','year','day','lon','lat','month')),'Herr Lrv',1995,TRUE)
ll[[41]]<-modfac(subset(lrv,var=='Arthropod',select=c('stdno','year','day','lon','lat','month')),'Lrv Arthro',1995,TRUE)
ll[[42]]<-modfac(subset(lrv,var=='Foraminifera',select=c('stdno','year','day','lon','lat','month')),'Foram',1995,TRUE)
ll[[43]]<-modfac(subset(lrv,var=='Chaetognath',select=c('stdno','year','day','lon','lat','month')),'Chaeto',1995,TRUE)
ll[[44]]<-modfac(subset(lrv,var=='Annelids',select=c('stdno','year','day','lon','lat','month')),'Lrv Annelid',1995,TRUE)
ll[[45]]<-modfac(subset(lrv,var=='Mollusc',select=c('stdno','year','day','lon','lat','month')),'Lrv Moll',1995,TRUE)
ll[[46]]<-modfac(subset(lrv,var=='Cnidarian',select=c('stdno','year','day','lon','lat','month')),'Lrv Cnid',1995,TRUE)
ll[[47]]<-modfac(subset(lrv,var=='Ctenophore',select=c('stdno','year','day','lon','lat','month')),'Lrv Cteno',1995,TRUE)
ll[[48]]<-modfac(subset(lrv,var=='Gadids',select=c('stdno','year','day','lon','lat','month')),'Lrv Gadids',1990,TRUE)
ll[[49]]<-modfac(subset(lrv,var=='Dragonfish',select=c('stdno','year','day','lon','lat','month')),'Lrv Dragon',1985,TRUE)
ll[[50]]<-modfac(subset(lrv,var=='Myctophids',select=c('stdno','year','day','lon','lat','month')),'Lrv Mycto',1995,TRUE)
ll[[51]]<-modfac(subset(lrv,var=='Flatfish',select=c('stdno','year','day','lon','lat','month')),'Lrv Flatf',1985,TRUE)
ll[[52]]<-modfac2(subset(cur,select=c('GS.Dist','year','GS.Dist.sd')),'GS Dist',1985,'GS Distance',TRUE)
ll[[53]]<-modfac2(subset(cur,select=c('SS.Dist','year','SS.Dist.sd')),'SS Dist',1985,'SS Distance',TRUE)
ll[[54]]<-modfac2(subset(her,select=c('r','year')),'Her Recruitment',1985,'Her Recruitment',TRUE)
ll[[55]]<-modfac2(subset(her,select=c('s','year')),'Her Survivorship',1985,'Her Survivorship',TRUE)
ll[[56]]<-modfac2(subset(her,select=c('f','year')),'Her F',1985,'Her F',TRUE)
ll[[57]]<-modfac(subset(plank,select=c('richness.z','year','day','lon','lat','month')),'Zoop Rich',1965,TRUE)
ll[[58]]<-modfac2(subset(her,select=c('obs.w','year','n.obs.inrange')),'Her CF obs',1985,'Her CF obs',TRUE)
ll[[59]]<-modfac2(subset(her,select=c('ssb','year')),'Her SSB',1985,'Her SSB',TRUE)
ll[[60]]<-modfac2(subset(her,select=c('wgt','year')),'Her Weight',1985,'Her Weight',TRUE)
ll[[61]]<-modfac2(subset(nao,select=c('nao','year','nao.sd')),'NAO',1985,'NAO',TRUE)
ll[[62]]<-modfac(subset(lrich,select=c('richness.lrv','year','day','lon','lat','month')),'Lrv Fish Rich',1985,TRUE)
ll[[63]]<-modfac(subset(lrich,select=c('shannon.lrv','year','day','lon','lat','month')),'Lrv Fish Div',1985,TRUE)
ll[[64]]<-modfac(subset(lrich,select=c('pe.lrv','year','day','lon','lat','month')),'Lrv Fish Even',1985,TRUE)
ll[[65]]<-modfac(subset(plank,select=c('shannon.z','year','day','lon','lat','month')),'Zoop Div',1995,TRUE)






library(akima)
library(corrplot)
pdat.fac<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), ll)#COMBINE
#names(pdat.fac)<-gsub('p\\.','',names(pdat.fac))
rownames(pdat.fac)<-pdat.fac$year
#pdat1<-pdat1[,-1]
#cr<-cor(pdat.fac,use='pairwise.complete.obs')
#pdat.lin<-pdat.lin[,-1]
names(pdat.fac)<-gsub('\\.',' ',names(pdat.fac))
names(pdat.fac)<-gsub('p  ',' ',names(pdat.fac))
names(pdat.fac)<-gsub('  ','',names(pdat.fac))


setwd(figsdir)
pdatt<-pdat.fac
cr<-cor(pdatt,use='pairwise.complete.obs')
pdf('corrplot.pdf',height=14,width=14)
par(mar=c(4,1,1,1))
corrplot(cr,method='ellipse',type='lower',diag=FALSE,order='AOE',tl.col='black')
dev.off()


#INTERPOLATES FOR VALUES WHERE ONLY 1 OR 2 YEARS MISSING
pdat<-subset(pdat.fac,year>=1978 & year<=2008)
#pdat<-pdat[,!(names(pdat) %in% c('CTZtotal','CTZoop','CTNutrients','phos','sil','nit','Prich1phyto','CFpred'))]

lnfun<-function(ddd){
names(ddd)[2]<-'y'
y<-unique(subset(ddd,is.na(y)==FALSE,select=c('year')))
y$dm<-rep(1,dim(y)[1])
years<-data.frame(year=seq(1980,2006,1))
y<-merge(years,y,by=c('year'),all.x=TRUE,all.y=FALSE)
y$dm<-ifelse(is.na(y$dm)==TRUE,-9,y$dm)
y <- rle(y$dm)
options(warn=-1)
y2<-max(y$lengths[y$values==-9])#GETS LONGEST STRETCH OF MISSING DATA
options(warn=0)
if(y2 %in% c(-Inf,Inf)){return(0)
}else {return(y2)}
}

vbls<-names(pdat[2:dim(pdat)[2]])
lt<-list()
setwd(figsdir)
pdf('interp_figs.pdf',height=10,width=10)
par(mfrow=c(7,5),mar=c(1.5,1,1,1),oma=c(1,1,1,1))
for(i in 1:length(vbls)){
    d<-subset(pdat,select=c('year',paste(vbls[i])))
    names(d)[2]<-'y'
    miss<-lnfun(na.omit(d))

    pd<-data.frame(year=seq(min(d$year),max(d$year),1))
    pd<-merge(pd,d,by=c('year'),all.x=TRUE,all.y=FALSE)
    pd1<-pd
    if(miss<=2){
        pd$p<-aspline(x=d$year,y=d$y,xout=pd$year)$y
        pd$mn<-mean(pd$p)
        pd$sdd<-sd(pd$p)
        pd$z<-(pd$p-pd$mn)/pd$sd
        pd<-subset(pd,year>1980 & year<=2005)
        pd1<-subset(pd1,year>1980 & year<=2005)
        pd$z[1]<-ifelse(is.na(pd$y)[1]==TRUE,NA,pd$z[1])
        pd$p[1]<-ifelse(is.na(pd$y)[1]==TRUE,NA,pd$p[1])
        ylb<-c(min(pd$p,na.rm=TRUE),max(pd$p,na.rm=TRUE))
        plot(pd$year,pd$p,pch=15,col='red',xlab='',ylab='',las=1,xlim=c(1980,2005),axes=FALSE)
        axis(1,at=c(1980,2005))
        axis(2,at=ylb,labels=FALSE,tick=TRUE)
        legend('topleft',vbls[i],bty='n')
        points(pd1$year,pd1$y,pch=15)
        pd<-subset(pd,select=c('year','z'))
    }    else {
        pd<-data.frame(year=2000,p=NA)
    }
    names(pd)[2]<-vbls[i]
    lt[[i]]<-pd
}
dev.off()
pdat.int<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), lt)#COMBINE
pdat.int<-pdat.int[,colSums(is.na(pdat.int)) != nrow(pdat.int)]#REMOVES COLUMNS THAT ARE ALL MISSING
#pdat.int<-pdat.int[,-1]
setwd(datadir)
write.csv(pdat.int,'mm_data.csv',row.names=FALSE)
write.csv(pdat.fac,'pdat.fac.csv',row.names=FALSE)



setwd(figsdir)
pdf('linear_univariate_models.pdf',height=10,width=10)
par(mfrow=c(7,5),mar=c(1,.5,1,1),oma=c(1,1,1,1))

pp<-pdat.fac[,names(pdat.fac) %in% c('Her CF obs','Her Recruitment','Her SSB','Her Survivorship','Her Weight','Herr Lrv')]
resp<-names(pp)
pp<-pdat.fac[,!(names(pdat.fac) %in% c('Her CF obs','Her Recruitment','Her SSB','Her Survivorship','Her Weight','Herr Lrv','year'))]
preds<-names(pp)


ii<-list()
for(i in 1:length(resp)){
    print(resp[i])
    jj<-list()
    for(j in 1:length(preds)){
        d<-subset(pdat.fac,select=c(resp[i],preds[j]))
        l<-dim(na.omit(d))[1]
        names(d)<-c('y','x')
        mod<-lm(y~x,data=d)
        s<-summary(mod)
        ttab<-data.frame(resp=resp[i],
                         pred=preds[j],
                         est=round(s$coef[2,1],digits=2),
                         estse=round(s$coef[2,2],digits=2),
                         pv=round(s$coef[2,4],digits=3))
ttab$sig<-ifelse(ttab$pv<=0.0001,'***','-')
ttab$sig<-ifelse(ttab$pv>0.0001 & ttab$pv<=0.001,'**',ttab$sig)
ttab$sig<-ifelse(ttab$pv>0.001 & ttab$pv<=0.05,'*',ttab$sig)
ttab$r2<-round(s$r.squared,digits=2)
ttab$n<-l

        xx<-seq(min(d$x,na.rm=TRUE),max(d$x,na.rm=TRUE),length.out=50)
        plot(d$x,d$y,pch=15,xlab=preds[j],ylab=resp[i],las=1,axes=FALSE)
        xlb<-c(min(d$x,na.rm=TRUE),max(d$x,na.rm=TRUE))
        ylb<-c(min(d$y,na.rm=TRUE),max(d$y,na.rm=TRUE))
        axis(1,at=xlb,labels=FALSE,tick=TRUE)
        axis(2,at=ylb,labels=FALSE,tick=TRUE)
        p<-predict(mod,newdata=data.frame(x=xx))
        lines(xx,p)
#        legend('topright',paste('r2 =',ttab$r2),bty='n')
        legend('topleft',paste(resp[i],'=',preds[j]),bty='n',col='red3')
        jj[[j]]<-ttab
    }
ii[[i]]<-data.frame(do.call('rbind',jj))
}
lmout<-data.frame(do.call('rbind.fill',ii))
dev.off()

lmout<-subset(lmout,is.na(pv)==FALSE & n>=10 & pv<=.05)
f<-function(d){    d<-d[order(d$r2,decreasing=TRUE),]}
lmout<-ddply(lmout,.(resp),.fun=f)
write.csv(lmout,'lin_uni_models.csv',row.names=FALSE)




#MULTIVARIATE MODEL SELECTION VIA BACKWARD AIC
pp<-pdat.int[,names(pdat.int) %in% c('Her CF obs','Her Recruitment','Her SSB','Her Survivorship','Her Weight','Herring Lrv')]
resp<-names(pp)
pp<-pdat.int[,!(names(pdat.int) %in% c('Her CF obs','Her Recruitment','Her SSB','Her Survivorship','Her Weight','Herr Lrv'))]
preds<-names(pp)


library(fields)
library(RColorBrewer)
library(heplots)
ii<-list()
for(i in 1:length(resp)){
print(resp[i])
d<-cbind(subset(pdat.int,select=c(resp[i])),pdat.int[,!(names(pdat.int) %in% resp)])
d<-d[,!(names(d) %in% c('year','Phyt Rich CT','Wind Seas max','Wind Dur','SST Amp','SST Dur'))]
names(d)[1]<-'y'
d<-subset(d,is.na(y)==FALSE)
d <- d[,colSums(is.na(d)) <=10]#EXCLUDE TS WHERE >20 MISSING
print(dim(d))
#        mod<-lm(y~.*.,data=d)       mod<-lm(y~.^2,data=d)
mod<-stepAIC(lm(y~.,data=d),direction='backward')
et<-data.frame(etasq(mod))
s<-summary(mod)
df<-data.frame(s$coef)
ttab<-data.frame(resp=resp[i],
                 pred=rownames(df),
                         est=round(df$Estimate,digits=2),
                         estse=round(df$Std..Error,digits=2),
                         pv=round(df$Pr...t..,digits=3))
ttab<-ttab[-1,]
ttab<-ttab[order(ttab$pv),]
ttab$sig<-ifelse(ttab$pv<=0.0001,'***','-')
ttab$sig<-ifelse(ttab$pv>0.0001 & ttab$pv<=0.001,'**',ttab$sig)
ttab$sig<-ifelse(ttab$pv>0.001 & ttab$pv<=0.05,'*',ttab$sig)
ttab$pr2<-round(data.frame(etasq(mod))[1:dim(ttab)[1],],digits=2)
ttab$r2<-round(s$r.squared,digits=2)

ii[[i]]<-ttab
}
mvmout<-data.frame(do.call('rbind.fill',ii))

write.csv(mvmout,'lin_multi_models.csv',row.names=FALSE)




################################################################3
#REDUNDANCY ANALYSIS OPERATES ON INTERPOLATED DATA - NEEDS NON-MISSING DATA
pdat.int<-subset(pdat.int,year>=1982)
prds<-pdat.int[,names(pdat.int) %in% preds]
rsp<-pdat.int[,names(pdat.int) %in% resp]
rownames(rsp)<-prds$year
rownames(prds)<-prds$year

names(prds)<-gsub(' ','_',names(prds))
srda<-rda(rsp~SST+ Wind+Strat+ Temp_50m+Sil+ Phos+Phyt_Div+  Phyt_Cells+Phyt_Rich+ Phyt_Even+ SST_Timing+  SST_Amp+ SST_Dur+ SST_Seas_max+Wind_Timing+ Wind_Amp+Wind_Dur+Wind_Seas_max +Phyt_TO+ Chl_CT+Phyt_Div_CT+ Phyt_Even_CT+Phyt_Rich_CT  +Phyt_Cells_CT+ GS_Dist_CT+SS_Dist_CT+GS_Dist+ SS_Dist+ Her_F+ NAO,data=prds,scale=TRUE)


RsquareAdj(srda)$r.squared
#biplot(srda,scaling='symmetric',type=c('text','points'))
#ordiellipse(srda,)

setwd(figsdir)
pdf('rdaplot.pdf',height=8,width=9)
xlm <-c(-1.6,1.5)
ylm<-c(-1.25,1.5)
clenv<-'gray10'
clresp<-'springgreen4'
yr<-ifelse(as.numeric(rownames(rsp))<=1988,'red3','green')
cls<-c(rev(heat.colors(dim(rsp)[1])),'darkred')

plot(0,0,xlim=xlm,ylim=ylm,las=1,xlab='RDA1',ylab='RDA2',xaxt='n',yaxt='n',bg='gray')
rect(xlm[1]-1,ylm[1]-1,xlm[2]+1,ylm[2]+1,col='gray80')
points(srda,display=c('sites'),scaling=2,xlim=xlm,ylim=ylm,cex=2,col=cls,pch=16)
points(srda,display=c('sites'),scaling=2,xlim=xlm,ylim=ylm,cex=2,col='black',pch=1,lwd=.1)
#plot(srda,dislay='sp',add=TRUE)
spe.sc <- scores(srda, choices=1:2, scaling=2, display="sp",xlim=xlm,ylim=ylm)
env.sc <- scores(srda, choices=1:2, scaling=2, display="sites",xlim=xlm,ylim=ylm)
#sc <- scores(srda, choices=1:2, scaling=2, display="bp",xlim=xlm,ylim=ylm)
#arrows(0, 0, sc[, 1], sc[, 2], length=.2, lty=1, col="black",xlim=xlm,ylim=ylm,lwd=.5)
text(srda,display=c('bp'),scaling=2,xlim=xlm,ylim=ylm,cex=.75,col=alpha(clenv,.5))
text(srda,display=c('sp'),scaling=2,xlim=xlm,ylim=ylm,cex=1,col=clresp)
axis(side=1,at=round(seq(xlm[1],xlm[2],.4),digits=1),tick=T,labels=T)
axis(side=2,at=seq(ylm[1],ylm[2],.4),tick=T,labels=T,las=1)
box(lwd=1)
#arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length=.2, lty=1, col="gray40",xlim=xlm,ylim=ylm,lwd=.5)
#arrows(0, 0, env.sc[, 1], env.sc[, 2], length=.2, lty=1, col="red",xlim=xlm,ylim=ylm,lwd=.5)
#text(srda,display=c('sites'),scaling=2,xlim=xlm,ylim=ylm,cex=.6,col=yr)
zr<-sort(as.numeric(rownames(rsp)))
colorbar.plot(1,1.15,zr,col=cls,horizontal=TRUE,strip.width=.04,strip.length=.3)
dev.off()






setwd(figsdir)
pdf('mda_plots.pdf',height=12,width=14)
par(mfrow=c(2,2),mar=c(4,4,1,1))
mdfun.cor<-function(yr1,yr2){
df<-subset(pdat.fac,year>=yr1 & year<=yr2)
df<-df[,!(names(df) %in% c('year','Phyt Rich','Phyt Rich Ct','Zoop Rich','Lrv Fish Rich Ct','Lrv Fish Rich'))]
n<-round(dim(df)[1]*0.5,digits=0)
df <- df[,colSums(is.na(df)) <=n]#EXCLUDE TS WHERE >20 MISSING
#colSums(is.na(df))
cldb<-data.frame(var=names(df))
cldb$cl<-ifelse(regexpr('Lrv', cldb$var)>0,'blue3','green3')
cldb$cl<-ifelse(regexpr('Her', cldb$var)>0,'black',cldb$cl)
cldb$cl<-ifelse(regexpr('Zoop', cldb$var)>0,'gold3',cldb$cl)
cldb$cl<-ifelse(regexpr('Dist', cldb$var)>0 |
                regexpr('Nitrate', cldb$var)>0 |
                regexpr('Wind', cldb$var)>0 |
                regexpr('Temp', cldb$var)>0 |
                regexpr('Strat', cldb$var)>0 |
                regexpr('SST', cldb$var)>0 |
                regexpr('Silicate', cldb$var)>0 |
                regexpr('Phosphate', cldb$var)>0 |
                regexpr('NAO', cldb$var)>0
               ,'red3',cldb$cl)
cldb$cl<-as.character(cldb$cl)

cldb$lab<-ifelse(regexpr('Lrv', cldb$var)>0,'Larval Fish','Phytoplankton')
cldb$lab<-ifelse(regexpr('Her', cldb$var)>0,'Herring',cldb$lab)
cldb$lab<-ifelse(regexpr('Zoop', cldb$var)>0,'Zooplankton',cldb$lab)
cldb$lab<-ifelse(regexpr('Dist', cldb$var)>0 |
                regexpr('Nitrate', cldb$var)>0 |
                regexpr('Wind', cldb$var)>0 |
                regexpr('Temp', cldb$var)>0 |
                regexpr('Strat', cldb$var)>0 |
                regexpr('SST', cldb$var)>0 |
                regexpr('Silicate', cldb$var)>0 |
                regexpr('Phosphate', cldb$var)>0 |
                regexpr('NAO', cldb$var)>0
               ,'Environment',cldb$lab)

#dmat<-1-(cor(df,use='pairwise.complete.obs'))
dmat<-1-(abs(cor(df,use='pairwise.complete.obs')))+.001
#dmat<-abs(cor(df,use='pairwise.complete.obs'))
dst<-as.dist(dmat)
#dst<-dist(df,method='euclidean')
ord<-metaMDS(dst,trymax=50,method='euclidean',trace=0)
plot(0,0,col='white',xlim=c(-.5,.6),ylim=c(-.5,.5),las=1)
text(ord, display = "sites", cex = 0.8, pch=21, col=cldb$cl, bg="yellow",adj=-.2)
points(ord, display = "sites", cex =1, col=alpha(cldb$cl,.5),pch=16)
points(ord, display = "sites", cex =1,,pch=1,lwd=.1)
cldb<-unique(subset(cldb,select=c('cl','lab')))
legend('topleft',cldb$lab,col=as.character(cldb$cl),pch=16,bty='n',pt.cex=2)
legend('bottomleft',paste(yr1,'-',yr2),bty='n')
}
par(mfrow=c(2,2))
mdfun.cor(1980,2006)
mdfun.cor(1965,2015)
mdfun.cor(1990,2015)


mdfun.dst<-function(yr1,yr2){
df<-subset(pdat.int,year>=yr1 & year<=yr2)
rownames(df)<-df$year
df<-df[,!(names(df) %in% c('year','Phyt Rich','Phyt Rich Ct','Zoop Rich','Lrv Fish Rich Ct','Lrv Fish Rich'))]

dst<-vegdist(df,method='euclidean')
ord<-metaMDS(dst,trymax=50,trace=0)
#ord<-metaMDS(decostand(df,method='standardize'),trace=0,distance='euclidean')
cl<-rev(heat.colors(dim(df)[1]))
ordiplot(ord,bg='gray',las=1)
rect(-8,-7,7,7,col='gray85')
clps<-ifelse(as.numeric(rownames(df))<=1990,'pre-collapse','post-collapse')
ordihull(ord,groups=clps,draw='lines',label=FALSE,col='black')
text(ord, display = "sites", cex = 1, pch=21, col=as.character(cl),adj=-.2)
points(ord, display = "sites", cex=1.5, col=alpha(cl,.5),pch=16)
points(ord, display = "sites", cex=1.5,,pch=1,lwd=.1)
zr<-as.numeric(rownames(df))
colorbar.plot(4,3.5,zr,col=as.character(cl),horizontal=TRUE,strip.width=.02,strip.length=.2)
legend('bottomleft',paste(yr1,'-',yr2),bty='n')
}
par(mfrow=c(2,2))
mdfun.dst(1982,2005)
mdfun.dst(1990,2015)
dev.off()





df<-subset(lrv,var=='Herring')
f<-function(d){
    names(d)[1]<-'time'
    print(dim(d))
    if(dim(d)[1]>100){
    mod<-gam(stdno~s(lon,lat,k=75),data=d)
    clonep<-seq(min(df$lon),max(df$lon),length.out=500)
    clatep<-seq(min(df$lat),max(df$lat),length.out=500)
    #PREDICTION DATA
    pdat2<-expand.grid(lon=clonep,lat=clatep)
    pdat2$p<-predict(mod,newdata=pdat2,type='response')
    dout<-subset(pdat2,p==max(pdat2$p))[1,]
    dout$n<-dim(d)[1]
    return(dout)
    } else NULL
}
com<-ddply(subset(df,select=c('tbin3','stdno','lon','lat')),.(tbin3),.fun=f,.progress='text')
comy<-ddply(subset(df,select=c('year','stdno','lon','lat')),.(year),.fun=f,.progress='text')

map('world',xlim=c(-70,-64),ylim=c(42,45))
points(comy$lon,comy$lat)
points(com$lon,com$lat)
plot(as.numeric(as.character(com$tbin3)),com$lat)
plot(as.numeric(as.character(com$tbin3)),com$lon)
plot(comy$year,comy$lat)
plot(as.numeric(as.character(com$tbin3)),com$lon)


f<-function(d){
    names(d)[1]<-'time'
    print(dim(d))
    if(dim(d)[1]>100){
    mod<-gam(stdno~s(lon,lat,k=75),data=d)
    clonep<-seq(min(df$lon),max(df$lon),length.out=500)
    clatep<-seq(min(df$lat),max(df$lat),length.out=500)
    #PREDICTION DATA
    pdat2<-expand.grid(lon=clonep,lat=clatep)
    pdat2$p<-predict(mod,newdata=pdat2,type='response')
    dout<-subset(pdat2,p==max(pdat2$p))[1,]
    dout$n<-dim(d)[1]
    return(dout)
    } else NULL
}
cos<-ddply(subset(df,select=c('tbin3','stdno','lon','lat')),.(tbin3),.fun=f,.progress='text')
cosy<-ddply(subset(df,select=c('year','stdno','lon','lat')),.(year),.fun=f,.progress='text')



library(reshape2)
library(zoo)
dat.m <- data.frame(yr=index(dat),value=melt(dat)$value)


dm<-na.omit(subset(her,select=c('year','s','ssb')))
s<-subset(dm,select=c('year','s'))
ssb<-subset(dm,select=c('year','ssb'))
ssb$year<-ssb$year+
s<-as.ts(s$s)
ssb<-as.ts(ssb$ssb)
par(mfrow=c(2,2))
plot(s$year,s$s,type='l')
plot(ss)
lines(lag(ss,k=-5),col='red')
plot(ssb)
lines(lag(ssb,k=-5),col='red')

s<-subset(dm,select=c('year','s'))
sts<-zoo(s$s,order.by=s$year,frequency=1)
plot(sts)
lines(lag(sts,k=5),col='red')

s<-subset(dm,select=c('year','s'))
sts<-zoo(s$s,order.by=s$year,frequency=1)
lt<-list()
for(i in 1:5){
slag<-lag(sts,k=i)
df<-data.frame(year=index(slag),
               surv=as.matrix(slag))
names(df)[2]<-gsub(' ','',paste(names(df)[2],i))
lt[[i]]<-df
}
lgdat<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), lt)#COMBINE

ssb<-subset(dm,select=c('year','ssb'))
lgdat<-merge(lgdat,ssb,by=c('year'),all.x=TRUE,all.y=FALSE)
lgdat<-subset(lgdat,year>=1965)
round(cor(lgdat,use='pairwise.complete.obs'),digits=2)


ccf(dm$s,dm$ssb,ylim=c(-.5,.5))
abline(v=0,lty=2)
plot(dm$s,dm$ssb)

dst<-vegdist(decostand(pdat.int,method='standardize'),method='euclidean')
ord<-metaMDS(pdat.int,trymax=50,distance='euclidean')
#ord<-isoMDS(dst)
plot(ord)
plot(0,0,col='white',xlim=c(-.5,.5),ylim=c(-.5,.5),las=1)
text(ord, display = "sites", cex = 0.8, pch=21, col=cldb$cl, bg="yellow",adj=0)
#text(ord, display = "sites", cex = 0.8, pch=21, col="blue", bg="yellow",adj=0)
points(ord, display = "sites", cex = 0.8, col=cldb$cl, bg="yellow",pch=15)
text(ord, display = "spec", cex=0.7, col="blue")

stressplot(ord,dst)
decostand(varespec, "standardize")

ord<-isoMDS(dst)
autoplot(ord,label=TRUE)
autoplot(ord,label=TRUE,colour=cldb$cl)
slm<-sammon(dst)
autoplot(slm,label=TRUE)


dd<-subset(pdat.fac,select=c('HerSSB','SST','Wind','Stratification'))
mod<-lm(HerSSB~.^3,data=dd)
mod<-lm(HerSSB~.*.,data=dd)

library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)
library(ggfortify)
pd<-pdat.int
pd<-pd[,!(names(pd) %in% c('Rich1phyto','NAO','SSTduration.10','SSTampitude','Windampitude','CTPtotal','CTSST','PrichnessCT','PPtotalCT','Windduration.10'))]
#era<-ifelse(pd$year<=1987,'pre','collapse')
#era<-ifelse(pd$year>1995,'post',era)
era<-ifelse(pd$year<=1988,'pre','post')
era<-as.factor(era)
pd<-pd[,-1]

lbls<-sort(unique(era))
#c('collapse','pose','pre')
cls<-colorRampPalette(c('green3','magenta4','royalblue'))
cls<-cls(length(lbls))
xlm<-c(-5,8)
ylm<-c(-5,4)
pc<-prcomp(pd,center=TRUE,scale=TRUE)
s<-summary(pc)
plot(seq(1980,2006,1),pc$x[,1],type='l')
plot(seq(1980,2006,1),pc$x[,2],type='l')
plot(seq(1980,2006,1),pc$x[,3],type='l')

setwd(figsdir)
pdf('pcaplot.pdf',height=7,width=7)
ggbiplot(pc, obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = FALSE,groups=era, var.axes=TRUE,size=5,shape=15)+
    scale_colour_manual(values=as.character(cls),breaks=lbls,guide=guide_legend(title=paste('Era')))+
    theme(legend.direction = 'horizontal', legend.position = 'top')+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='top',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.2, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=xlm,labels=xlm,limits=xlm)+
scale_y_continuous(expand=c(0,0),breaks=ylm,labels=ylm,limits=ylm)+
    coord_equal()+
    coord_cartesian(ylim=ylm,xlim=xlm)
dev.off()


names(pd)<-gsub('>','',names(pd))
rownames(pd)<-seq(1980,2006,1)

pd2<-pd[,!(names(pd) %in% c('HerSSB','HerWeight','HerRecruitment','HerSurvivorship','HerF'))]
pd2<-pd[,(names(pd) %in% c('SST','WInd','Stratification','Temperature50m'))]
pd2<-pd[,(names(pd) %in% c('ShannonP','TotPhyto','RIchphyto','Pielouphyto','PhytoTO','Richphyto'))]

pd2<-pdat.fac[,!((names(pdat.fac) %in% c('HerSSB','HerWeight','HerRecruitment','HerSurvivorship','HerF')))]

nms<-c('HerSSB','HerWeight','HerRecruitment','HerSurvivorship','HerF')
pd2<-(pdat.fac[,names(pdat.fac) %in% c('HerSSB','HerWeight','HerRecruitment','HerSurvivorship','HerF','LHerring')])[5:54,]



p1<-autoplot(kmeans(pd2,2),data=pd,label=TRUE,label.size=5,,alpha=.1,frame=TRUE)
p2<-autoplot(pam(pd2, 3), frame = TRUE, frame.type = 'norm',label=TRUE)
grid.arrange(p2,p3,ncol=2)
ag<-(agnes(pdat.fac,diss=FALSE,metric='manhattan',stand=TRUE))



#PLOTS ORDIPLOTS
ordfun<-function(pd2,n){
p1<-autoplot(kmeans(pd2,n),data=pd2,label=TRUE,label.size=3,,alpha=.1,frame=TRUE)
p2<-autoplot(pam(pd2, n), frame = TRUE, frame.type = 'norm',label=TRUE)
p3<-autoplot(clara(pd2,n),label=TRUE,alpha=.1,label.size=3,frame=TRUE)
pc2<-prcomp(pd2,center=TRUE,scale=TRUE)
p4<-autoplot(pc2,label=TRUE,alpha=.1,label.size=3,frame=TRUE)

dmat<-1-cor(pd2,use='pairwise.complete.obs')
dst<-as.dist(dmat)
grid.arrange(p1,ncol=1)
grid.arrange(p2,ncol=1)
grid.arrange(p3,ncol=1)
grid.arrange(p4,ncol=1)
}



setwd(figsdir)
pdf('ordiplots_herring.pdf',width=10,height=8)
par(mar=c(8,4,8,1))
psidat<-unique(subset(bpdat,var %in% c('Her F','Her SSB','Her CF obs','Her Weight','L Herring')))
her$ssb2<-her$ssb/1000
plot(her$year,her$ssb2,type='l',las=1,xlab='Year',ylab='SSB [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
mx<-subset(her,ssb2==max(her$ssb2,na.rm=TRUE))
mn<-subset(her,ssb2==min(her$ssb2,na.rm=TRUE))
lines(c(mx$year-1,mx$year),c(mx$ssb2,mx$ssb2))
text(mx$year-1.75,mx$ssb2,round(mx$ssb2,digits=0))
text(mn$year,mn$ssb2+40,round(mn$ssb2,digits=0))

plot(her$year,her$ssb2,type='l',las=1,xlab='Year',ylab='SSB [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
cl1<-'gray80'
cl2<-'gray50'
cl3<-'gray10'
mx<-max(her$ssb2,na.rm=TRUE)+100
rect(1960,-100,1988.5,mx,col=alpha(cl1,.75),border=NA)
rect(1988.5,-100,1995,mx,col=alpha(cl2,.75),border=NA)
rect(1995,-100,2010,mx,col=alpha(cl3,.75),border=NA)
cl<-'red3'
rug(psidat$bpt,lwd=2,col=cl)
text(psidat$bpt[1],0,'Herring Condition Factor declines',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[2],0,'F increase',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[3]-.5,0,'SSB collapse',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[4],0,'Herring weights decline',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[5]+.5,0,'Herring larvae decline',srt=40,adj=0,col=cl,cex=.75)
text(1975,750,'Pre-collapse')
text(1991,750,'Collapse')
text(2001,750,'Post-collapse')
mx<-subset(her,ssb2==max(her$ssb2,na.rm=TRUE))
mn<-subset(her,ssb2==min(her$ssb2,na.rm=TRUE))
lines(c(mx$year-1,mx$year),c(mx$ssb2,mx$ssb2))
lines(her$year,her$ssb2,col=cl,lwd=2)
text(mx$year-1.75,mx$ssb2,round(mx$ssb2,digits=0))
text(mn$year,mn$ssb2+40,round(mn$ssb2,digits=0))

her$r2<-her$r/10000
plot(her$year,her$r2,type='l',las=1,xlab='Year',ylab='R [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
mx<-subset(her,r2==max(her$r2,na.rm=TRUE))
mn<-subset(her,year==2004)
lines(c(mx$year-1,mx$year),c(mx$r2,mx$r2))
text(mx$year-1.75,mx$r2,round(mx$r2,digits=0))
text(mn$year+1,mn$r2+10,round(mn$r2,digits=0))
lines(her$year,her$ssb2,type='l',col='red')

par(mar=c(6,4,4,2))
pd2<-pdat.fac[,names(pdat.fac) %in% c('Her SSB','Her Weight','Her Recruitment','Her Survivorship','Her F')][5:54,]
names(pd2)<-gsub(' ','_',names(pd2))
ag<-(agnes(pd2,diss=FALSE,metric='euclidean',stand=TRUE))
agh<-as.hclust(ag)
labelColors = c('red3','royalblue','forestgreen','gold3')
clusMember = cutree(agh, 4)
# function to get color labels
colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}
# using dendrapply
clusDendro = dendrapply(as.dendrogram(agh), colLab)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7)

dn<-(diana(pd2,diss=FALSE,metric='euclidean',stand=TRUE))
dnh<-as.hclust(dn)
labelColors = c('red3','royalblue','forestgreen','gold3')
clusMember = cutree(dnh, 2)
clusDendro = dendrapply(as.dendrogram(dnh), colLab)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7,type='triangle')

par(mar=c(8,4,8,2))
plot(clusDendro, main = "",las=1,horiz=FALSE,cex=.7)

#par(mfrow=c(2,2),mar=c(4,4,1,1))
ordfun(na.omit(pd2),3)
dev.off()



ordfun(c('ShannonP','TotPhyto','RIchphyto','Pielouphyto','PhytoTO','Richphyto'))
ordfun(c('SST','WInd','Stratification','Temperature50m','Silicate','Phosphate','GSDistance','SSDistance'))
ordfun(c('SSTtiming','SSTduration10','SSTseasonalmax','Windtiming','Windduration10','Windseasonalmax','ChlCT','PdiversityCT','PevennessCT','GSDistCT','SSDistCD'))
ordfun(names(pd))
dev.off()


autoplot(kmeans(pd,2),data=pd,label=TRUE,label.size=5,,alpha=.1,frame=TRUE)
autoplot(pam(pd, 2), frame = TRUE, frame.type = 'norm')
autoplot(clara(pd,3),label=TRUE,alpha=.1,label.size=5,frame=TRUE)

autoplot(prcomp(pd2),loadings=TRUE,colour='era',data=pd,loadings.colour='blue',loadings.label=TRUE,size=3,shape=15)


df<-pd
k<-3
dmat<-1-cor(df,use='pairwise.complete.obs')
dst<-as.dist(dmat)
ff<-fanny(dst,k,maxit=5000,diss=T)

dum<-c('red3','magenta3','forestgreen','gold','cornflowerblue','darkblue')
plot(silhouette(ff),col=dum[1:k],main='')#silhouette plot

dc.pcoa<-cmdscale(dst)
dc.scores<-scores(dc.pcoa,choices=c(1,2))
spefuz.g<-ff$clustering
a<-data.frame(var=as.character(sort(unique(names(df)))),
              clusters=ff$clustering)
aa<-data.frame(ff$membership)
aa$var<-rownames(a)

par(mar=c(1,1,1,8),oma=c(1,1,1,1))
plot(scores(dc.pcoa),asp=1,type='n',xlim=c(-.75,1.1),ylim=c(-1,1.2),las=1,axes=FALSE,xlab='',ylab='')
stars(ff$membership,location=scores(dc.pcoa),draw.segments=T,add=T,scale=F,len=.1,col.segments=alpha(c(dum[1:k]),.75),byt='n',labels=NULL,xlim=c(-1.1,1.4),ylim=c(-1,.5),lwd=.0001,xpd=TRUE,border=NULL,radius=FALSE,col.radius=alpha('white',.1))
for(i in 1:k){ cl<-dum[i]
    gg<-dc.scores[spefuz.g==i,]
    hpts<-chull(gg)
    hpts<-c(hpts,hpts[1])
    lines(gg[hpts,],col=cl,lwd=3,xlim=c(-1.1,1.4),ylim=c(-1,.5))
}


cx <- data.frame(cell=aa$var,
                 cx=apply(aa[, 1:k], 1, max))
cx$cxx <- rescale(cx$cx,newrange=c(.05,.5))










autoplot(fanny(pd,k,memb.exp=1.01,maxit=5000,diss=F),frame=TRUE)

autoplot(ff)




vbls<-names(pdat1)
m<-matrix(0,dim(pdat1)[1],dim(pdat1)[2])
for(i in 1:length(vbls)){
    d<-subset(pdat1,select=c(paste(vbls[i])))
    names(d)[1]<-'y'
    d$mn<-mean(d$y,na.rm=TRUE)
    d$sdd<-sd(d$y,na.rm=TRUE)
    d$z<-(d$y-d$mn)/d$sdd
    m[,i]<-d$z
}





detach("package:raster", unload=TRUE)
library(RColorBrewer)








pc<-prcomp(pdat1,center=TRUE,scale=TRUE,na.action=na.omit)


pdat1<-rbind.fill(l)



round(cor(subset(strat,select=c('s25','s50','temp25','temp50','sal25','sal50','temp1')),use='pairwise.complete.obs'),digits=2)

zz<-dlply(data,.(db),.fun=modf)








modf<-function(d){
    mod<-gam(log10(chl+.1)~as.factor(year) + s(day,bs='cc',k=5) + s(lon,lat,k=20),data=d,gamma=1.4)
    lon<-mean(d$lon)
    lat<-mean(d$lat)
    day<-300
    pdat<-data.frame(year=sort(unique(d$year)))
    p<-predict(mod,newdata=data.frame(year=pdat$year,lon=lon,lat=lat,day=day),se.fit=TRUE)
    pdat$p<-(10^p$fit)-.1
    pdat$se<-(10^p$se.fit)-.1
    pdat$upr<-pdat$p+(1.96*pdat$se)
    pdat$lwr<-pdat$p-(1.96*pdat$se)
#    plot(mod,select=1)
    plot(pdat$year,pdat$p,pch=15,main=as.character(unique(d$db)))
    f<-function(dd){
        lines(c(dd$year,dd$year),c(dd$lwr,dd$upr))
    }
    zz<-dlply(pdat,.(year),.fun=f)
}
par(mfrow=c(2,3))
zz<-dlply(data,.(db),.fun=modf)


modf<-function(d){
    d<-subset(d,day>=213)
    mod<-gam(log10(chl+.1)~as.factor(year) + s(day,k=4) + s(lon,lat,k=20),data=d,gamma=1.4)
    mods<-gam(log10(chl+.1)~s(year,k=5) + s(day,k=4) + s(lon,lat,k=20),data=d,gamma=1.4)
    lon<-mean(d$lon)
    lat<-mean(d$lat)
    day<-300
    pdat<-data.frame(year=sort(unique(d$year)),
                     n=tapply(d$chl,d$year,length))
    p<-predict(mod,newdata=data.frame(year=pdat$year,lon=lon,lat=lat,day=day),se.fit=TRUE)
    pdat$p<-(10^p$fit)-.1
    pdat$se<-(10^p$se.fit)-.1
    pdat$upr<-pdat$p+(1.96*pdat$se)
    pdat$lwr<-pdat$p-(1.96*pdat$se)
#    plot(mod,select=1)
    plot(pdat$year,pdat$p,pch=15,main=as.character(unique(d$db)),cex=rescale(pdat$n,newrange=c(.75,4)))
    f<-function(dd){  lines(c(dd$year,dd$year),c(dd$lwr,dd$upr))    }
#    zz<-dlply(pdat,.(year),.fun=f)
    pdat$db<-unique(d$db)
    return(pdat)
    pdat<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
    p<-predict(mods,newdata=data.frame(year=pdat$year,lon=lon,lat=lat,day=day),se.fit=TRUE)
    pdat$p<-(10^p$fit)-.1
    pdat$se<-(10^p$se.fit)-.1
    pdat$upr<-pdat$p+(1.96*pdat$se)
    pdat$lwr<-pdat$p-(1.96*pdat$se)
    lines(pdat$year,pdat$p)
    lines(pdat$year,pdat$upr,lty=2)
    lines(pdat$year,pdat$lwr,lty=2)
}
par(mfrow=c(2,3))
zz<-ddply(data,.(db),.fun=modf)

plot(0,0,xlim=c(1960,2016),ylim=c(0,3))
points(subset(zz,db=='SW')$year,subset(zz,db=='SW')$p,pch=15,col='red')
points(subset(zz,db=='MO')$year,subset(zz,db=='MO')$p,pch=15,col='blue')
points(subset(zz,db=='ME')$year,subset(zz,db=='ME')$p,pch=15,col='green')
points(subset(zz,db=='CZ')$year,subset(zz,db=='CZ')$p,pch=15,col='purple')
points(subset(zz,db=='INSITU')$year,subset(zz,db=='INSITU')$p,pch=15,col='gold')
points(subset(zz,db=='CPR')$year,subset(zz,db=='CPR')$p,pch=15,col='black')

a<-subset(zz,year==2003)
f<-function(d){
    return(data.frame(delt=subset(a,db=='SW')$p-d$p))
}
aa<-ddply(a,.(db),.fun=f)

zz<-merge(zz,aa,by=c('db'),all.x=TRUE,all.y=FALSE)
zz$delt<-ifelse(is.na(zz$delt)==TRUE,0,zz$delt)
zz$p2<-zz$p+zz$delt

zz$scl<-rescale(zz$n,newrange=c(1,3))
plot(0,0,xlim=c(1960,2016),ylim=c(.8,3),las=1)
nms<-c('INSITU','CPR','SW','MO','ME','CZ')
cls<-c('red3','blue','green','purple','gold','orange')
for(i in 1:length(nms)){
    d<-subset(zz,db==nms[i])
    points(d$year,d$p2,pch=16,col=alpha(cls[i],.5),cex=d$scl)
}
legend('topleft',nms,col=cls,pch=16,bty='n')


f<-function(d){
map('world',xlim=c(-70,-63),ylim=c(43,45))
crd<-unique(subset(d,select=c('lon','lat')))
points(crd$lon,crd$lat,pch=16,col='red')
}
par(mfrow=c(2,3))
zz<-dlply(data,.(db),.fun=f)

a<-subset(data,month>=8)










######################################################


######################################################
#SUMMARY PLOTS OF HERRING ASSESSMENTS

setwd(figsdir)
pdf('her_ssb_timetrend_timeline.pdf',height=5,width=10)
psidat<-unique(subset(bpdat,var %in% c('Her F','Her SSB','Her CF obs','Her Weight','L Herring')))

her$ssb2<-her$ssb/1000
plot(her$year,her$ssb2,type='l',las=1,xlab='Year',ylab='SSB [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
mx<-subset(her,ssb2==max(her$ssb2,na.rm=TRUE))
mn<-subset(her,ssb2==min(her$ssb2,na.rm=TRUE))
lines(c(mx$year-1,mx$year),c(mx$ssb2,mx$ssb2))
text(mx$year-1.75,mx$ssb2,round(mx$ssb2,digits=0))
text(mn$year,mn$ssb2+40,round(mn$ssb2,digits=0))

plot(her$year,her$ssb2,type='l',las=1,xlab='Year',ylab='SSB [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
rect(1965,-100,1988.5,7e+06,col=alpha('green3',.1),border=NA)
rect(1988.5,-100,1995,7e+06,col=alpha('gold',.1),border=NA)
rect(1995,-100,2006,7e+06,col=alpha('blue',.1),border=NA)
cl<-'red3'
rug(psidat$bpt,lwd=2,col=cl)
text(psidat$bpt[1],0,'Herring Condition Factor declines',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[2],0,'F increase',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[3]-.5,0,'SSB collapse',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[4],0,'Herring weights decline',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[5]+.5,0,'Herring larvae decline',srt=40,adj=0,col=cl,cex=.75)
text(1975,750,'Pre-collapse',col='green3')
text(1991,750,'Collapse',col='gold3')
text(2001,750,'Post-collapse',col='blue3')
mx<-subset(her,ssb2==max(her$ssb2,na.rm=TRUE))
mn<-subset(her,ssb2==min(her$ssb2,na.rm=TRUE))
lines(c(mx$year-1,mx$year),c(mx$ssb2,mx$ssb2))
text(mx$year-1.75,mx$ssb2,round(mx$ssb2,digits=0))
text(mn$year,mn$ssb2+40,round(mn$ssb2,digits=0))

her$r2<-her$r/10000
plot(her$year,her$r2,type='l',las=1,xlab='Year',ylab='R [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
mx<-subset(her,r2==max(her$r2,na.rm=TRUE))
mn<-subset(her,year==2004)
lines(c(mx$year-1,mx$year),c(mx$r2,mx$r2))
text(mx$year-1.75,mx$r2,round(mx$r2,digits=0))
text(mn$year+1,mn$r2+10,round(mn$r2,digits=0))
lines(her$year,her$ssb2,type='l',col='red')
dev.off()





library(scales)
library(akima)
setwd('/scratch/dboyce/spera/data/stagingdat')
her<-read.csv('herring_assess_ssb_r_f_w_spera.csv',header=TRUE)

#AVERAGE F PER YEAR
f<-function(d){
    return(data.frame(f=mean(d$fage2,d$fage3,d$fage4,d$fage5,d$fage6,d$fage7,d$fage8,d$fage9,d$fage10,d$fage11,na.rm=TRUE)))
}
d2<-ddply(her,.(year),.fun=f)
her<-merge(her,d2,by=c('year'),all.x=TRUE,all.y=FALSE)






d<-subset(her,select=c('ssb','year'))
d<-subset(her,select=c('r','year'))


setwd(figsdir)
pdf('herring_yearly_tsplots.pdf',height=12,width=12)
par(mfrow=c(4,3),mar=c(4,4,1,1))
yssb<-modf2(subset(her,select=c('ssb','year')),'SSB',1995,'Tons [x10000]')
yr<-modf2(subset(her,select=c('r','year')),'Recruitment',1995,'Recruitment')
ys<-modf2(subset(her,select=c('s','year')),'Survivorship',1995,'Larval survivorship')
yf<-modf2(subset(her,select=c('f','year')),'F',1995,'Average F')
ylrv<-modf2(subset(her,select=c('larv','year','larvse')),'Larvae',1995,'Larvae [n m3]')
ycf<-modf2(subset(her,select=c('pred.w','year','total.n')),'CF pred',1995,'Predicted condition factor')
ycfo<-modf2(subset(her,select=c('obs.w','year','n.obs.inrange')),'CF obs',1995,'Observed condition factor')
dev.off()







setwd(figsdir)
pdf('herring_stock_breakpointtrends.pdf',height=9,width=8)
par(mfrow=c(3,2),mar=c(2,4,1,1))


nms<-names(her)[19:29]
#cls<-colorRampPalette(brewer.pal(name='Blues',8))(length(nms))
cls<-matlab.like(length(nms))
plot(0,0,col='white',xlim=c(1965,2015),ylim=c(0,.5),las=1,ylab='Weight')
l<-list()
for(i in 1:length(nms)){
    cl<-cls[i]
    d<-na.omit(d)
    d<-subset(her,select=c(paste(nms[i]),'year'))
    names(d)[1]<-'y'
lin.mod <- lm(y~year,data=d)
s<-summary(lin.mod)
segmod <- segmented(lin.mod, seg.Z = ~year)
s2<-summary(segmod)
summary(segmod)
psi<-segmod$psi[2]
out<-data.frame(slope(segmod)$year)
pdat<-data.frame(xx=seq(min(d$year),max(d$year),length.out=1000))
p<-predict(segmod,newdata=data.frame(year=pdat$xx),se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat$upr<-pdat$p+(1.96*pdat$se)
pdat$lwr<-pdat$p-(1.96*pdat$se)

lines(pdat$xx,pdat$p,col=cl)
points(d$year,d$y,pch=16,col=alpha(cl,.5),cex=.75)

#ADD POINT WHERE BREAK IS
pdat$xx<-round(pdat$xx,digits=1)
bp<-subset(pdat,xx==round(psi,digits=1))[1,]
    if((AIC(lin.mod)-AIC(segmod))>4){
        points(psi,bp$p,col=cl,pch=16,cex=2)
        points(psi,bp$p,col='black',pch=1,cex=2,lwd=.1)
} else NULL


out$psi<-segmod$psi[2]
out$var<-nms[i]
out$segtrue<-ifelse((AIC(lin.mod)-AIC(segmod))>4,1,0)
out$r2<-round(s2$adj.r.squared,digits=2)
out$bp<-rownames(out)
l[[i]]<-out
}
ss<-data.frame(do.call('rbind',l))
legend('topright',nms,col=cls,lwd=2,bty='n',ncol=2)


#SAME AS ABOVE BUT EXTRACTS STANDARDIZED SLOPES
nms<-names(her)[19:29]
cls<-matlab.like(length(nms))
plot(0,0,col='white',xlim=c(1965,2015),ylim=c(0,100),las=1,ylab='Weight')
l<-list()
for(i in 1:length(nms)){
    cl<-cls[i]
    d<-na.omit(d)
    d<-subset(her,select=c(paste(nms[i]),'year'))
    names(d)[1]<-'y'
    d$y<-(d$y/max(d$y))*100
lin.mod <- lm(y~year,data=d)
s<-summary(lin.mod)
segmod <- segmented(lin.mod, seg.Z = ~year)
s2<-summary(segmod)
summary(segmod)
psi<-segmod$psi[2]
out<-data.frame(slope(segmod)$year)
pdat<-data.frame(xx=seq(min(d$year),max(d$year),length.out=1000))
p<-predict(segmod,newdata=data.frame(year=pdat$xx),se.fit=TRUE)
pdat$p<-p$fit
lines(pdat$xx,pdat$p,col=cl)

#ADD POINT WHERE BREAK IS
pdat$xx<-round(pdat$xx,digits=1)
bp<-subset(pdat,xx==round(psi,digits=1))[1,]
points(psi,bp$p,col=cl,pch=16,cex=2)
out$psi<-segmod$psi[2]
out$var<-nms[i]
out$segtrue<-ifelse((AIC(lin.mod)-AIC(segmod))>4,1,0)
out$r2<-round(s2$adj.r.squared,digits=2)
out$bp<-rownames(out)
out$pval<-round(s$coef[2,4],digits=4)
l[[i]]<-out
}
s.wgt<-data.frame(do.call('rbind',l))


#BREAKPOINT YEAR BY AGE
dum<-unique(subset(s.wgt,select=c('psi','var','segtrue')))
dum$age<-seq(1,11,1)
dum<-subset(dum,segtrue==1)
plot(dum$age,dum$psi,las=1,pch=15,ylim=c(1980,1990),xlab='Age',ylab='Breakpoint year')
mod<-lm(psi~age,data=dum)
pdat<-data.frame(xx=seq(min(dum$age),max(dum$age),length.out=100))
p<-predict(mod,newdata=data.frame(age=xx),se.fit=TRUE)
pdat$p<-p$fit
pdat$upr<-pdat$p+(1.96*p$se.fit)
pdat$lwr<-pdat$p-(1.96*p$se.fit)
lines(pdat$xx,pdat$p)
polygon(c(pdat$xx,pdat$xx[length(pdat$xx):1]),c(pdat$upr,pdat$lwr[length(pdat$lwr):1]),col=alpha('lightgray',.75),border=NA)


#POST-BREAKPOINT SLOPE BY AGE
dum<-unique(subset(s.wgt,bp=='slope2',select=c('psi','var','bp','Est.','CI.95...l','CI.95...u','segtrue')))
names(dum)<-c('psi','var','bp','slp','lwr','upr','segtrue')
dum$age<-seq(1,11,1)
dum<-subset(dum,segtrue==1)
plot(dum$slp,dum$age,pch=15,yaxt='n',xlim=c(-1.4,2.4),xlab='Change in weight after breakpoint [% yr]',ylab='Age',xaxt='n')
axis(1,seq(-1.4,2.4,.8),las=1)
axis(2,seq(1,11,1),las=1)
abline(v=0,lty=2)
f<-function(d){ lines(c(d$lwr,d$upr),c(d$age,d$age))}
z<-dlply(dum,.(age),.fun=f)
text(2.35,dum$age,labels=gsub(' ','',paste('(',round(dum$psi,digits=0),')')),cex=.8)
dev.off()




plot(her$year,her$r,type='l')

setwd(figsdir)
pdf('herring_stock_trends.pdf',height=9,width=10)
par(mfrow=c(3,2),mar=c(2,4,1,1))
plot(0,0,col='white',ylim=c(0,100),xlim=c(1965,2016),las=1,xlab='Year',ylab='Percent of max')
f<-function(d,cl){
    names(d)[1]<-'y'
    d$y<-(d$y/max(d$y,na.rm=TRUE))*100
    lines(d$year,d$y,lwd=3,col=alpha(cl,.5))
}
f(subset(her,select=c('r','year')),'red')
f(subset(her,select=c('ssb','year')),'royalblue')
f(subset(her,select=c('s','year')),'forestgreen')
f(subset(her,select=c('larv','year')),'purple')
legend('topright',c('r','ssb','s','larvae'),col=c('red3','royalblue','forestgreen','purple'),lwd=2,bty='n')


plot(0,0,col='white',ylim=c(0,100),xlim=c(1965,2016),las=1,xlab='Year',ylab='Percent of max')
f<-function(d,cl){
    d<-na.omit(d)
    names(d)[1]<-'y'
    d$y<-(d$y/max(d$y,na.rm=TRUE))*100
    xo<-seq(min(d$year),max(d$year),length.out=1000)
    asp<-aspline(d$year,d$y,xout=xo)
    lines(asp$x,asp$y,lwd=3,col=alpha(cl,.5))
}
f(subset(her,select=c('r','year')),'red')
f(subset(her,select=c('ssb','year')),'royalblue')
f(subset(her,select=c('s','year')),'forestgreen')
f(subset(her,select=c('larv','year')),'purple')
legend('topright',c('r','ssb','s','larvae'),col=c('red3','royalblue','forestgreen','purple'),lwd=2,bty='n')



#cls<-colorRampPalette(c('magenta4','blue3','green','yellow','red3'))
#cls<-(cls(length(lbls)))
library(RColorBrewer)
nms<-names(her)[5:15]
#cls<-colorRampPalette(brewer.pal(name='Blues',8))(length(nms))
cls<-matlab.like(length(nms))
plot(0,0,col='white',xlim=c(1965,2005),ylim=c(0,4.5),las=1,ylab='F')
for(i in 1:length(nms)){
    d<-na.omit(d)
    d<-subset(her,select=c(paste(nms[i]),'year'))
    names(d)[1]<-'y'
    xo<-seq(min(d$year),max(d$year),length.out=1000)
    asp<-aspline(d$year,d$y,xout=xo)
    lines(asp$x,asp$y,col=alpha(cls[i],.75),lwd=2)
}
legend('topleft',nms,col=cls,lwd=2,bty='n',ncol=2)

plot(0,0,col='white',xlim=c(1965,2005),ylim=c(.01,4.5),las=1,log='y',ylab='F')
for(i in 1:length(nms)){
    d<-na.omit(d)
    d<-subset(her,select=c(paste(nms[i]),'year'))
    names(d)[1]<-'y'
    lines(d$year,d$y,col=alpha(cls[i],.75),lwd=2)
}



nms<-names(her)[19:29]
cls<-matlab.like(length(nms))
plot(0,0,col='white',xlim=c(1965,2015),ylim=c(0,.55),las=1,ylab='Weight')
for(i in 1:length(nms)){
    d<-na.omit(d)
    d<-subset(her,select=c(paste(nms[i]),'year'))
    names(d)[1]<-'y'
    xo<-seq(min(d$year),max(d$year),length.out=1000)
    asp<-aspline(d$year,d$y,xout=xo)
    lines(asp$x,asp$y,col=alpha(cls[i],.75),lwd=2)
}
legend('topright',nms,col=cls,lwd=2,bty='n',ncol=2)

plot(0,0,col='white',xlim=c(1965,2015),ylim=c(0,.55),las=1,ylab='Weight')
for(i in 1:length(nms)){
    d<-na.omit(d)
    d<-subset(her,select=c(paste(nms[i]),'year'))
    names(d)[1]<-'y'
    lines(d$year,d$y,col=alpha(cls[i],.75),lwd=2)
}
dev.off()






map.axes()
points(a2$lon,a2$lat,pch=16,col='red')
