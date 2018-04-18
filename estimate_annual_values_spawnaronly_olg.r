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

setwd(datadir)
hdata<-read.csv('herring_landings_historical_lotze.csv',header=TRUE)
hdata<-subset(hdata,is.na(year)==FALSE)
plot(hdata,pch=16)
par(mfrow=c(2,3))
plot(hdata$year,hdata$her.land.50,type='b')
plot(hdata$year,hdata$her.land.51,type='b')
plot(hdata$year,hdata$her.land.52,type='b')
plot(hdata$year,hdata$her.land.53,type='b')
plot(hdata$year,hdata$her.land.50.53,type='b')
plot(hdata$year,hdata$herjuv.land.50.53,type='b')


par(mfrow=c(2,3))
plot(hdata$year,hdata$her.land.50.53,type='b',xlim=c(1860,1945))
lines(hdata$year,hdata$herjuv.land.50.53,type='b',xlim=c(1860,1945),col='red')
lines(hdata$year,hdata$her.land.t50.53,type='b',xlim=c(1860,1945),col='green')

par(mfrow=c(2,3))
plot(hdata$year,hdata$her.land.t50,type='b')
plot(hdata$year,hdata$her.land.t51,type='b')
plot(hdata$year,hdata$her.land.t52,type='b')
plot(hdata$year,hdata$her.land.t53,type='b')
#plot(hdata$year,hdata$her.land.t50.53,type='b')


setwd(datadir)
load('SPERA_andata.RData')
a<-subset(hdata,select=c('year','her.land.t50.53'))
dat<-merge(data,a,by=c('year'),all=TRUE)
plot(dat$year,dat$her.land.t50.53*10,type='b',col='red',ylim=c(20000,2000000),pch=15)
lines(dat$year,dat$her.land,type='b',pch=16)

plot(dat$her.land,dat$her.land.t50.53,pch=16)
summary(dat$her.land.t50.53*100)
150*1000
summary(dat$her.land.t50.53)

#LON/LAT FOR GERMAN BANK
gblon<--66.3
gblat<-43.3


#CALCULATE SUM OF COUNTS (NUMBER OF INDIVIDUALS) PER SAMPLE; RUNS MORE THAN ONCE SO NEED TO SOURCE
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



###############################################################

###############################################################
#SPATIAL STRUCTURE ANALYSES FOR LARVAL AND ADULT HERRING: TAKE AWHILE TO RUN SO RUNS ONCE AND SAVES
#JUVENILE HERRING
setwd(datadir)
load("BoF_larval_herring_counts_spera_spawnar.RData")
#SPAWNING ON GERMAN BANK AUG-SEPT, EARLIER ON LURCHER
larv<-subset(dat,bathy<0 & spec_stage %in% c('J') & month>=8 & gearprocdesc=='Sawtooth Oblique')
#STANDARDIZE PER VOLUME STRAINED
larv$stdno<-(larv$totno/larv$volm3)*abs(larv$bathy)#NUMBER PER M2 - INTEGRATED OVER DEPTH
#larv$stdno<-(larv$totno/larv$volm3)
larv<-subset(larv,is.na(stdno)==FALSE )
larv$phylum<-tolower(larv$phylum)
larv$order<-tolower(larv$order)
larv$kingdom<-tolower(larv$kingdom)
#ADD ID FOR SAMPLES
larv$sid<-gsub(' ','',paste(larv$cruise,'_',larv$setno,'_',larv$sample))

#SUBSET ONLY STATIONS THAT ARE SAMPLED IN >=20 OF 23 YEARS
larv$station<-as.character(larv$station)
dm<-data.frame(station=sort(unique(larv$station)),
       nyr=tapply(larv$year,larv$station,function(x) length(unique(x))))
dm2<-subset(dm,nyr>=19)
larv<-subset(larv,station %in% dm2$station)

jher<-ddply(subset(larv,comname %in% c('herring/capelin like','herring(atlantic)')),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
jher$occ<-factor(ifelse(jher$stdno>0,TRUE,FALSE))




#HERRING LARVAE
setwd(datadir)
load("BoF_larval_herring_counts_spera_spawnar.RData")
#SPAWNING ON GERMAN BANK AUG-SEPT, EARLIER ON LURCHER
larv<-subset(dat,bathy<0 & spec_stage %in% c('L') & month>=8 & gearprocdesc=='Sawtooth Oblique')
#STANDARDIZE PER VOLUME STRAINED
larv$stdno<-(larv$totno/larv$volm3)*abs(larv$bathy)#NUMBER PER M2 - INTEGRATED OVER DEPTH
#larv$stdno<-(larv$totno/larv$volm3)
larv<-subset(larv,is.na(stdno)==FALSE )
larv$phylum<-tolower(larv$phylum)
larv$order<-tolower(larv$order)
larv$kingdom<-tolower(larv$kingdom)
#ADD ID FOR SAMPLES
larv$sid<-gsub(' ','',paste(larv$cruise,'_',larv$setno,'_',larv$sample))

#SUBSET ONLY STATIONS THAT ARE SAMPLED IN >=20 OF 23 YEARS
larv$station<-as.character(larv$station)
dm<-data.frame(station=sort(unique(larv$station)),
       nyr=tapply(larv$year,larv$station,function(x) length(unique(x))))
dm2<-subset(dm,nyr>=19)
larv<-subset(larv,station %in% dm2$station)

lher<-ddply(subset(larv,comname %in% c('herring/capelin like','herring(atlantic)')),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
lher$occ<-factor(ifelse(lher$stdno>0,TRUE,FALSE))
#lher$stdno<-lher$stdno+.01


##############################################################
#ESTIMATE SPATIAL STRUCTURE FOR EACH YEAR- TAKES AWHILE TO RUN SO RUN ONCE AND SAVE, THEN RE-IMPORT
#USE HERRING OVER ENTIRE SURVEY DOMIAN FOR SPATIAL ANALYSES
setwd(datadir)
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw<-subset(rvw,month %in% c(6,7,8))
rvw$strat<-as.numeric(as.character(rvw$strat))
rvw<-subset(rvw,strat>=480 & strat<=495)


#SUBSET ONLY STRATA THAT ARE SAMPLED IN >45 OF 47 YEARS
rvw$strat<-as.character(rvw$strat)
dm<-data.frame(strat=sort(unique(rvw$strat)),
       nyr=tapply(rvw$year,rvw$strat,function(x) length(unique(x))))
dm2<-subset(dm,nyr>=45)
rvw<-subset(rvw,strat %in% dm2$strat)



ff<-function(d){
names(d)[1]<-'y'
d$y<-log10(d$y+.01)
print(unique(d$year))
f<-function(a){return(data.frame(y=median(a$y)))}
d<-ddply(d,.(lon,lat),.fun=f)

if(dim(d)[1]>10 & length(unique(d$y))>1){
crds<-as.matrix(subset(d,select=c('lon','lat')))
sp<-spatialProcess(crds,d$y)
spcv<-(sd(d$y)/(mean(d$y)+10))*100#SPATIAL VARIANCE RE. LITZOW 2008
out<-data.frame(nug=sp$sigma.MLE,#NUGGET
                spvar=sp$rho.MLE,#PROCESS VARIANCE
                rng=sp$theta.MLE,
                spcv=spcv,
                n=dim(d)[1])#RANGE PARAMETER
#set.panel(2,2)#plot(sp)
return(out)
} else NULL
}
spjuv<-ddply(subset(jher,select=c('stdno','lon','lat','year')),.(year),.fun=ff,.progress='text')
names(spjuv)<-c('year','herjuv.spnug','herjuv.spvar','herjuv.sprng','herjuv.spcv','herjuv.sp.n')
splarv<-ddply(subset(lher,select=c('stdno','lon','lat','year')),.(year),.fun=ff,.progress='text')
names(splarv)<-c('year','herlrv.spnug','herlrv.spvar','herlrv.sprng','herlrv.spcv','herlrv.sp.n')
spher<-ddply(subset(rvw,select=c('totwgt','lon','lat','year')),.(year),.fun=ff,.progress='text')
names(spher)<-c('year','her.spnug','her.spvar','her.sprng','her.spcv','her.sp.n')


library(rgeos)
library(rgdal)
mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast<-readOGR('N:/data/shapefiles/naturalearthdata_ne_10m_land_poly',layer='ne_10m_land')#works for Alfcoast<-fortify(coast)
coast.mc<-crop(coast,extent(-180,180,-90,90),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#fcoast.mc<-fortify(coast.mc)


#CREATE POLYGON AROUND RV SURVEY DOMAIN
aa<-unique(subset(rvw,select=c('lon','lat')))
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$lon,y=aa$lat)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=.2, byid=TRUE)
plgrv <- gUnaryUnion(buf1)
#buf <-  gBuffer(buf2, width=0)
plgrv<-erase(plgrv,coast.mc)
#plot(plgrv)


#CREATE POLYGON AROUND LARVAL HERRING SURVEY DOMAIN
aa<-unique(subset(lher,select=c('lon','lat')))
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$lon,y=aa$lat)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=.22, byid=TRUE)
plgls <- gUnaryUnion(buf1)
plgls<-erase(plgls,coast.mc)
#plot(plgls)

#ESTIMATE GEOGRAPHIC RANGE CHANGES
f<-function(d){
nm<-names(d)[1]
names(d)[1]<-'y'
d$y2<-ifelse(d$y>0,TRUE,FALSE)
if(dim(d)[1]>50){
if(nm=='stdno'){
mod<-gam(y2~s(lon,lat,k=50),data=d,family=binomial)
pdat<-expand.grid(lon=seq(min(lher$lon),max(lher$lon),length.out=100),lat=seq(min(lher$lat),max(lher$lat),length.out=100))
} else {
mod<-gam(y2~s(lon,lat,k=25),data=d,family=binomial)
pdat<-expand.grid(lon=seq(min(rvw$lon),max(rvw$lon),length.out=100),
                 lat=seq(min(rvw$lat),max(rvw$lat),length.out=100))
}
p<-predict(mod,newdata=pdat,type='response')
pdat$p<-p

#RESTRICT TO SURVEY DOMAIN
crds<-SpatialPoints(data.frame(lon=pdat$lon,lat=pdat$lat),proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
if(nm=='stno'){
pdat$ov<-over(crds,plgls)
} else {
pdat$ov<-over(crds,plgrv)
}

pdat<-subset(pdat,ov==1)
pdat2<-subset(pdat,p>.95)
out<-data.frame(rng=(dim(pdat2)[1]/dim(pdat)[1])*100,
                n=dim(d)[1])
return(out)
} else NULL
}
rnglher<-ddply(subset(lher,select=c('stdno','lon','lat','year')),.(year),.fun=f,.progress='text')
names(rnglher)<-c('year','herlrv.georng','herlrv.georng.n')
rngher<-ddply(subset(rvw,select=c('totno','lon','lat','year')),.(year),.fun=f,.progress='text')
names(rngher)<-c('year','her.georng','her.georng.n')

herspat2<-merge(splarv,spher,by=c('year'),all=TRUE)
herspat2<-merge(herspat2,spjuv,by=c('year'),all=TRUE)
herspat2<-merge(herspat2,rnglher,by=c('year'),all=TRUE)
herspat2<-merge(herspat2,rngher,by=c('year'),all=TRUE)
#save(herspat2,file='herspat2.RData')
#load('herspat2.RData')



#############################################################

# READS IN DATA THAT IS SPATIO-TEMPORAL AND CALCULATES ANNUAL AVERAGES AND VARIANCES FROM GAM (ANDATA)

##############################################################
setwd(datadir)
load('phytoplankton_biomass_all_spera_spawnar.RData');phyt<-dat
path<-load('pathfinder_sst_daily_1981_2012_spawnar.RData');path<-dat
nit<-read.csv("Atlas1999nuts_nit_spawnar.csv",header=TRUE)
phos<-read.csv("Atlas1999nuts_phos_spawnar.csv",header=TRUE)
sil<-read.csv("Atlas1999nuts_sil_spawnar.csv",header=TRUE)
strat<-read.csv("physical_stratification_spera_spawnar.csv",header=TRUE)
#rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
#rvw<-subset(rvw,month %in% c(6,7,8))
zto<-read.csv("zoop_turnover_spera_spawnar.csv",header=TRUE)
names(zto)[5]<-'zp.to.cpr'
pto<-read.csv("phyto_turnover_spera_spawnar.csv",header=TRUE)
names(pto)[5]<-'ph.to.cpr'
plank<-read.csv('plank_cpr_richness_spera_spawnar.csv',header=TRUE)
names(plank)[8:19]<-c('zp.rich.cpr','zp.rich1.cpr','zp.div.cpr','zp.pe.cpr','zp.mn.cpr','zp.nspec.cpr','zp.cop.cpr','ph.rich.cpr','ph.rich1.cpr','ph.div.cpr','ph.pe.cpr','ph.mn.cpr')


###############################################################
#IMPORT AND FORMAT LARVAL HERRING OBS FOR FUTURE AVERAGING;
#EXPORT LARVAL HERRING DATA FOR THERMAL ENVELOPE CALCS
load("BoF_larval_herring_counts_spera_spawnar.RData")
larv<-subset(dat,bathy<0 & spec_stage=='L')
#STANDARDIZE PER VOLUME STRAINED
larv$stdno<-(larv$totno/larv$volm3)*abs(larv$bathy)#NUMBER PER M2 - INTEGRATED OVER DEPTH
#larv$stdno<-(larv$totno/larv$volm3)
larv<-subset(larv,is.na(stdno)==FALSE)
larv<-subset(larv,gearprocdesc=='Sawtooth Oblique')
#SPAWNING ON GERMAN BANK AUG-SEPT, EARLIER ON LURCHER
larv<-subset(larv,month>=7)
larv$phylum<-tolower(larv$phylum)
larv$order<-tolower(larv$order)
larv$kingdom<-tolower(larv$kingdom)

#ADD ID FOR SAMPLES
larv$sid<-gsub(' ','',paste(larv$cruise,'_',larv$setno,'_',larv$sample))

#OUTPUT HERRING LARVAE DATA
herlrv.mn.bof<-ddply(subset(larv,month>5 & comname=='herring(atlantic)'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')

#OUTPUT HERRING LARVAE DATA FOR THERMAL ENVELOPE FIGURES
lher<-ddply(subset(larv,comname=='herring(atlantic)'),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
lher$occ<-factor(ifelse(lher$stdno>0,TRUE,FALSE))
lher$stdno<-lher$stdno+.01
lher<-subset(lher,is.na(surftemp)==FALSE)
#save(lher,file='lher.RData')



#CALCULATE RICHNESS, EVENNESS, FOR LARVAE
#NEED TO REMOVE 0'S OR CAN'T CALCULATE SHANNON INDEX
load("BoF_larval_herring_counts_spera_spawnar.RData")
larv<-subset(dat,bathy<0 & month>6 & spec_stage=='L')
#STANDARDIZE PER VOLUME STRAINED
larv$stdno<-(larv$totno/larv$volm3)*abs(larv$bathy)#NUMBER PER M2 - INTEGRATED OVER DEPTH
#larv$stdno<-(larv$totno/larv$volm3)
larv<-subset(larv,is.na(stdno)==FALSE)
larv<-subset(larv,gearprocdesc=='Sawtooth Oblique')
#SPAWNING ON GERMAN BANK AUG-SEPT, EARLIER ON LURCHER
larv$phylum<-tolower(larv$phylum)
larv$order<-tolower(larv$order)
larv$kingdom<-tolower(larv$kingdom)
#ADD ID FOR SAMPLES
larv$sid<-gsub(' ','',paste(larv$cruise,'_',larv$setno,'_',larv$sample))

#################################################################
#SUM OF COUNTS (NUMBER OF INDIVIDUALS) PER SAMPLE FOR DIFFERENT LARVAL CATEGORIES
avfun1<-function(d){
dout<-unique(subset(d,select=c('tday','tbin3','tbin5','tbin10','surftemp','bottemp','airtemp','sid')))
dout$stdno<-sum(d$stdno)

dout$var<-'T zoop'
return(dout)
}
#PREF PREY OF JUVENILLE HERRING (STEVENSON 2005 AND REFS THEREIN)
zp1<-ddply(subset(larv,order %in% c('calanoida','cyclopoida','decapoda') | phylum %in% c('mollusca')),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun1,.progress='text')
zp1$var<-'herjuv.prey.bof'

#PREF PREY OF ADULT HERRING (STEVENSON 2005 AND REFS THEREIN)
zp2<-ddply(subset(larv,order %in% c('euphausiacea','calanoida','cyclopoida') | phylum %in% c('chaetognatha')),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun1,.progress='text')
zp2$var<-'her.prey.bof'

#TOTAL PREY OF ALL HERRING (LARVAE AND ZOOP)
zp<-ddply(larv,.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun1,.progress='text')
zp$var<-'her.preytot.bof'
#predprey.bof<-rbind.fill(zp1,zp2,zp)
predprey.bof<-rbind.fill(zp1,zp)




########################################################
#CALCULATE RICHNESS/DIVERISTY/EVENNESS OF ALL LARVAE EXCLUDING HERRING
larv2<-subset(larv,stdno>0 & !(comname %in% c('herring(atlantic)','herring/capelin like')))
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
lrich.bof<-ddply(subset(larv2,kingdom=='animalia',select=c('year','month','day','tday','lon','lat','tbin3','tbin5','tbin10','valid_name','bathy','stdno')),.(year,month,day,tday,tbin3,tbin5,tbin10,bathy,lon,lat),.fun=richfun,.progress='text')
names(lrich.bof)[11:16]<-gsub(' ','',paste(names(lrich.bof)[11:16],'.lrv'))



##########################################
#CPR ZOOPLANKTON: RETAIN LARGEST SPECIES
zoop<-read.csv("zoop_cpr_counts_spera_spawnar.csv",header=TRUE)
zoop<-subset(zoop,scientificname %in% c('Calanoida','Oithona','Calanus', 'Calanus finmarchicus','Calanus glacialis','Calanus hyperboreus','Candacia','Candacia armata','Centropages','Centropages hamatus','Centropages typicus','Pseudocalanus elongatus','Temora longicornis','Copepoda'))

#SUM COUNTS (INDIVIDUALS) PER SAMPLE
f<-function(d){  return(data.frame(zoop=sum(d$counts)))}
zoop<-ddply(zoop,.(lon,lat,year, month,djul,day,tbin3,tbin5,tbin10),.fun=f,.progress='text')



####################################################################
##MODELING FUNCTIONS: FIRST FUNCTION TAKES DATA AND FIRST ESTIMATES ANNUAL AVERAGES, SECOND JUST TAKES ANNUAL AVERAGES
 modf.ave<-function(d,nm){
    options(warn=-1)
    nmm<-names(d)[1]
    d<-subset(d,year>1960)
    names(d)[1]<-'y'
    d<-na.omit(d)
    #ESTIMATE ANNUAL AVERAGES
    if(nm %in% c('herlr')){
        print('lognormal')
        d$y<-d$y+.1
    mod<-gam(y~as.factor(year)+s(lon,lat,k=30),data=d,gamma=1.4,family='gaussian'(link='log'))
    } else {
        mod<-gam(y~as.factor(year) + s(day,bs='cs',k=5) + s(lon,lat,k=30),data=d,gamma=1.4)
    }

    #PREDICT ANNUAL AVERAGES AND CIS
    pdat<-data.frame(year=sort(unique(d$year)),
                     lon=gblon,
                     lat=gblat,
                     day=250)
    p<-predict(mod,newdata=pdat,se.fit=TRUE,type='response')
    pdat$p<-p$fit
    pdat$se<-p$se.fit
    pdat<-subset(pdat,p!=Inf)

    dum<-data.frame(year=sort(unique(d$year)),
    n=tapply(d$y,d$year,length),
    nmonths=tapply(d$month,d$year,function(x) length(unique(x))))
    pdat<-merge(pdat,dum,by=c('year'),all.x=TRUE,all.y=FALSE)
#    if(nm %in% c('Nit','Sil','Phos','Lrv Fish','Strat','Temp 50m','Zoop Cells')){
#            pdat<-subset(pdat,nmonths>0)
#    } else { pdat<-subset(pdat,nmonths>1)
#    }

pdat<-subset(pdat,select=c('year','p','se'))
    names(pdat)[2]<-nm
    names(pdat)[3]<-gsub(' ','.',paste(nm,'se'))
    return(pdat)
    options(warn=0)
}

l<-list()
l[[1]]<-modf.ave(subset(path,select=c('sst','year','day','lon','lat','month')),'sst.mn')
l[[2]]<-modf.ave(subset(path,select=c('wind','year','day','lon','lat','month')),'wnd.mn')
l[[3]]<-modf.ave(subset(strat,select=c('s25','year','day','lon','lat','month')),'strt')
l[[4]]<-modf.ave(subset(strat,select=c('temp50','year','day','lon','lat','month')),'t50')#R WITH SST: 76
l[[5]]<-modf.ave(subset(predprey.bof,var=='herjuv.prey.bof',select=c('stdno','year','day','lon','lat','month')),'herjuv.prey.bof')
#l[[6]]<-modf.ave(subset(predprey.bof,var=='her.prey.bof',select=c('stdno','year','day','lon','lat','month')),'her.prey.bof')
l[[6]]<-modf.ave(subset(herlrv.mn.bof,select=c('stdno','year','day','lon','lat','month')),'herlrv.mn.bof')
l[[7]]<-modf.ave(subset(predprey.bof,var=='her.preytot.bof',select=c('stdno','year','day','lon','lat','month')),'her.preytot.bof')
l[[8]]<-modf.ave(subset(lrich.bof,select=c('shannon.lrv','year','day','lon','lat','month')),'lrv.rich.bof')
l[[9]]<-modf.ave(subset(lrich.bof,select=c('pe.lrv','year','day','lon','lat','month')),'lrv.pe.bof')
l[[10]]<-modf.ave(subset(lrich.bof,select=c('totalabundance.lrv','year','day','lon','lat','month')),'lrv.mn.bof')
l[[11]]<-modf.ave(subset(zoop,select=c('zoop','year','day','lon','lat','month')),'zp.cop.cpr')
l[[12]]<-modf.ave(subset(nit,select=c('nit','year','day','lon','lat','month')),'nit')
l[[13]]<-modf.ave(subset(sil,select=c('sil','year','day','lon','lat','month')),'sil')#
l[[14]]<-modf.ave(subset(phos,select=c('phos','year','day','lon','lat','month')),'phos')
l[[15]]<-modf.ave(subset(plank,select=c('ph.div.cpr','year','day','lon','lat','month')),'ph.div.cpr')
l[[16]]<-modf.ave(subset(plank,select=c('ph.mn.cpr','year','day','lon','lat','month')),'ph.mn.cpr')
l[[17]]<-modf.ave(subset(plank,select=c('ph.rich.cpr','year','day','lon','lat','month')),'ph.rich.cpr')
l[[18]]<-modf.ave(subset(plank,select=c('zp.pe.cpr','year','day','lon','lat','month')),'zp.pe.cpr')
l[[19]]<-modf.ave(subset(plank,select=c('ph.pe.cpr','year','day','lon','lat','month')),'ph.pe.cpr')
l[[20]]<-modf.ave(subset(plank,select=c('zp.mn.cpr','year','day','lon','lat','month')),'zp.mn.cpr')
l[[21]]<-modf.ave(subset(plank,select=c('zp.rich.cpr','year','day','lon','lat','month')),'zp.rich.cpr')
l[[22]]<-modf.ave(subset(plank,select=c('zp.div.cpr','year','day','lon','lat','month')),'zp.div.cpr')
andata<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), l)#COMBINE


#LOAD ADDITIIONAL DATA AND MERGE TO ONE
setwd(datadir)
load('herspat2.RData')
load('plank_standrews.RData')

#load('sstpath.RData')
andata<-merge(andata,herspat2,by=c('year'),all=TRUE)
andata<-merge(andata,plank_standrews,by=c('year'),all=TRUE)
#andata<-merge(andata,sstpath,by=c('year'),all=TRUE)
andata<-list(andata)




##############################################################

##############################################################
#################################################
setwd(datadir)
nao<-read.csv('nao.csv',header=TRUE)
names(nao)[3]<-'nao.se'

#POSITION OF SLOPE WATERS DURING FALL
setwd(datadir)
cur<-read.csv('GS_SS_distances.csv',header=TRUE)
f<-function(d){
    return(data.frame(gs.dist=mean(d$GS.Dist),
               ss.dist=mean(d$SS.Dist),
               gs.dist.se=sd(d$GS.Dist),
               ss.dist.se=sd(d$SS.Dist)))
}
cur<-ddply(subset(cur,month>=8 & month<12),.(year),.fun=f,.progress='text')

#ACCUMULATOR
andata<-c(andata,list(merge(cur,nao,by=c('year'),all=TRUE)))



#PREDICT ANNUAL AVERAGES FOR HERRING LARVAE LENGTH
setwd(datadir)
load("BoF_larval_herring_lengths_spera_spawnar.RData");larvs<-dat
larvs<-subset(larvs,gearprocdesc=='Sawtooth Oblique' & spec_stage=='L' & month>8)

#SUBSET ONLY STATIONS THAT ARE SAMPLED IN >=20 OF 23 YEARS
larvs$station<-as.character(larvs$station)
dm<-data.frame(station=sort(unique(larvs$station)),
       nyr=tapply(larvs$year,larvs$station,function(x) length(unique(x))))
dm2<-subset(dm,nyr>=19)
larvs<-subset(larvs,station %in% dm2$station)

#load("BoF_larval_herring_counts_spera_spawnar.RData")
#dm<-unique(subset(larv,select=c('lon','lat','bathy')))
#larvs<-merge(larvs,dm,by=c('lon','lat'),all.x=TRUE,all.y=FALSE)
#larvs<-subset(larvs,bathy<0)
#STANDARDIZE PER VOLUME STRAINED
#larvs$clen<-larvs$clen*abs(larvs$bathy)#NUMBER PER M2 - I
modf<-gam(length~as.factor(year) + s(lon,lat,k=20) + s(day,k=5),weights=larvs$clen,data=larvs,gamma=1.4)
herlrvlen<-data.frame(year=sort(unique(larvs$year)),
                 lon=gblon,
                 lat=gblat,
                 day=230)
p<-predict(modf,newdata=herlrvlen,se.fit=TRUE)
herlrvlen$herlrv.len<-p$fit
herlrvlen$herlrv.len.se<-p$se.fit

andata<-c(andata,list(subset(herlrvlen,select=c('year','herlrv.len','herlrv.len.se'))))



##############################################
#AGE DIVERSITY OF THE CATCH
library(tidyr)
setwd(datadir)
cdat<-read.csv('herring_assess2016_catchatage_pct.csv',header=TRUE)
cdat<-cdat[,1:12]
a<- cdat %>% gather(age, y, X1:X11)
names(a)[1]<-'year'
a$y<-as.numeric(as.character(a$y))

richfun<-function(d){
d<-subset(d,is.na(age)==FALSE & y>0)
dout<-unique(subset(d,select=c('year')))

#RICHNESS WITH 1% CUTOFF LOCALLY
tcells<-sum(d$y)#TOTAL LANDINGS

#SHANNON I
f2<-function(d3){ return(data.frame(p=d3$y/tcells)) }
sdat<-ddply(subset(d,select=c('age','y')),.(age),.fun=f2)
sdat$lp<-log(sdat$p)
sdat$pp<--1*(sdat$p*sdat$lp)
sdat<-subset(sdat,is.na(pp)==FALSE)

dout$shannon<-sum(sdat$pp,na.rm=TRUE)
dout$pe<-sum(sdat$pp)/(log(length(unique(d$age))))
return(dout)
}
her.adiv.cat<-ddply(a,.(year),.fun=richfun,.progress='text')

#HERRING AVERAGE AGE DIVERSITY AND EVENNESS FROM CATCH
her.adiv.cat<-subset(her.adiv.cat,select=c('year','shannon'))
her.adiv.cat$her.agediv.cat<-her.adiv.cat$shannon

#ACCUMULATOR
andata<-c(andata,list(her.adiv.cat))








######################################################
#HERRING DATA FROM ASSESSMENTS: BIOMASS, WEIGHT, F, ETC...
setwd(datadir1)
her<-read.csv('herring_assess_ssb_r_f_w_spera.csv',header=TRUE)
names(her)<-tolower(names(her))
#LANDINGS ARE FOR ALL AGES(1-11) AND SSB IS ONLY FOR AGES 4-8: THIS IS WHY EXPLOITATION RATE EXCEEDS 1 IN MANY CASES
her$expr<-(her$landings/10)/her$ssb
her$expr4wx<-her$x4wx.stock.nominal.landings/her$ssb

#PREDICT SSB FROM ACOUSTIC- ADJUST ACOUSTIC DATA DOWNWARD TO MATCH VPA
mod<-lm(ssb~acoustic,data=her)
her$acoustic2<-predict(mod,newdata=data.frame(acoustic=her$acoustic))

#PREDICT ACOUSTIC FROM SSB- ADJUST SSB DATA UPWARD TO MATCH ACOUSTIC
mod<-lm(acoustic~ssb,data=her)
her$ssb2<-predict(mod,newdata=data.frame(ssb=her$ssb))


#AVERAGE F PER YEAR
f<-function(d){
d<-na.omit(d)
return(data.frame(f=mean(c(d$fage2,d$fage3,d$fage4,d$fage5,d$fage6,d$fage7,d$fage8,d$fage9,d$fage10,d$fage11),na.rm=TRUE)))
}
d2<-ddply(subset(her,is.na(year)==FALSE,select=c('year','fage1','fage3','fage4','fage5','fage6','fage7','fage8','fage9','fage10','fage11')),.(year),.fun=f)
her<-merge(her,d2,by=c('year'),all.x=TRUE,all.y=FALSE)

#AVERAGE WEIGHT PER YEAR
f<-function(d){
 return(data.frame(wgt=mean(c(d$wage7,d$wage8,d$wage9,d$wage10,d$wage11),na.rm=TRUE)))
}
d2<-ddply(her,.(year),.fun=f)
her<-merge(her,d2,by=c('year'),all.x=TRUE,all.y=FALSE)

plot(her$year,her$f,pch=15)
her<-subset(her,is.na(year)==FALSE,select=c('year','r','ssb','s','obs.w','pred.w','f','wgt','expr4wx','x4wx.stock.nominal.landings'))
names(her)<-c('year','her.rec1','her.ssb','herlrv.surv','her.cf.rv','her.cf.cat','her.f','her.waa','her.expr','her.land')
her<-subset(her,select=c('year','her.rec1','her.ssb','herlrv.surv','her.cf.rv','her.cf.cat','her.waa','her.land'))



########################################################
#ADDS FISHING RATE (CATCH/TOTAL BIOMASS), AND SURPLUS PRODUCTION
setwd(datadir)
hdat<-read.csv('herring_assess2006_ssb.csv',header=TRUE)
a<-na.omit(subset(her,select=c('year','her.land')))
aa<-subset(hdat,select=c('year','totbiotons','ssbtons'))
b<-merge(a,aa,by=c('year'),all=FALSE)
b$deltat<-c(diff(b$totbiotons),NA)
b$her.prod<-b$deltat+b$her.land
b$her.expr<-b$her.land/b$totbiotons

b<-subset(b,select=c('year','totbiotons','her.prod','her.expr'))
names(b)[2]<-c('her.tbio')
her<-merge(her,b,by=c('year'),all=TRUE)

#ACCUMULATOR
andata<-c(andata,list(her))






#######################################
#ABUNDANCE AT AGE FROM LANDINGS
setwd(datadir)
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
#plot(tdat$age,tdat$beta,pch=15,las=1,xlab='Age',ylab='Trend in number landed over time',xaxt='n')
#axis(1,seq(1,11,1))


#TREND IN PROPORTION OF OLD HERRING IN CATCH OVER TIME
f<-function(d){
    a<-subset(d,age>=5)
    return(data.frame(her.age5.cat=(sum(a$catch)/sum(d$catch))*100))
}
tdat2<-ddply(hera,.(year),.fun=f)
#plot(tdat2$year,tdat2$her.age5.cat,pch=15,type='b',las=1,xlab='Year',ylab='Proportion of age 4+ herring in catch')


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

dout<-data.frame(her.agepe.cat=sum(sdat$pp)/(log(length(unique(d2$age)))))
return(dout)
}
tdat3<-ddply(hera,.(year),.fun=f)
plot(tdat3$year,tdat3$age.pe,pch=15,type='b',xlab='Year',ylab='Evenness of catch age composition',las=1)

#ABUNDANCE WEIGHTED AGE OVER TIME
mod<-lm(age~as.factor(year),weights=hera$catch,data=hera)
pdat<-data.frame(year=sort(unique(hera$year)))
p<-predict(mod,newdata=data.frame(year=pdat$year))
pdat$her.age.cat<-p
tdat4<-pdat
#plot(tdat4$year,tdat4$her.age.cat,type='b',pch=15,las=1)

agedat<-merge(tdat2,tdat3,by=c('year'),all=FALSE)
agedat<-merge(agedat,tdat4,by=c('year'),all=FALSE)

#ADD TO ACCUMULATOR
andata<-c(andata,list(agedat))



########################################################
#AVERAGE LENGTH OF HERRING FROM RV
#CORRECTION AS GIVEN BY MIKE MCMAHON: 2016 SWITCH FROM CM TO MM AND FROM FORK LENGTH TO TOTAL LENGTH
setwd(datadir)
rvl<-read.csv("herring_lengths_RV_survey_spera_spawnar.csv",header=TRUE)
rvl$strat<-as.numeric(as.character(rvl$strat))
rvl<-subset(rvl,strat>=480 & strat<=495)
rvl$flen<-ifelse(rvl$mission=='NED2016016',rvl$flen/10,rvl$flen)
rvl$flen<-ifelse(rvl$mission!='NED2016016',(rvl$flen*1.0866)+0.95632,rvl$flen)
#plot(log10(rvl$flen),log10(rvl$fwt))
rvl<-subset(rvl,log10(rvl$flen)>=.5 & month %in% c(6,7,8))#REMOVE OUTLIERS

#HERRING AVERAGE LENGTH FROM RV SURVEY
mod<-gam(flen~as.factor(year) + s(lon,lat,k=50),weights=rvl$clen,data=rvl,gamma=1)
hlen<-data.frame(year=sort(unique(rvl$year)),
                 lon=gblon,
                 lat=gblat)
p<-predict(mod,newdata=hlen,se.fit=TRUE)
hlen$her.len.rv<-p$fit
hlen$her.len.rv.se<-p$se.fit
#plot(hlen$year,hlen$her.len.rv,pch=15)

#ACCUMULATOR
andata<-c(andata,list(subset(hlen,select=c('year','her.len.rv','her.len.rv.se'))))




####################################################################

#DIVERISTY AND EVENNESS OF ADULT HERRING SIZE
setwd(datadir)
#rvl<-read.csv("herring_lengths_RV_survey_spera_spawnar.csv",header=TRUE)
rvl<-read.csv("herring_lengths_RV_survey_spera_allar.csv",header=TRUE)
rvl$strat<-as.numeric(as.character(rvl$strat))
rvl<-subset(rvl,strat>=480 & strat<=495)
rvl$flen<-ifelse(rvl$mission=='NED2016016',rvl$flen/10,rvl$flen)
rvl$flen<-ifelse(rvl$mission!='NED2016016',(rvl$flen*1.0866)+0.95632,rvl$flen)
plot(log10(rvl$flen),log10(rvl$fwt))
rvl<-subset(rvl,log10(rvl$flen)>=.5 & month %in% c(6,7,8))#REMOVE OUTLIERS
rvl$flencat<-(floor(rvl$flen))+.5
rvl$sid<-gsub(' ','',paste(rvl$lon,'_',rvl$lat,'_',rvl$mission,'_',rvl$setno))

richfun<-function(d){
d<-subset(d,is.na(clen)==FALSE)
dout<-unique(subset(d,select=c('year','month','day','lon','lat')))

#RICHNESS WITH 1% CUTOFF LOCALLY
tcells<-sum(d$clen)#TOTAL NUMBER OF INDIVIDUALS

#SHANNON I
f2<-function(d3){ return(data.frame(p=sum(d3$clen)/tcells)) }
sdat<-ddply(subset(d,select=c('flencat','clen')),.(flencat),.fun=f2)
sdat$lp<-log(sdat$p)
sdat$pp<--1*(sdat$p*sdat$lp)
sdat<-subset(sdat,is.na(pp)==FALSE)

dout$shannon<-sum(sdat$pp,na.rm=TRUE)
dout$pe<-sum(sdat$pp)/(log(length(unique(d$flencat))))
dout$srange<-max(d$flen)-min(d$flen)
return(dout)
}
her.sdiv<-ddply(rvl,.(sid),.fun=richfun,.progress='text')


#HERRING AVERAGE SIZE DIVERSITY AND EVENNESS
mod<-gam(shannon~as.factor(year) + s(lon,lat,k=100),data=her.sdiv,gamma=1)
mod2<-gam(pe~as.factor(year) + s(lon,lat,k=100),data=her.sdiv,gamma=1)
pdat<-data.frame(year=sort(unique(her.sdiv$year)),
                 lon=gblon,
                 lat=gblat)
p<-predict(mod,newdata=pdat,se.fit=TRUE)
pdat$her.szdiv.rv<-p$fit
pdat$her.szdiv.rv.se<-p$se.fit

pdat2<-data.frame(year=sort(unique(subset(her.sdiv,is.na(pe)==FALSE)$year)),
lon=gblon,
lat=gblat)
p2<-predict(mod2,newdata=pdat2,se.fit=TRUE)
pdat2$her.szpe.rv<-p2$fit
pdat2$her.szpe.rv.se<-p2$se.fit

pdat<-merge(pdat,pdat2,by=c('year'),all=TRUE)

#ACCUMULATOR
andata<-c(andata,list(subset(pdat,select=c('year','her.szdiv.rv','her.szdiv.rv.se','her.szpe.rv','her.szpe.rv.se'))))



##############################################################
#ESTIMATE CHANGES IN DEPTH OF OCCUPANCY - NONE
setwd(datadir)
rvw<-read.csv("herring_weights_RV_survey_spera_allar.csv",header=TRUE)
rvw<-subset(rvw,month %in% c(6,7,8) & is.na(dmin)==FALSE)
rvw$depth<-(rvw$dmin+rvw$dmax)/2
rvw$strat<-as.numeric(as.character(rvw$strat))
rvw<-subset(rvw,strat>=480 & strat<=495)

f<-function(d){
    if((max(d$time)-min(d$time))>15){
    mod<-gam(totwgt~s(time,bs='cc',k=4)+s(lon,lat,k=20),data=d,gamma=1,family='nb'(link='log'))
s<-summary(mod)
    pdat<-data.frame(time=seq(min(d$time),max(d$time),1),
                     lon=median(d$lon),
                     lat=median(d$lat))
    pdat$p<-predict(mod,newdata=pdat,type='response')
    pdat$pz<-(pdat$p-mean(pdat$p,na.rm=TRUE))/sd(pdat$p,na.rm=TRUE)
return(data.frame(her.durng.rv=max(pdat$pz,na.rm=TRUE)[1]-min(pdat$pz,na.rm=TRUE)[1],her.dumx.rv=subset(pdat,p==max(pdat$p)[1])$time[1]))
} else NULL
}
ot<-ddply(rvw,.(year),.fun=f,.progress='text')

#ACCUMULATOR
andata<-c(andata,list(ot))


#HERRING AVERAGE LENGTH FROM RV SURVEY
#mod<-gam(depth~as.factor(year) + s(lon,lat,k=100) + s(time,bs='cc',k=5),weights=rvw$totno,data=rvw,gamma=1)
mod<-gam(depth~as.factor(year) + s(time,bs='cc',k=5),weights=rvw$totwt,data=rvw,gamma=1)
hdep<-data.frame(year=sort(unique(rvw$year)),
                 lon=gblon,
                 lat=gblat,
                 time=1200)
p<-predict(mod,newdata=hdep,se.fit=TRUE)
hdep$her.dep.rv<-p$fit
hdep$her.dep.rv.se<-p$se.fit
#plot(hdep$year,hdep$her.dep.rv,pch=15,ylim=c(60,120))
f<-function(d){    lines(c(d$year,d$year),c(d$her.dep.rv+(1.96*d$her.dep.rv.se),d$her.dep.rv-(1.96*d$her.dep.rv.se)))}
#ll<-dlply(hdep,.(year),.fun=f)


#ACCUMULATOR
andata<-c(andata,list(subset(hdep,select=c('year','her.dep.rv','her.dep.rv.se'))))


############################################################

#CALCULATE MASS-SPECIFIC METABOLIC RATE
#rvw<-read.csv("herring_weights_RV_survey_spera_allar.csv",header=TRUE)
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw<-subset(rvw,month %in% c(6,7,8) & is.na(dmin)==FALSE & rvw$totwgt>0)
rvw$strat<-as.numeric(as.character(rvw$strat))
rvw<-subset(rvw,strat>=480 & strat<=495)
rvw$depth<-(rvw$dmin+rvw$dmax)/2
rvw$totno<-ifelse(rvw$totno==0,1,rvw$totno)
rvw$M<-rvw$totwgt/rvw$totno
rvw$temperature<-(rvw$surface_temperature+rvw$bottom_temperature)/2
rvw<-subset(rvw,is.na(temperature)==FALSE)
rvw$kelvins<-rvw$temperature+273.15
k<-0.00008617#BOLTZMANN CONSTANT
rvw$metai<-(rvw$M^.25)*(exp(-1*(1/(k*rvw$kelvins))))
mod<-gam(metai~as.factor(year) + s(lon,lat,k=100),data=rvw,gamma=1)
met<-data.frame(year=sort(unique(rvw$year)),
                 lon=gblon,
                 lat=gblat)
p<-predict(mod,newdata=met,se.fit=TRUE)
met$her.metai.rv<-p$fit
met$her.metai.rv.se<-p$se.fit


#ACCUMULATOR
andata<-c(andata,list(subset(met,select=c('year','her.metai.rv','her.metai.rv.se'))))





##################################################################

#HADDOCK WEIGHT AND ABUNDANCE
setwd(datadir)
rvw<-read.csv("haddock_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw<-subset(rvw,month %in% c(6,7,8) & is.na(dmin)==FALSE)
rvw$strat<-as.numeric(as.character(rvw$strat))
rvw<-subset(rvw,strat>=480 & strat<=495)

#rvw$totno<-log10(rvw$totno+1)
#rvw$totwgt<-log10(rvw$totwgt+1)
mod<-gam(totwgt~as.factor(year)+s(lon,lat,k=100) + s(time,k=5,bs='cc'),data=rvw,gamma=1.4,family=nb)
mod2<-gam(totno~as.factor(year)+s(lon,lat,k=100) + s(time,k=5,bs='cc'),data=rvw,gamma=1.4,family=nb)
pdat<-data.frame(year=sort(unique(rvw$year)),
                 lon=gblon,
                 lat=gblat,
                 time=1200)
p<-predict(mod,newdata=pdat,se.fit=TRUE,type='response')
p2<-predict(mod2,newdata=pdat,se.fit=TRUE,type='response')
pdat$had.totwgt.rv<-p$fit
pdat$had.totwgt.rv.se<-p$se.fit
pdat$had.totno.rv<-p2$fit
pdat$had.totno.rv.se<-p2$se.fit


#HADDOCK DEPTH OF OCCUPANCY
rvw<-read.csv("haddock_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw$depth<-(rvw$dmax+rvw$dmin)/2
rvw<-subset(rvw,month %in% c(6,7,8) & is.na(dmin)==FALSE & rvw$totwgt>0)
rvw$strat<-as.numeric(as.character(rvw$strat))
rvw<-subset(rvw,strat>=480 & strat<=495)
mod<-gam(depth~as.factor(year)+s(lon,lat,k=100) + s(time,k=5,bs='cc'),weights=rvw$totno,data=rvw,gamma=1.4)
mod<-gam(depth~as.factor(year)+ s(time,k=5,bs='cc'),weights=rvw$totno,data=rvw,gamma=1.4)
p<-predict(mod,newdata=pdat,se.fit=TRUE)
pdat$had.dep.rv<-p$fit
pdat$had.dep.rv.se<-p$se.fit


rvl<-read.csv("haddock_lengths_RV_survey_spera_spawnar.csv",header=TRUE)
rvl<-subset(rvl,month %in% c(6,7,8))
rvl$strat<-as.numeric(as.character(rvl$strat))
rvl<-subset(rvl,strat>=480 & strat<=495)
mod<-gam(flen~as.factor(year)+s(lon,lat,k=100),weights=rvl$clen,data=rvl,gamma=1.4)
p<-predict(mod,newdata=pdat,se.fit=TRUE)
pdat$had.len.rv<-p$fit
pdat$had.len.rv.se<-p$se.fit

#ACCUMULATOR
andata<-c(andata,list(subset(pdat,select=c('year','had.totwgt.rv','had.totno.rv','had.dep.rv','had.len.rv','had.totwgt.rv.se','had.totno.rv.se','had.dep.rv.se','had.len.rv.se'))))



####################################################

#HADDOCK PREDATION INTENSITY
setwd(datadir)
load('haddock_predation_index.RData')
#ACCUMULATOR
andata<-c(andata,list(dat))



############################################################

#CALCULATE TREND IN AVERAGE WEIGHT OF HERRING FROM RV
#rvw<-read.csv("herring_weights_RV_survey_spera_allar.csv",header=TRUE)
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw<-subset(rvw,month %in% c(6,7,8) & is.na(dmin)==FALSE & rvw$totwgt>0)
rvw$strat<-as.numeric(as.character(rvw$strat))
rvw<-subset(rvw,strat>=480 & strat<=495)
rvw$totno<-ifelse(rvw$totno==0,1,rvw$totno)
rvw$M<-rvw$totwgt/rvw$totno
mod<-gam(M~as.factor(year) + s(lon,lat,k=100),data=rvw,gamma=1)
hmass<-data.frame(year=sort(unique(rvw$year)),
                 lon=gblon,
                 lat=gblat)
p<-predict(mod,newdata=hmass,se.fit=TRUE)
hmass$her.fmass.rv<-p$fit
hmass$her.fmass.rv.se<-p$se.fit

#ACCUMULATOR
andata<-c(andata,list(subset(hmass,select=c('year','her.fmass.rv','her.fmass.rv.se'))))





#########################################################

#HERRING AVERAGE ABUNDANCE FROM RV SURVEY
#rvw<-read.csv("herring_weights_RV_survey_spera_allar.csv",header=TRUE)
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
#rvw<-subset(rvw,month %in% c(6,7,8) & is.na(dmin)==FALSE & rvw$totwgt>0)
rvw<-subset(rvw,month %in% c(6,7,8) & is.na(dmin)==FALSE)
rvw$strat<-as.numeric(as.character(rvw$strat))
rvw<-subset(rvw,strat>=480 & strat<=495)
rvw$no<-log10(rvw$totno+1)
rvw$wgt<-log10(rvw$totwgt+1)
rvw$id<-gsub(' ','',paste(rvw$mission,'_',rvw$setno))

dm<-data.frame(year=sort(unique(rvw$year)),
               ntows=tapply(rvw$id,rvw$year,function(x) length(unique(x))))
rvw<-merge(rvw,dm,by=c('year'),all.x=TRUE,all.y=FALSE)

mod<-gam(totno~as.factor(year) + s(lon,lat,k=100) + s(time,k=5,bs='cc'),data=rvw,gamma=1,family=nb)
mod2<-gam(totwgt~as.factor(year) + s(lon,lat,k=100) + s(time,k=5,bs='cc'),data=rvw,gamma=1,family=nb)
hno<-data.frame(year=sort(unique(rvw$year)),
                 lon=gblon,
                 lat=gblat,
                time=1200)
p<-predict(mod,newdata=hno,se.fit=TRUE,type='response')
p2<-predict(mod2,newdata=hno,se.fit=TRUE,type='response')
hno$her.totno.rv<-p$fit
hno$her.totno.rv.se<-p$se.fit
hno$her.totwgt.rv<-p2$fit
hno$her.totwgt.rv.se<-p2$se.fit


#ACCUMULATOR
andata<-c(andata,list(subset(hno,select=c('year','her.totwgt.rv','her.totwgt.rv.se','her.totno.rv','her.totno.rv.se'))))



##########################################################
########################################################
#CALCULATES DETAILED PHENOLOGY CHANGES FOR SST/WIND
#LOADS PATHFINDER
setwd(datadir)
load("pathfinder_sst_daily_1981_2012_spawnar.RData")
dat$cell<-gsub(' ','',paste(dat$lon,'_',dat$lat))


#ESTIMATES PHENOLOGY CHANGE AVERAGED OVER CORE AREA
f<-function(d){
print(unique(d$year))
if(length(unique(d$month))>=10){
    mod<-gam(sst~s(day,bs='cc',k=7) + s(lon,lat,k=20),data=d,gamma=1)
pdat<-data.frame(lon=-66.3,
                 lat=43.3,
                 day=seq(1,365,1),
                 year=mean(d$year))
 p<-predict(mod,newdata=pdat,se.fit=TRUE)
 pdat$p<-p$fit
 pdat$se<-p$se.fit

#TEMPERATURES ABOVE 18
pdat0<-subset(pdat,p>=12.3)
dur12<-dim(pdat0)[1]
tim12<-min(pdat0$day)

#LOOK AT FALL ONLY
pdatfall<-subset(pdat,day>200 & day<300)
mnfall<-mean(pdatfall$p)
mxfall<-max(pdatfall$p)
pdatfall2<-subset(d,day>200 & day<300)
sdfall<-sd(pdatfall2$sst)

mx<-subset(pdat,p==max(pdat$p))
mn<-subset(pdat,p==min(pdat$p))

out<-data.frame(sst.tmax=mx$day,
                sst.tmin=mn$day,
                sst.amp=diff(range(pdat$p)),
                sst.max=mx$p,
                sst.min=mn$p,
                sst.dur12=dur12,
                sst.t12=tim12,
                sst.fmin=mnfall,
                sst.fmax=mxfall,
                sst.fsd=sdfall)
 return(out)
} else NULL
}
sstpath<-ddply(subset(dat,year>1981),.(year),.fun=f,.progress='text')



#ESTIMATES PHENOLOGY CHANGE AVERAGED OVER CORE AREA
f<-function(d){
print(unique(d$year))
if(length(unique(d$month))>=10){
    mod<-gam(wind~s(day,bs='cc',k=7) + s(lon,lat,k=20),data=d,gamma=1)
pdat<-data.frame(lon=-66.3,
                 lat=43.3,
                 day=seq(1,365,1),
                 year=mean(d$year))
 p<-predict(mod,newdata=pdat,se.fit=TRUE)
 pdat$p<-p$fit
 pdat$se<-p$se.fit

#LOOK AT FALL ONLY
pdatfall<-subset(pdat,day>200 & day<300)
mnfall<-mean(pdatfall$p)
mxfall<-max(pdatfall$p)
pdatfall2<-subset(d,day>200 & day<300)
sdfall<-sd(pdatfall2$sst)

mx<-subset(pdat,p==max(pdat$p))
mn<-subset(pdat,p==min(pdat$p))

wind.tmax<-subset(mx,select=c('day'))
wind.tmax$day<-wind.tmax$day+180
wind.tmax$day<-ifelse(wind.tmax$day>365,wind.tmax$day-365,wind.tmax$day)

out<-data.frame(wind.tmax=wind.tmax$day[1],
                wind.tmin=mn$day[1],
                wind.amp=diff(range(pdat$p))[1],
                wind.max=mx$p[1],
                wind.min=mn$p[1],
                wind.fmin=mnfall[1],
                wind.fmax=mxfall[1],
                wind.fsd=sdfall[1])
 return(out)
} else NULL
}
windpath<-ddply(subset(dat,year>1981),.(year),.fun=f,.progress='text')


######################################################################
#CALCULATES PROPORTION OF AREA THAT QUALIFIES AS SUITABLE THERMAL HABITAT (SST<12.3) DURING THE FALL (DAY 200-300)
f<-function(d){
if(length(unique(d$month))>=10){
mod<-gam(sst~s(day,bs='cc',k=6),data=d,gamma=1.4)
pdat<-data.frame(day=seq(1,365,1))

p<-predict(mod,newdata=data.frame(day=pdat$day),se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat<-subset(pdat,day>=200 & day<=300)

pdat$lon<-unique(d$lon)
pdat$lat<-unique(d$lat)
return(pdat)
} else NULL
}
sst.path2<-ddply(dat,.(cell,year),.fun=f,.progress='text')

f2<-function(d){
d2<-subset(d,p<=12.3)
if(dim(d2)[1]<=0){
out<-data.frame(sst.sha=0)
} else {
out<-data.frame(sst.sha=(dim(d2)[1]/dim(d)[1])*100)
}
return(out)
}
sst.path3<-ddply(sst.path2,.(year,day),.fun=f2,.progress='text')

f3<-function(d){
    return(data.frame(sst.sha=mean(d$sst.sha),
                      sst.sha.se=sd(d$sst.sha)))
}
sst.path4<-ddply(sst.path3,.(year),.fun=f3)

sstpath<-merge(sstpath,sst.path4,by=c('year'),all=TRUE)
sstpath<-merge(sstpath,windpath,by=c('year'),all=TRUE)

#ACCUMULATOR
andata<-c(andata,list(sstpath))


#CENTRAL TENDANCY PHENOLOGY SERIES
setwd(datadir)
load('ct_data.Rdata')
andata<-c(andata,list(ctdat))



#############################################################
#SST FROM AZMP STATIONS: HALIFAX, PRINCE, STANDREWS, LURCHER, BOOTHBAY, GEORGES
setwd(datadir)
load('azmp_temp_series.RData')
andata<-c(andata,list(a))
#data<-merge(data,a,by=c('year'),all=TRUE)


#JELLYFISH OBS FROM BOFLARVAL
setwd(datadir)
load("BoF_larval_herring_counts_spera_spawnar.RData")
larv<-subset(dat,bathy<0 & !(spec_stage %in% c('A','E','J','Y','I')))
larv$stdno<-(larv$totno/larv$volm3)*abs(larv$bathy)#
larv<-subset(larv,is.na(stdno)==FALSE)
larv<-subset(larv,gearprocdesc=='Sawtooth Oblique')

a<-subset(larv,phylum=='Cnidaria' & is.na(stdno)==FALSE)
mod<-gam(stdno~as.factor(year) + s(lon,lat,k=20) + s(day,k=5,bs='cc'),family=gaussian,data=a)
pdat<-data.frame(year=sort(unique(a$year)),
                 lon=gblon,
                 lat=gblat,
                 day=250)
p<-predict(mod,newdata=pdat,type='response',se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit

jelly<-subset(pdat,select=c('year','p','se'))
names(jelly)<-c('year','her.jf.bof','her.jf.bof.se')

andata<-c(andata,list(jelly))




#DEPTH CHANGES FOR LARVAL HERRING
setwd(datadir)
load("BoF_larval_herring_counts_spera_spawnar.RData")
larv<-subset(dat,bathy<0 & spec_stage %in% c('L') & month>8 & gearprocdesc=='Sawtooth Oblique')
#STANDARDIZE PER VOLUME STRAINED
larv$stdno<-(larv$totno/larv$volm3)#NUMBER PER M2 - INTEGRATED OVER DEPTH
larv<-subset(larv,is.na(stdno)==FALSE )
larv$phylum<-tolower(larv$phylum)
larv$order<-tolower(larv$order)
larv$kingdom<-tolower(larv$kingdom)
larv$family<-tolower(larv$family)
#ADD ID FOR SAMPLES
larv$sid<-gsub(' ','',paste(larv$cruise,'_',larv$setno,'_',larv$sample))

larv$depth<-(larv$startdepth+larv$enddepth)/2
larv<-subset(larv,stdno>0)
mod<-gam(depth~as.factor(year) + s(lon,lat,k=20) + s(day,k=4,bs='cc'),weights=larv$stdno,data=larv)
lrvdep<-data.frame(year=sort(unique(larv$year)),
                 lon=gblon,
                 lat=gblat,
                 day=250)
p<-predict(mod,newdata=lrvdep,type='response',se.fit=TRUE)
lrvdep$herlrv.dep.bof<-p$fit
lrvdep$herlrv.dep.bof.se<-p$se.fit
lrvdep<-subset(lrvdep,select=c('year','herlrv.dep.bof','herlrv.dep.bof.se'))

#ACCUMULATOR
andata<-c(andata,list(lrvdep))




###################################
#COMBINE ALL INTO A SINGLE DATABASE
data<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), andata)#COMBINE


#VARIABLES ARE COLLINEAR WITH OTHERS - OMIT
#data<-data[,!(names(data) %in% c('zp.cop.cpr','herlrv.sp.n','herlrv.georng.n','sst.max','sst.fmin','wind.fsd','zp.cop.cpr.se','her.georng.n','her.georng','her.sp.n','sst.min','wind.min','herlrv.spvar'))]




###########################################################

##CALCULATE HERRING 'STATE' VARIABLE

########################################################

#SAME AS ABOVE BUT EXCLUDES VPA SSB SERIES
d<-subset(data,is.na(year)==FALSE,select=c('year','her.cf.rv','her.len.rv','her.spcv','her.waa','her.metai.rv','her.fmass.rv','her.szpe.rv','her.spvar'))
d$her.spcv<-d$her.spcv*-1
d$her.spvar<-d$her.spvar*-1

#STANDARDIZE COLUMNS TO UNIT VARIANCE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
d<-data.frame(cbind(subset(d,select='year'),apply(d[,2:dim(d)[2]],2,f)))

fmn<-function(x){mean(x,na.rm=TRUE)}
fsd<-function(x){sd(x,na.rm=TRUE)}
dd<-data.frame(cbind(subset(d,select=c('year')),apply(d[,2:dim(d)[2]],1,fmn),apply(d[,2:dim(d)[2]],1,fsd)))
names(dd)<-c('year','her.state','her.state.se')
d<-merge(d,dd,by=c('year'),all=TRUE)

cr<-cor(d,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
cr.t<-subset(cr.t,!(Var2 %in% c('year','her.state.se')) & !(Var1 %in% c('year','her.state.se')))

#AVERAGE CORRELATION OF 'STATE' VARIABLES
cc<-subset(cr.t,Var1!='her.state' & Var2!='her.state')
#mean(cc$Freq);
#mdr<-round(median(cc$Freq),digits=2)
#MEAN=0.54; MEDIAN=0.56


setwd(figsdir)
pdf('herring_state_derivation.pdf',height=9,width=11)
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(b$year,b$her.ssb,type='l',xlim=c(1965,2015),lwd=2,las=1,xlab='Year',ylab='Z-score',col='magenta')
lines(b$year,b$her.cf.rv,type='l',col='green',lwd=2)
lines(b$year,b$her.len.rv,type='l',col='lightblue',lwd=2)
lines(b$year,b$her.spcv,type='l',col='gold',lwd=2)
lines(b$year,b$her.waa,type='l',col='orange',lwd=2)
lines(b$year,b$her.mass.rv,type='l',col='darkblue',lwd=2)
lines(b$year,b$her.metai.rv,type='l',col='yellow',lwd=2)
lines(b$year,b$her.szpe,type='l',col='pink',lwd=2)
lines(b$year,b$her.spvar,type='l',col='black',lwd=2)
lines(b$year,b$her.state,type='l',col='red',lwd=3)
legend('topright',legend=paste('r = ',mdr),bty='n')
legend(1998,2.4,legend=names(b)[2:11],bty='n',col=c('magenta','green','lightblue','gold','orange','darkblue','yellow','pink','black','red'),lwd=2)

b$upr99<-b$her.state+(b$her.state.se)
b$lwr99<-b$her.state-(b$her.state.se)
plot(b$year,b$her.state,type='l',pch=15,ylim=c(-1.5,1.75),las=1,xlab='Year',ylab='Herring state',xlim=c(1965,2016))
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr99,b$lwr99[length(b$lwr99):1]),col=alpha('lightskyblue1',.75),border=NA)
lines(b$year,b$her.state,ylim=c(-1.5,1.75),xlim=c(1965,2016),col='dodgerblue4')

d$upr99<-d$her.state+(d$her.state.se)
d$lwr99<-d$her.state-(d$her.state.se)
plot(d$year,d$her.state,type='l',pch=15,ylim=c(-1.5,1.75),las=1,xlab='Year',ylab='Herring state (excluding VPA)',xlim=c(1965,2016))
polygon(c(d$year,d$year[length(d$year):1]),c(d$upr99,d$lwr99[length(d$lwr99):1]),col=alpha('lightskyblue1',.75),border=NA)
lines(d$year,d$her.state,ylim=c(-1.5,1.75),xlim=c(1965,2016),col='dodgerblue4')
dev.off()


########################################################
#ADD TO DATA
b<-subset(b,select=c('year','her.state','her.state.se'))
data<-merge(data,b,by=c('year'),all.x=TRUE,all.y=FALSE)
data<-subset(data,is.na(year)==FALSE & year>=1965)

#CREATE LAGGED VERSIONS OF HERRING STATE
data<-slide(data,Var='her.state',slideBy=-3,NewVar='her.state.tplus3')
data<-slide(data,Var='her.state',slideBy=-4,NewVar='her.state.tplus4')





#DERIVATION OF SST STATE VARIABLE
b<-(subset(data,select=c('year','sst.sha','sst.fmax','sst.t12','sst.dur12','sst.amp','sst.mn','t50','sst.stalt','sst.grg','sst.prin','sst.lur')))
b$sst.t12<-ifelse(b$sst.t12=='Inf',NA,b$sst.t12)
b$sst.sha<-b$sst.sha*-1
b$sst.t12<-b$sst.t12*-1
plot(b,pch=16)

#STANDARDIZE COLUMNS TO UNIT VARIANCE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
b<-data.frame(cbind(subset(b,select='year'),apply(b[,2:dim(b)[2]],2,f)))

fmn<-function(x){mean(x,na.rm=TRUE)}
fsd<-function(x){sd(x,na.rm=TRUE)}
bb<-data.frame(cbind(subset(b,select=c('year')),apply(b[,2:dim(b)[2]],1,fmn),apply(b[,2:dim(b)[2]],1,fsd)))
names(bb)<-c('year','sst.state','sst.state.se')
b<-merge(b,bb,by=c('year'),all=TRUE)


cr<-cor(b,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
cr.t<-subset(cr.t,!(Var2 %in% c('year','sst.state.se')) & !(Var1 %in% c('year','sst.state.se')))
cr.t<-cr.t[order(cr.t$Freq,decreasing=TRUE),]

#AVERAGE CORRELATION OF 'STATE' VARIABLES
cc<-subset(cr.t,Var1!='sst.state' & Var2!='sst.state')
#mean(cc$Freq);
mdr<-round(median(cc$Freq),digits=2)
#MEAN=0.7; MEDIAN=0.76

setwd(figsdir)
pdf('sst_state_derivation.pdf',height=9,width=11)
par(mfrow=c(2,2),mar=c(4,4,1,1))
cls<-c('black','green','blue','pink','gold','deeppink','gray','orange','lightblue','gold3','darkred','aquamarine')

plot(b$year,b$sst.sha,type='l',xlim=c(1920,2015),lwd=3,las=1,xlab='Year',ylab='Z-score',ylim=c(-3,3.5),col=alpha(cls[1],.6))
lines(b$year,b$sst.fmax,type='l',col=alpha(cls[2],.6),lwd=3)
lines(b$year,b$sst.t12,type='l',col=alpha(cls[3],.6),lwd=3)
lines(b$year,b$sst.dur12,type='l',col=alpha(cls[4],.6),lwd=3)
lines(b$year,b$sst.amp,type='l',col=alpha(cls[5],.6),lwd=3)
lines(b$year,b$sst.mn,type='l',col=alpha(cls[6],.6),lwd=3)
lines(b$year,b$t50,type='l',col=alpha(cls[7],.6),lwd=3)
lines(b$year,b$sst.stalt,col=alpha(cls[8],.6),lwd=2)
lines(b$year,b$sst.grg,col=alpha(cls[9],.6),lwd=2)
lines(b$year,b$sst.prin,col=alpha(cls[10],.6),lwd=2)
lines(b$year,b$sst.bbay,col=alpha(cls[11],.6),lwd=2)
lines(b$year,b$sst.lur,col=alpha(cls[12],.6),lwd=2)
lines(b$year,b$sst.state,type='l',col='red',lwd=3)
legend('topright',legend=paste('r = ',mdr),bty='n')
legend('topleft',legend=names(b)[2:12],col=cls,lwd=2,bty='n',ncol=2)

b<-b[order(b$year,decreasing=FALSE),]
b$sst.state.se<-ifelse(is.na(b$sst.state.se)==TRUE,0,b$sst.state.se)
b$upr99<-b$sst.state+(b$sst.state.se)
b$lwr99<-b$sst.state-(b$sst.state.se)
plot(b$year,b$sst.state,type='l',pch=15,ylim=c(-2.5,2.5),las=1,xlab='Year',ylab='SST state',xlim=c(1920,2016))
abline(h=0,lty=3,col='gray')
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr99,b$lwr99[length(b$lwr99):1]),col=alpha('lightskyblue1',.75),border=NA)
lines(b$year,b$sst.state,type='l',col='dodgerblue4')
dev.off()

b<-subset(b,select=c('year','sst.state','sst.state.se'))
data<-merge(data,b,by=c('year'),all=TRUE)





########################################################

#DERIVATION OF NUTRIENT STATE VARIABLE
b<-(subset(data,select=c('year','nit','sil','phos')))
plot(b,pch=16)

#STANDARDIZE COLUMNS TO UNIT VARIANCE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
b<-data.frame(cbind(subset(b,select='year'),apply(b[,2:dim(b)[2]],2,f)))

#THIS ESTIMATES BY MIXED MODEL - ALMOST SAME AS JUST AVERAGING
bb<-b %>% gather(var,y,nit:phos)
bb<-na.omit(bb)
mod<-gamm(y~as.factor(year),data=bb,random=list(var=~1))
pdat<-data.frame(year=sort(unique(bb$year)))
p<-predict(mod$gam,newdata=pdat,se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit



fmn<-function(x){mean(x,na.rm=TRUE)}
fsd<-function(x){sd(x,na.rm=TRUE)}
bb<-data.frame(cbind(subset(b,select=c('year')),apply(b[,2:dim(b)[2]],1,fmn),apply(b[,2:dim(b)[2]],1,fsd)))
names(bb)<-c('year','nut.state','nut.state.se')
b<-merge(b,bb,by=c('year'),all=TRUE)


cr<-cor(b,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
cr.t<-subset(cr.t,!(Var2 %in% c('year','nut.state.se')) & !(Var1 %in% c('year','nut.state.se')))
cr.t<-cr.t[order(cr.t$Freq,decreasing=TRUE),]

#AVERAGE CORRELATION OF 'STATE' VARIABLES
cc<-subset(cr.t,Var1!='nut.state' & Var2!='nut.state')
#mean(cc$Freq);
mdr<-round(median(cc$Freq),digits=2)
#MEAN=0.75; MEDIAN=0.83

setwd(figsdir)
pdf('nut_state_derivation.pdf',height=9,width=11)
par(mfrow=c(2,2),mar=c(4,4,1,1))
cls<-c('green','blue','deeppink','red3')

plot(b$year,b$nit,type='l',xlim=c(1965,2015),lwd=3,las=1,xlab='Year',ylab='Z-score',ylim=c(-2,3.5),col=alpha(cls[1],.6))
lines(b$year,b$sil,type='l',col=alpha(cls[2],.6),lwd=3)
lines(b$year,b$phos,type='l',col=alpha(cls[3],.6),lwd=3)
lines(b$year,b$nut.state,type='l',col='red',lwd=3)
legend('topright',legend=paste('r = ',mdr),bty='n')
legend('topleft',legend=names(b)[2:5],col=cls,lwd=2,bty='n',ncol=1)

b<-b[order(b$year,decreasing=FALSE),]
b$nut.state.se<-ifelse(is.na(b$nut.state.se)==TRUE,0,b$nut.state.se)
b$upr99<-b$nut.state+(b$nut.state.se)
b$lwr99<-b$nut.state-(b$nut.state.se)
b$upr99<-ifelse(is.na(b$upr99)==TRUE,0,b$upr99)
b$lwr99<-ifelse(is.na(b$lwr99)==TRUE,0,b$lwr99)
plot(b$year,b$nut.state,type='l',pch=15,ylim=c(-2.5,2.5),las=1,xlab='Year',ylab='NUT state',xlim=c(1965,2016))
abline(h=0,lty=3,col='gray')
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr99,b$lwr99[length(b$lwr99):1]),col=alpha('lightskyblue1',.75),border=NA)
lines(b$year,b$nut.state,type='l',col='dodgerblue4')
dev.off()

b<-subset(b,select=c('year','nut.state','nut.state.se'))
data<-merge(data,b,by=c('year'),all=TRUE)






########################################################

#DERIVATION OF CENTRAL TENDANCY STATE VARIABLE
b<-(subset(data,select=c('year','sst.t12','chl.ct','strt.ct','t50.ct','wind.tmin')))
b$sst.t12<-ifelse(b$sst.t12==Inf,NA,b$sst.t12)
plot(b,pch=16)

#STANDARDIZE COLUMNS TO UNIT VARIANCE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
b<-data.frame(cbind(subset(b,select='year'),apply(b[,2:dim(b)[2]],2,f)))

fmn<-function(x){mean(x,na.rm=TRUE)}
fsd<-function(x){sd(x,na.rm=TRUE)}
bb<-data.frame(cbind(subset(b,select=c('year')),apply(b[,2:dim(b)[2]],1,fmn),apply(b[,2:dim(b)[2]],1,fsd)))
names(bb)<-c('year','ct.state','ct.state.se')
b<-merge(b,bb,by=c('year'),all=TRUE)


cr<-cor(b,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
cr.t<-subset(cr.t,!(Var2 %in% c('year','ct.state.se')) & !(Var1 %in% c('year','ct.state.se')))
cr.t<-cr.t[order(cr.t$Freq,decreasing=TRUE),]

#AVERAGE CORRELATION OF 'STATE' VARIABLES
cc<-subset(cr.t,Var1!='ct.state' & Var2!='ct.state')
#mean(cc$Freq);
mdr<-round(median(cc$Freq,na.rm=TRUE),digits=2)
#MEAN=0.75; MEDIAN=0.83

setwd(figsdir)
pdf('ct_state_derivation.pdf',height=9,width=11)
par(mfrow=c(2,2),mar=c(4,4,1,1))
cls<-c('green','blue','deeppink','gold3','black','brown','red3')

plot(b$year,b$strt.ct,type='l',xlim=c(1965,2015),lwd=3,las=1,xlab='Year',ylab='Z-score',ylim=c(-2,3.5),col=alpha(cls[1],.6))
lines(b$year,b$chl.ct,type='l',col=alpha(cls[2],.6),lwd=3)
lines(b$year,b$t50.ct,type='l',col=alpha(cls[3],.6),lwd=3)
lines(b$year,b$wind.tmin,type='l',col=alpha(cls[4],.6),lwd=3)
lines(b$year,b$sst.t12,type='l',col=alpha(cls[5],.6),lwd=3)
lines(b$year,b$nut.ct,type='l',col=alpha(cls[6],.6),lwd=3)
lines(b$year,b$ct.state,type='l',col='red',lwd=3)
legend('topright',legend=paste('r = ',mdr),bty='n')
legend('topleft',legend=names(b)[2:5],col=cls,lwd=2,bty='n',ncol=1)

b<-b[order(b$year,decreasing=FALSE),]
b$nut.state.se<-ifelse(is.na(b$nut.state.se)==TRUE,0,b$nut.state.se)
b$upr99<-b$nut.state+(b$nut.state.se)
b$lwr99<-b$nut.state-(b$nut.state.se)
b$upr99<-ifelse(is.na(b$upr99)==TRUE,0,b$upr99)
b$lwr99<-ifelse(is.na(b$lwr99)==TRUE,0,b$lwr99)
plot(b$year,b$nut.state,type='l',pch=15,ylim=c(-2.5,2.5),las=1,xlab='Year',ylab='NUT state',xlim=c(1965,2016))
abline(h=0,lty=3,col='gray')
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr99,b$lwr99[length(b$lwr99):1]),col=alpha('lightskyblue1',.75),border=NA)
lines(b$year,b$nut.state,type='l',col='dodgerblue4')
dev.off()

b<-subset(b,select=c('year','nut.state','nut.state.se'))
data<-merge(data,b,by=c('year'),all=TRUE)



















#CREATE LAGGED VERSIONS OF HERRING STATE
data<-slide(data,Var='sst.state',slideBy=-3,NewVar='sst.state.tplus3')
data<-slide(data,Var='sst.state',slideBy=-4,NewVar='sst.state.tplus4')

setwd(datadir)
data.spawnar<-data
save(data.spawnar,file='SPERA_andata_spawnar.RData')

herl<-read.csv('herring_historical_landings.csv',header=FALSE,col.names=c('year','her.totfishers','her.landvalue','her.land','her.valuepp'))
herl$year<-as.numeric(as.character(herl$year))
herl$her.totfishers<-as.numeric(as.character(herl$her.totfishers))
herl$her.landvalue<-as.numeric(as.character(herl$her.landvalue))
herl$her.land<-as.numeric(as.character(herl$her.land))
herl$her.valuepp<-as.numeric(as.character(herl$her.valuepp))
herl<-subset(herl,is.na(year)==FALSE)

plot(herl,pch=15)
plot(herl$her.valuepp,log10(herl$her.land),pch=15)
cor(herl$her.valuepp,log10(herl$her.land),use='pairwise.complete.obs')

plot(data$herlrv.len,log(data$herlrv.mn.bof),pch=15)

plot(data$sst.t12,data$her.state,pch=15)

plot(data$year,data$her.state,pch=15,type='b',col='red',yaxt='n')
par(new=TRUE)
plot(data$year,data$her.expr,pch=15,type='b')
plot(data$year,data$herlrv.len,pch=15,type='b')
plot(data$year,data$sst.state*-1,pch=15,type='b')
par(new=TRUE)
plot(data$year,data$sst.sha,pch=15,type='b',col='blue')
















#####################################################

#TRY TO ESTABLISH STATE VARIABLE FOR HERRING LARVAE


names(data)
aa<-subset(data,select=c('her.ssb','her.ssb.tplus3','her.ssb.tplus4','her.ssb.tplus5','her.rec1','herlrv.len','herlrv.mn.bof'))
aa$herlrv.mn.bof<-log10(aa$herlrv.mn.bof)
aa$her.rec1<-log10(aa$her.rec1+.1)
plot(aa,pch=15)



b<-na.omit(subset(data,select=c('year','herlrv.mn.bof','herlrv.spnug','herlrv.sprng','herlrv.spcv','herlrv.len','her.rec1','herlrv.surv','her.state.tplus3')))
b$herlrv.sprng<-log10(b$herlrv.sprng+.1)
b$herlrv.spnug<-log10(b$herlrv.spnug+1)
b$herlrv.mn.bof<-log10(b$herlrv.mn.bof)
b$her.rec1<-log10(b$her.rec1+1)
b$herlrv.surv<-log10(b$herlrv.surv+1)
#b$herlrv.sprng<-b$herlrv.sprng*-1
#b$herlrv.len<-b$herlrv.len*-1
plot(b,pch=15)


cor(data$herlrv.len,data$her.state.tplus3,use='pairwise.complete.obs',method='spearman')
cor(data$herlrv.len,data$herlrv.surv,use='pairwise.complete.obs',method='spearman')

b$herlrv.mn.bof<-scale(b$herlrv.mn.bof,center=TRUE,scale=TRUE)
b$herlrv.spnug<-scale(b$herlrv.spnug,center=TRUE,scale=TRUE)
b$herlrv.spvar<-scale(b$herlrv.spvar,center=TRUE,scale=TRUE)
b$herlrv.sprng<-scale(b$herlrv.sprng,center=TRUE,scale=TRUE)
b$herlrv.spcv<-scale(b$herlrv.spcv,center=TRUE,scale=TRUE)
b$herlrv.len<-scale(b$herlrv.len,center=TRUE,scale=TRUE)

f<-function(d){
d$herlrv.state<-mean(c(d$herlrv.mn.bof,d$herlrv.spnug,d$herlrv.spvar,d$herlrv.sprng,d$herlrv.spcv,d$herlrv.len),na.rm=TRUE)
d$herlrv.state.se<-sd(c(d$herlrv.mn.bof,d$herlrv.spnug,d$herlrv.spvar,d$herlrv.sprng,d$herlrv.spcv,d$herlrv.len),na.rm=TRUE)
    return(d)
}
b<-ddply(b,.(year),.fun=f)

cr<-cor(b,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
cr.t<-subset(cr.t,!(Var2 %in% c('year','herlrv.state.se')) & !(Var1 %in% c('year','herlrv.state.se')))


#AVERAGE CORRELATION OF 'STATE' VARIABLES
cc<-subset(cr.t,Var1!='herlrv.state' & Var2!='herlrv.state')
#mean(cc$Freq);
mdr<-round(median(cc$Freq),digits=2)
#MEAN=0.5; MEDIAN=0.44

setwd(figsdir)
pdf('herlarvae_state_derivation.pdf',height=10,width=8)
par(mfrow=c(2,1),mar=c(4,4,1,1))
dm<-data.frame(year=seq(min(b$year),max(b$year),1))
b<-merge(b,dm,by=c('year'),all=TRUE)
plot(b$year,b$herlrv.mn.bof,type='l',xlim=c(1975,2010),lwd=3,las=1,xlab='Year',ylab='Z-score')
lines(b$year,b$herlrv.spnug,type='l',col='green',lwd=3)
lines(b$year,b$herlrv.spvar,type='l',col='blue',lwd=3)
lines(b$year,b$herlrv.sprng,type='l',col='pink',lwd=3)
lines(b$year,b$herlrv.spcv,type='l',col='gold',lwd=3)
lines(b$year,b$herlrv.len,type='l',col='orange',lwd=3)
lines(b$year,b$herlrv.state,type='l',col='red',lwd=3)
legend('topright',legend=paste('r = ',mdr),bty='n')

plot(b$year,b$sst.state,type='l',pch=15,ylim=c(-1.5,2),las=1,xlab='Year',ylab='Sstring state',xlim=c(1980,2016))
lines(b$year,b$sst.state+(1*b$sst.state.se),lty=3,lwd=2)
lines(b$year,b$sst.state-(1*b$sst.state.se),lty=3,lwd=2)
dev.off()

b<-subset(b,select=c('year','sst.state','sst.state.se'))
data<-merge(data,b,by=c('year'),all.x=TRUE,all.y=FALSE)
data<-subset(data,is.na(year)==FALSE & year>=1965)




plot(b,pch=15)

plot(bb,pch=15)

subset(data,select=c('year','her.state','her.state.tplus3'))















nms<-names(data)
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.sta',nms)==FALSE]
nms<-nms[grepl('year',nms)==FALSE]
a<-data[,names(data)%in% nms]
a$herlrv.mn.bof<-log10(a$herlrv.mn.bof)

hist(a$herlrv.surv,breaks=15)
cor(a$her.ssb.tplus3,(a$herlrv.len),use='pairwise.complete.obs',method='spearman')

cr<-cor(a,use='pairwise.complete.obs',method='kendall')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
#cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
cutoff<-.25

aa<-subset(cr.t,abs(Freq)>.9)

f<-function(vr1,vrs2){
d<-subset(cr.t,Var1==vr1 & abs(Freq)>=cutoff & !(Var2 %in% vrs2))
d<-d[order(abs(d$Freq),decreasing=TRUE),]
print(d)
return(d)
}
cher4<-f('her.ssb.tplus4',c('her.ssb.tplus3','her.ssb.tplus4','her.ssb.tplus5','her.ssb.tplus6','her.ssb'))
cher4<-f('her.ssb',c('her.ssb.tplus3','her.ssb.tplus4','her.ssb.tplus5','her.ssb.tplus6'))
ccf<-f('her.cf.rv',c('her.ssb.tplus3','her.ssb.tplus4','her.ssb.tplus5','her.ssb.tplus6'))

clrvlen<-f('herlrv.len',c('her.ssb.tplus3','her.ssb.tplus4','her.ssb.tplus5','her.ssb.tplus6'))
crec<-f('her.rec1',c('her.ssb.tplus3','her.ssb.tplus4','her.ssb.tplus5','her.ssb.tplus6'))




a<-subset(data,select=c('her.cf.rv','her.ssb','her.ssb.tplus3','her.ssb.tplus4','her.ssb.tplus5'))
plot(a,pch=15)
round(cor(a,use='pairwise.complete.obs'),digits=2)

csurv<-subset(cr.t,Var1=='herlrv.surv' & abs(Freq)>=cutoff & !(Var2 %in% c('her.ssb.tplus3','her.ssb.tplus4','her.ssb.tplus5','her.ssb.tplus6')))
csurv<-csurv[order(abs(csurv$Freq),decreasing=TRUE),]


a<-subset(data,select=c('her.ssb.tplus3','her.ssb.tplus4','her.ssb.tplus5','her.ssb.tplus6','herlrv.surv','herlrv.mn.bof','herlrv.len'))
a$herlrv.mn.bof<-log10(a$herlrv.mn.bof)
plot(a,pch=15)
plot(data$herlrv.surv,data$her.ssb.tplus4)
plot(data$herlrv.mn.bof,data$herlrv.surv)

plot(data$herlrv.len,data$herlrv.surv)
plot(log10(data$herlrv.len),log10(data$herlrv.surv+1))
cor(data$herlrv.len,data$herlrv.surv,use='pairwise.complete.obs')
cor(log10(data$herlrv.len),log10(data$herlrv.surv+1),use='pairwise.complete.obs')


plot(a,pch=15)
round(cor(a,use='pairwise.complete.obs'),digits=2)
hist(data$herlrv.len,breaks=100)

plot(data$herlrv.len,data$her.ssb)

slide(data,Var='her.ssb',slideBy=-5))
plot(data$year,data$herlrv.len)

d<-subset(data,select=c('her.ssb.tmin4','herlrv.mn.bof','herlrv.len','her.len.rv.t1','her.f.t1'))
d<-na.omit(d)
mod<-stepAIC(lm(her.ssb.tmin4~.,data=d,na.action=na.omit))
mod<-stepAIC(lm(her.ssb-4~herlrv.mn.bof + herlrv.len + her.len.rv1 + her.cf.rv1 +her.f1,data=data,na.action=na.omit))









r
day<-subset(rvw,time>700 & time<2000)
ngt<-subset(rvw,time<=700 | time>=2000)
mod<-gam(totwgt~as.factor(year)+s(lon,lat,k=100) + s(time,k=4),data=day,gamma=1.4,family=nb)
pdat.d<-data.frame(year=sort(unique(rvw$year)),
                 lon=gblon,
                 lat=gblat,
                 time=1200)
p<-predict(mod,newdata=pdat.d,type='response')
pdat.d$pd<-p

mod<-gam(totwgt~as.factor(year)+s(lon,lat,k=100) + s(time,k=4),data=ngt,gamma=1.4,family=nb)
pdat.n<-data.frame(year=sort(unique(rvw$year)),
                 lon=gblon,
                 lat=gblat,
                   time=2300)
p<-predict(mod,newdata=pdat.n,type='response')
pdat.n$pn<-p
a<-merge(pdat.d,pdat.n,by=c('year'),all=FALSE)
a$dvm<-a$pd/a$pn
par(mfrow=c(2,2))
plot(a$year,a$dvm,ylim=c(0,100),pch=15,type='b')
plot(a$year,log10(a$dvm),pch=15,type='b',ylim=c(0,3))
plot(a$year,a$pd,pch=15)
plot(a$year,a$pn,pch=15)
plot(a$year,log10(a$pd+1),pch=15,type='b',ylim=c(-0,2))
points(a$year,log10(a$pn+1),pch=15,type='b',col='red')
