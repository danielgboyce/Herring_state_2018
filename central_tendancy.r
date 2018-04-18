library(fossil)
library(marmap)
library(segmented)
library(splines)
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
codedir<-'C:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/code/analysis'
codedir<-'C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/code/analysis'

setwd(codedir)
source('helper_functions.r')

setwd(datadir)
load('phytoplankton_biomass_all_spera_spawnar.RData');phyt<-dat
#load('phytoplankton_biomass_all_spera_allar.RData');phyt<-dat
path<-load('pathfinder_sst_daily_1981_2012_spawnar.RData');path<-dat
nit<-read.csv("Atlas1999nuts_nit_spawnar.csv",header=TRUE)
phos<-read.csv("Atlas1999nuts_phos_spawnar.csv",header=TRUE)
sil<-read.csv("Atlas1999nuts_sil_spawnar.csv",header=TRUE)
sz<-read.csv("phyto_size_monthly_1997_2007_spera_spawnar.csv",header=TRUE)
fg<-read.csv("phyto_functional_monthly_1997_2010_spera_spawnar.csv",header=TRUE)
strat<-read.csv("physical_stratification_spera_spawnar.csv",header=TRUE)
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvl<-read.csv("herring_lengths_RV_survey_spera_spawnar.csv",header=TRUE)
#zto<-read.csv("zoop_turnover_spera_spawnar.csv",header=TRUE)
#pto<-read.csv("phyto_turnover_spera_spawnar.csv",header=TRUE)
plank<-read.csv('plank_cpr_richness_spera_spawnar.csv',header=TRUE)
cur<-read.csv('GS_SS_distances.csv',header=TRUE)
cur$lon<-NA
cur$lat<-NA

wnd<-read.csv('petrie.wind.phen.csv',header=TRUE)
wnd$lon<-NA
wnd$lat<-NA
wnd$day<-(wnd$month*30)-15


#HERRING LARV
#larv<-read.csv("BoF_larval_herring_counts_spera_spawnar.csv",header=TRUE)
load("BoF_larval_herring_counts_spera_spawnar.RData")
larv<-subset(dat,bathy<0 & month>=7 & gearprocdesc=='Sawtooth Oblique')
#STANDARDIZE PER VOLUME STRAINED
larv$stdno<-(larv$totno/larv$volm3)*abs(larv$bathy)#NUMBER PER M2 - INTEGRATED OVER DEPTH
#larv$stdno<-(larv$totno/larv$volm3)
larv$stdno<-ifelse(is.na(larv$stdno)==TRUE,0,larv$stdno)
#SPAWNING ON GERMAN BANK AUG-SEPT, EARLIER ON LURCHER
larv$phylum<-tolower(larv$phylum)
larv$order<-tolower(larv$order)
larv$kingdom<-tolower(larv$kingdom)
larv$sid<-gsub(' ','',paste(larv$cruise,'_',larv$setno,'_',larv$sample))



#SUM OF COUNTS FOR LARVAL PREY SPECIES
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
#OUTPUT HERRING LARVAE DATA
herlrv.mn.bof<-ddply(subset(larv,comname %in% c('herring(atlantic)','herring/capelin like')),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
lrv.mn.bof<-ddply(subset(larv,!(comname %in% c('herring(atlantic)','herring/capelin like'))),.(lon,lat,year,month,day,cruise,setno,sample),.fun=avfun,.progress='text')
herlrv.mn.bof$var<-'herlrv.mn.bof'
lrv.mn.bof$var<-'lrv.mn.bof'
lrv<-rbind.fill(zp1,zp2,zp,herlrv.mn.bof,lrv.mn.bof)



#ZOOPLANKTON
zoop<-read.csv("zoop_cpr_counts_spera_spawnar.csv",header=TRUE)
zoop<-subset(zoop,scientificname %in% c('Calanoida','Oithona','Calanus', 'Calanus finmarchicus','Calanus glacialis','Calanus hyperboreus','Candacia','Candacia armata','Centropages','Centropages hamatus','Centropages typicus','Pseudocalanus elongatus','Temora longicornis','Copepoda'))
f<-function(d){
    return(data.frame(zoop=sum(d$counts)))
}
zoop<-ddply(zoop,.(lon,lat,year, month,djul,day,tbin3,tbin5,tbin10),.fun=f,.progress='text')
zoop$var<-'zp.cpr'




f<-function(d){
    names(d)[1:2]<-c('y','x')
    plot(d$x,d$y)
}
par(mfrow=c(4,5))
f(subset(wnd,select=c('wnd.pctabove10','day')))
f(subset(wnd,select=c('wnd.spd323','day')))
f(subset(wnd,select=c('wnd.strss323','day')))
f(subset(plank,select=c('richness.z','day')))
f(subset(plank,select=c('copepods.z','day')))
f(subset(plank,select=c('shannon.p','day')))
f(subset(path,select=c('wind','day')))
f(subset(path,select=c('sst','day')))
f(subset(strat,select=c('s25','day')))
f(subset(strat,select=c('temp50','day')))
f(subset(zoop,select=c('zoop','day')))
f(subset(larv,select=c('stdno','day')))
f(subset(phyt,select=c('chl','day')))
f(subset(nit,select=c('nit','day')))
f(subset(phos,select=c('phos','day')))
f(subset(lrv,var=='Mollusc',select=c('stdno','day')))



#FOR OCEANOGRAPHIC VARS, THAT ARE UNIMODAL BUT DON'T PEAK IN THE FALL, CENTER SO THAT THE PEAK OCCURS IN THE CENTER OF THE ANNULAL SERIES (DAY=182.5)

d<-subset(wnd,select=c('wnd.pctabove10','day','year','month','lon','lat'))
var<-'wnd'
tm<-'all'

centerfun<-function(d,var,tm){
    #DETERMINE DAY OF MAX
    names(d)[1]<-'y'
    mod<-gam(y~s(day,bs='cs',k=5),data=d,gamma=1.4)
    pdat<-data.frame(xx=seq(min(d$day),max(d$day),1))
    pdat$p<-predict(mod,newdata=data.frame(day=pdat$xx))
    mday<-subset(pdat,p==max(pdat$p))$xx[1]
    #CENTER OF TIMING OF MAX
    shft<-(365/2)-mday
    d$day2<-d$day+shft
    d$day2<-ifelse(d$day2>365,d$day2-365,d$day2)
    d$day2<-ifelse(d$day2<0,d$day2+365,d$day2)
    d$day2<-ceiling(d$day2)
    d$seas<-tm
    d$var<-var
    if(var=='chl'){
        dout<-subset(d,is.na(y)==FALSE,select=c('year','month','day','day2','y','var','seas','db','cl','lon','lat'))
    } else {
        dout<-subset(d,is.na(y)==FALSE,select=c('year','month','day','day2','y','var','seas','lon','lat'))
    }

    if(var=='nit'){
        dout$db<-'NIT'
        dout$var<-'nut'
    } else if(var=='phos'){
        dout$db<-'PHOS'
        dout$var<-'nut'
    } else if(var=='sil'){
        dout$db<-'SIL'
        dout$var<-'nut'
    } else NULL

    return(dout)
}
l<-list()
l[[1]]<-centerfun(subset(lrv,var=='herjuv.prey.bof' & day>250,select=c('stdno','day','year','month','tbin3','tbin5','tbin10','var','lon','lat')),'herjuv.prey.bof','fall')
l[[2]]<-centerfun(subset(lrv,var=='her.prey.bof' & day>250,select=c('stdno','day','year','month','tbin3','tbin5','tbin10','var','lon','lat')),'her.prey.bof','fall')
l[[3]]<-centerfun(subset(lrv,var=='her.preytot.bof' & day>250,select=c('stdno','day','year','month','tbin3','tbin5','tbin10','var','lon','lat')),'her.preytot.bof','fall')
l[[4]]<-centerfun(subset(lrv,var=='herlrv.mn.bof' & day>250,select=c('stdno','day','year','month','tbin3','tbin5','tbin10','var','lon','lat')),'herlrv.mn.bof','fall')
l[[5]]<-centerfun(subset(lrv,var=='lrv.mn.bof' & day>250,select=c('stdno','day','year','month','tbin3','tbin5','tbin10','var','lon','lat')),'lrv.mn.bof','fall')
l[[6]]<-centerfun(subset(zoop,day>75,select=c('zoop','day','year','month','tbin3','tbin5','tbin10','var','lon','lat')),'zp.cpr','fall')
l[[7]]<-centerfun(subset(strat,day>=175,select=c('temp50','day','year','month','tbin3','tbin5','tbin10','lon','lat')),'t50','all')
l[[8]]<-centerfun(subset(cur,select=c('GS.Dist','day','year','tbin3','tbin5','tbin10','month','lon','lat')),'gs','all')
l[[9]]<-centerfun(subset(cur,select=c('SS.Dist','day','year','tbin3','tbin5','tbin10','month','lon','lat')),'ss','all')
l[[10]]<-centerfun(subset(plank,day>210,select=c('totalabundance.p','day','year','month','tbin3','tbin5','tbin10','lon','lat')),'ph.mn.cpr','fall')
l[[11]]<-centerfun(subset(plank,day>210,select=c('totalabundance.z','day','year','month','tbin3','tbin5','tbin10','lon','lat')),'zp.mn.cpr','fall')
l[[12]]<-centerfun(subset(plank,day>210,select=c('copepods.z','day','year','month','tbin3','tbin5','tbin10','lon','lat')),'cp.cop.cpr','fall')
l[[13]]<-centerfun(subset(phyt,day>210,select=c('chl','day','year','month','tbin3','tbin5','tbin10','lon','lat','db','cl')),'chl','fall')
l[[14]]<-centerfun(subset(nit,select=c('nit','day','year','month','tbin3','tbin5','tbin10','lon','lat')),'nit','all')
l[[15]]<-centerfun(subset(phos,select=c('phos','day','year','month','tbin3','tbin5','tbin10','lon','lat')),'phos','all')
l[[16]]<-centerfun(subset(sil,select=c('sil','day','year','month','tbin3','tbin5','tbin10','lon','lat')),'sil','all')
l[[17]]<-centerfun(subset(strat,select=c('s25','day','year','month','tbin3','tbin5','tbin10','lon','lat')),'strt','all')
l[[18]]<-centerfun(subset(phyt,day<210 & day>50,select=c('chl','day','year','month','tbin3','tbin5','tbin10','lon','lat','db','cl')),'chl','spring')
l[[19]]<-centerfun(subset(wnd,select=c('wnd.pctabove10','day','year','month','lon','lat')),'wnd','all')
dat<-rbind.fill(l)

#MAKE UNIQUE DB FOR EACH VAR
dat$db<-ifelse(!(dat$var %in% c('nut','chl')),dat$var,dat$db)

#MERGE BATHYMETRY DATA FROM NOAA
crds<-na.omit(unique(subset(dat,select=c('lon','lat'))))
crds$id<-gsub(' ','',paste(crds$lon,'_',crds$lat))
atl<-getNOAA.bathy(lon1=min(crds$lon),lon2=max(crds$lon),lat1=min(crds$lat),lat2=max(crds$lat), resolution=.1)
bathy<-fortify.bathy(atl)
names(bathy)<-c('lon','lat','bathy')


bfun<-function(d){
bathy2<-bathy
bathy2$dist<-deg.dist(long1=d$lon,lat1=d$lat,long2=bathy2$lon,lat2=bathy2$lat)
dd<-subset(bathy2,dist==min(bathy2$dist))[1,]
d$bathy<-dd$bathy
d$dist<-dd$dist
d$n<-dim(dd)[1]
d$nunique<-length(unique(dd$bathy))
return(d)
}
crdso2<-dlply(crds,.(id),.fun=bfun,.progress='text')
crdso2<-data.frame(do.call('rbind',crdso2))
#crdso<-dlply(crds,.(id),.fun=bfun,.progress='text')
#crdso<-data.frame(do.call('rbind',crdso))
#save(crdso,file='bathy_ct.RData')
load('N:/cluster_2017/scratch/spera/data/finaldat_v2/bathy_ct.RData')

a<-subset(crdso,select=c('lon','lat','bathy','dist'))
dat<-merge(dat,a,by=c('lon','lat'),all.x=TRUE,all.y=FALSE)

dat$y<-ifelse(dat$db %in% c('SW','CZ','MO','ME','INSITU'),log10(dat$y+1),dat$y)


#BECAUSE SOME YEARS HAVE FEW OBSERVATIONS DURING PEAK, NEED TO RESTRICT DATA TO THE PEAK WINDOW WE'RE INTERESTED IN OR ELSE PHENOLOGY CALCS WON'T WORK
f<-function(d){
if(unique(d$var) %in% c('nut')){
    d<-subset(d,day2>=100 & day2<=300)
} else if (unique(d$var) %in% c('chl') & unique(d$seas) %in% c('fall')){
    d<-subset(d,day2<=250)
} else if (unique(d$var) %in% c('chl') & unique(d$seas) %in% c('spring')){
    d<-subset(d,day2>=175)
} else if (unique(d$var) %in% c('strt')){
    d<-subset(d,day2>=100 & day2<=250)
} else if (unique(d$var) %in% c('t50')){
    d<-subset(d, day2<=225)
} else NULL
return(d)
}
dat<-ddply(dat,.(var,seas),.fun=f,.progress='text')

f<-function(d){
    plot(d$day2,d$y)
}
par(mfrow=c(4,5))
f(subset(dat,var=='wnd'))
f(subset(dat,var=='ph.mn.cpr'))
f(subset(dat,db=='INSITU' & seas=='fall'))
f(subset(dat,db=='INSITU' & seas=='spring'))

##########################################################3
#CENTRAL TENDANCY FOR NON PHYTOPLANKTON
#FIRST DETERMINE APPROPRIATE DATA AVAILABILITY
dat$id<-gsub(' ','',paste(dat$db,'_',dat$seas,'_',dat$year))
f<-function(d){
#    d<-subset(d,is.na(bathy)==FALSE)
    return(data.frame(n=dim(d)[1],
                      nmo=length(unique(d$month)),
                      nlon=length(unique(d$lon)),
                      nlat=length(unique(d$lat)),
                      ny=length(unique(d$y))))
}
dm<-ddply(dat,.(id),.fun=f,.progress='text')
dm<-subset(dm,nmo>=2 & n>=10 & ny>1)



d<-subset(dat,id=="CPR_spring_1979")

ctfun<-function(d){
#d<-subset(d,is.na(bathy)==FALSE)
d0<-subset(d,is.na(bathy)==FALSE)
print(unique(d$id))
    names(d)[3]<-'time'
    nmo<-ifelse(unique(d$seas)=='fall',2,6)
    d$y<-d$y+(min(d$y)+.1)

#CTFUNCTION
fn<-function(dd){
#CT USING AREA UNDER CURVE
dum<-data.frame(day=sort(unique(dd$day2)),
                sump=tapply(dd$y,dd$day2,mean))
dum$ct1<-dum$day*dum$sump
ct<-sum(dum$ct1)/sum(dum$sump)
#CT USING WEIGHTED MODEL APPRAOCH
mod<-gam(day2~bathy,weights=dd$y,data=dd)
s<-summary(mod)
mod2<-gam(day2~s(lon,lat,k=6),weights=dd$y,data=dd)
s2<-summary(mod2)
return(data.frame(ct=ct,
                  ctm1=s$p.table[1,1],
                  ctm2=s2$p.table[1,1],
                  ctw=weighted.mean(dd$day,w=dd$y)))
}

fn2<-function(dd){
#CT USING AREA UNDER CURVE
dum<-data.frame(day=sort(unique(dd$day2)),
                sump=tapply(dd$y,dd$day2,mean))
dum$ct1<-dum$day*dum$sump
ct<-sum(dum$ct1)/sum(dum$sump)
return(data.frame(ct=ct,
                  ctm1=NA,
                  ctm2=NA,
                  ctw=weighted.mean(dd$day,w=dd$y)))
}

if((unique(d$db) %in% c('gs','ss','wnd') | length(unique(d$lon))<6) | dim(d0)[1]==0){ dout<-fn2(d)
} else { dout<-fn(d) }

dout$id<-unique(d$id)
dout$time<-unique(d$time)
dout$nmonth<-length(unique(d$month))
dout$n<-dim(d)[1]
dout$var<-unique(d$var)
dout$db<-unique(d$db)
dout$seas<-unique(d$seas)
return(dout)
}

pdat<-dlply(subset(dat,id %in% dm$id),.(id),.fun=ctfun)
pdat<-rbind.fill(pdat)
pdat$id<-gsub(' ','',paste(pdat$db,'_',pdat$seas))

f<-function(d){
    d$ctav<-mean(c(d$ct,d$ctm2),na.rm=TRUE)
    d$ctav.sd<-sd(c(d$ct,d$ctm2),na.rm=TRUE)
return(d)
}
pdat<-ddply(pdat,.(db,time),.fun=f,.progress='text')


xyplot(ctav~time | id,data=pdat,pch=15)
xyplot(ctm1~time | id,data=pdat,pch=15)
xyplot(ctm2~time | id,data=pdat,pch=15)
xyplot(ctw~time | id,data=pdat,pch=15)

#REMOVE INSTANCES WHERE TREND IS COMPLETELY FLAT- NO PURPOSE
pdat<-subset(pdat,!(var %in% c('herjuv.prey.bof','herlrv.mn.bof','lrv.mn.bof','gs','her.prey.bof','her.preytot.bof','ph.mn.cpr')))

xyplot(ct~time | var,data=pdat,pch=15,xlim=c(1960,2017))
xyplot(ct~time | var,data=subset(pdat,db=='all'),pch=15,xlim=c(1960,2017))


f<-function(d){
out<-data.frame(year=sort(unique(d$time)),
                ct=tapply(d$ct,d$time,mean))
names(out)[2]<-gsub(' ','',paste(unique(d$var),'.ct'))
return(out)
}
pdat2<-dlply(pdat,.(var),.fun=f)
ctdat<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), pdat2)#COMBINE


setwd(datadir)
save(ctdat,file='ct_data.RData')








#########################################################################

#########################################################################

####PHENOLOGY
#dat<-subset(data,year>=2000)

f1<-function(d0){

f<-function(d){
day<-unique(d$day)
d2<-subset(d0,select=c('day'))#GET ALL DAYS IN ID2 OBJECT
d2$day2<-d2$day-(day)
d2$day2<-ifelse(d2$day2<=0,d2$day2+(365),d2$day2)
    return(data.frame(day=day,
                   span=diff(range(min(d2$day2),max(d2$day2)))))
}
z<-ddply(d0,.(db,day),.fun=f)
z<-subset(z,span==max(z$span))
return(subset(z,day==min(z$day)))
}
dayshiftcell<-dlply(unique(subset(dat,select=c('db','day'))),.(db),.fun=f1,.progress='text',.parallel=FALSE)
dayshiftcell<-rbindlist(dayshiftcell)


f2<-function(d){
d2<-subset(dayshiftcell,db==unique(d$db))
day<-d2$day#NUMER OF DAYS TO SHIFT
span<-d2$span
d$day2<-d$day-(day)
d$day2<-ifelse(d$day2<=0,d$day2+365,d$day2)
d$celloffset<-day
d$span<-span
return(d)
}
dat<-dlply(dat,.(db),.fun=f2,.progress='text',.parallel=FALSE)
dat<-rbindlist(dat)



##############################################
#USE CALCULATIONS ABOVE TO ADJUST DAY VALUES FOR EACH CELL AND DB TO MAXIMIZE PHENOLOGY SPAN

d<-subset(dat,db=='CPR')

phenfun<-function(d){
#print(unique(d$cell))
options(warn=-1)

if(length(unique(d$year))>1 & dim(d)[1]>=10 & length(unique(d$month))>=8){
#FIT MODEL - THE WORKHORSE
mod<-gam(log10(chl+.1)~s(day,by=factor(year),bs='cc',k=6,id=1) + factor(year) + s(lon,lat,k=4),data=d,select=TRUE)
#mod<-try(gam(log10(chl+.1)~s(day2,bs='cc',k=6,id=1) + as.factor(year),data=d,select=TRUE,gamma=1.4))
s<-summary(mod)

#FORMATS DATA TO PREDICT AT
    datout<-expand.grid(day=seq(1,365,1),
                        year=sort(unique(d$year),decreasing=TRUE))

    datout$lon<-median(d$lon)
    datout$lat<-median(d$lat)


#PREDICT
pred<-predict(mod,newdata=datout,na.action=na.pass,se.fit=TRUE)
datout$pl<-pred$fit
datout$p<-10^datout$p-.1
datout$sel<-pred$se.fit
datout$se<-10^pred$se.fit-.1


#Z-STANDARDIZE FOR EACH DB
f<-function(d){
d$pz<-round((d$p-mean(d$p))/sd(d$p),digits=4)
d$plz<-round((d$pl-mean(d$pl))/sd(d$pl),digits=4)
d$sez<-round((d$se-mean(d$se))/sd(d$se),digits=4)
d$selz<-round((d$sel-mean(d$sel))/sd(d$sel),digits=4)
d$p<-round(d$p,digits=4)
d$pl<-round(d$pl,digits=4)
return(d)
}
datout<-f(datout)

cls<-rainbow(length(unique(datout$year)))
plot(0,0,xlim=c(0,365),ylim=c(-2,3),main=unique(d$db),col='white')
    yearss<-sort(unique(datout$year))
    for(i in 1:length(yearss)){
    dd<-subset(datout,year==yearss[i])
    lines(dd$day,dd$pz,col=cls[i],lwd=2)
    }
    legend('topleft',legend=yearss,col=cls,pch=15,bty='n')

    if(length(unique(datout$p))<=2){nunique<-1
    } else nunique<-2

datout$id<-unique(d$id)
datout$db<-unique(d$db)
datout$lat<-median(d$lat)
datout$lon<-median(d$lon)
datout$aic<-AIC(mod)
datout$dev<-s$dev.expl
datout$r2<-s$r.sq
datout$celloffset<-unique(d$celloffset)
datout$nunique<-nunique


#GETS DATA AVAILABILITY FOR ID2
datout$n.id<-dim(d)[1]
datout$nmonths.id<-length(unique(d$month))
datout$span.id<-max(d$day2,na.rm=TRUE)-min(d$day2,na.rm=TRUE)

return(datout)

}else NULL
}

setwd(figsdir)
pdf('spera_pheno_trends.pdf')
z<-dlply(dat,.(db),.fun=phenfun,.parallel=FALSE,.progress='text')
z<-rbindlist(z)
z<-subset(z,nunique==2)#OMIT INSTANCES WHERE PHENOLOGY IS FLAT
dev.off()

#ADD MONTHS
date<-strptime(gsub(' ','',paste(2000,'-',ceiling(z$day))),"%Y-%j")
z$month<-month(date)


d<-subset(z,id==sort(unique(z$id))[50] & year==2012)
d<-subset(z,id==sort(unique(z$id))[50])
d<-subset(d,month>4)
dum<-data.frame(month=sort(unique(d$month)),
                sump=tapply(d$p,d$month,mean))

dum$ct1<-dum$month*dum$sump
sum(dum$ct1)/sum(dum$sump)


d<-subset(z,id==sort(unique(z$id))[50])
d<-subset(d,month>4)


f<-function(d){
    dd<-subset(d,month>=5)
    return(data.frame(nmonth=length(unique(d$month)),
                      n=dim(d)[1],
                      wind=length(unique(dd$month))))
}
dum<-ddply(data,.(id),.fun=f,.progress='text')
dum<-subset(dum, wind>=7)
dat2<-subset(dat,id %in% dum$id & day>=184)



ctfun<-function(d){
dum<-data.frame(day=sort(unique(d$day)),
                sump=tapply(d$chl,d$day,mean))

dum$ct1<-dum$day*dum$sump
ct<-sum(dum$ct1)/sum(dum$sump)
return(data.frame(year=unique(d$year),
                  n=dim(d)[1],
                  db=unique(d$db),
                  ct=ct,
                  cl=unique(d$cl)))
}
ctdat<-ddply(dat2,.(id),.fun=ctfun,.progress='text')


plot(ctdat$year,ctdat$ct,pch=15,col=as.character(ctdat$cl),las=1,ylim=c(184,335))
abline(h=mean(ctdat$ct,na.rm=TRUE),lty=2)
hist(ctdat$ct,breaks=50)
plot(ctdat$year,ctdat$ct,pch=16,col=as.character(ctdat$cl),cex=rescale(ctdat$n,newrange=c(.75,3)))
mod<-lm(ct~year,data=ctdat)
abline(mod)
abline(mod2)
mod2<-lm(ct~year,data=subset(ctdat,db!='CPR'))

plot(d$day,d$p)
abline(v=242)


lnfunn<-function(ddd){
y<-unique(subset(ddd,select=c('day2')))
y$dm<-rep(1,dim(y)[1])
days<-data.frame(day2=seq(1,366,1))
y<-merge(days,y,by=c('day2'),all.x=TRUE,all.y=FALSE)
y$dm<-ifelse(is.na(y$dm)==TRUE,-9,y$dm)
y <- rle(y$dm)
options(warn=-1)
y2<-max(y$lengths[y$values==-9])#GETS LONGEST STRETCH OF MISSING DATA
options(warn=0)
if(y2 %in% c(-Inf,Inf)){return(0)
}else {return(y2)}
}

dum<-ddply(dat,.(db,year),.fun=lnfunn,.progress='text')


a<-subset(dat,db=='SW' & year==2002)
plot(a$chl,a$day)










a<-subset(dat,db=='ME' & year==2009)
f<-function(aa){
map('world',xlim=c(-70,-63),ylim=c(42,45))
points(a$lon,a$lat,pch=16,col='red')
}
par(mfrow=c(3,4))
dlply(a,.(month),.fun=f)
plot(a$day,a$chl,pch=15)


#PUTS TIMESERIES IN PROPER FORMAT TO CALCULATE CORRELATIONS
fun3<-function(d){
    nar<-unique(d$cell)
    dt<-subset(d,select=c('day','pz'))
    dt<-dt[order(dt$day),]
    names(dt)[2]<-as.character(nar)
    return(dt)
}
dat2<-dlply(z,.(cell),.fun=fun3,.progress='text')
q<-Reduce(function(x, y) merge(x, y, by=c('day'),all=TRUE), dat2)#COMBINE
q<-q[,colSums(is.na(q)) != nrow(q)]#REMOVES COLUMNS THAT ARE ALL MISSING
rownames(q)<-q$day
q<-q[,-1]









df<-q
k<-4#NUMBER OF CLUSTERS
dmat<-1-cor(df,use='pairwise.complete.obs')
dst<-as.dist(dmat)
ff<-fanny(dst,k,maxit=5000,diss=T)


setwd(figsdir)
pdf('chl_phen_cluser_spera.pdf',width=7,height=5)
#par(mfrow=c(2,1),mar=c(4,4,1,1),oma=c(5,1,5,1))
par(mar=c(4,4,4,4),oma=c(4,1,1,1))
dum<-c('red3','magenta3','forestgreen','gold','cornflowerblue','darkblue')
#dum<-c('red3','cornflowerblue','forestgreen','gold','magenta3','darkblue')
#dum<-c('red3','green3','purple','gold3','orangered2','blue3')
plot(silhouette(ff),col=dum[1:k],main='')#silhouette plot

dc.pcoa<-cmdscale(dst)
dc.scores<-scores(dc.pcoa,choices=c(1,2))
spefuz.g<-ff$clustering
a<-data.frame(cell=as.character(sort(unique(names(df)))),
              clusters=ff$clustering)
aa<-data.frame(ff$membership)
aa$cell<-rownames(a)

par(mar=c(1,1,1,8),oma=c(1,1,1,1))
plot(scores(dc.pcoa),asp=1,type='n',xlim=c(-1.1,1),ylim=c(-1.2,1),las=1,axes=FALSE,xlab='',ylab='')
stars(ff$membership,location=scores(dc.pcoa),draw.segments=T,add=T,scale=F,len=.1,col.segments=alpha(c(dum[1:k]),.25),byt='n',labels=NULL,xlim=c(-1.1,1.1),ylim=c(-1.2,1),lwd=.0001,xpd=TRUE,border=NULL,radius=FALSE,col.radius=alpha('white',.1))
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
par(mfrow=c(1,4),mar=c(2,.5,1,.5),oma=c(14,2,14,.5))

mb<-seq(1,k,1)
l<-list()
for(i in 1:length(mb)){
print(mb[i])
one<-subset(a,clusters==mb[i],select=c('cell'))#CELLS IN THIS CLUSTER
cx2<-subset(cx,cx>=0.5)#GETS INSTANCES WHERE CLUSTER PROBABILITY>0.5
data5<-subset(df,select=c(as.character(one$cell)))#TS FOR CLUSTER
cx2<-subset(cx2,cell %in% names(data5))
data5<-subset(df,select=c(as.character(cx2$cell)))#TS FOR CLUSTER
x<-rownames(data5)
t<-data.frame(day=as.numeric(x),mn=rowMeans(data5,na.rm=F))
t2<-data.frame(day=as.numeric(x),mn=rowMeans(data5,na.rm=T))

cl<-dum[i]
plot(0,0,pch=16,cex=.01,xlim=c(0,365),ylim=c(-2,3),main='',xaxt='n',las=1,axes=FALSE,xlab='',ylab='')
abline(v=182.5,lwd=.5,col='gray50')
if(i==1){axis(side=2,at=seq(-2,3,1),las=1,lwd=.001,cex.axis=.75)
}else NULL
axis(side=1,at=c(0,182.5,365),labels=c(-182.5,0,182.5),cex.axis=.75)
    for(j in 1:length(data5[1,])){
        try(dat<-data5[,j])
        try(trnsp<-subset(cx,cell==as.character(names(data5[j])))$cxx)
        try(lines(x,dat,cex=.5,ylim=c(0,1),lwd=1.25,col=alpha(cl,rescale(trnsp,newrange=c(.01,.25)))))
    }
points(as.numeric(as.character(t$day)),t$mn,pch=16,col='gold3',cex=.7)
dm<-subset(t2,mn==max(t2$mn,na.rm=TRUE),select=c('day'))
dm$cluster=mb[i]
names(dm)<-c('clusterday','clusters')
l[[i]]<-dm
}


ll<-data.frame(do.call('rbind',l))
dev.off()


#######################################
#PLOTS SPATIAL DISTRIBUTION OF CLUSTERS
#CALCULATES AVERAGE LOCATION OF EACH STRATUM
cx2<-merge(cx,a,by=c('cell'),all=FALSE)

cls<-data.frame(clusters=sort(unique(cx2$clusters)),
                cl=dum[1:length(unique(cx2$clusters))])

cx2<-merge(cx2,cls,by=c('clusters'),all.x=TRUE,all.y=FALSE)

#ADDS COORDINATES
crds<-unique(subset(z,select=c('cell','lon','lat')))
cx2$cell<-as.factor(cx2$cell)
crds$cell<-as.factor(crds$cell)
cx2<-merge(cx2,crds,by=c('cell'),all.x=TRUE,all.y=FALSE)


##GGPLOT MAP OF CLUSTERS
ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }


dum<-unique(subset(cx2,select=c('clusters','cl')))
cx2$cxr<-rescale(cx2$cx,newrange=c(.05,.95))


xlm1<--68
xlm2<--64
ylm1<-42
ylm2<-46

setwd(figsdir)
pdf('chl_phen_cluser_map_spera.pdf',width=5,height=5)
ggplot()+
geom_tile(data=cx2, aes(x=lon, y=lat,fill=as.factor(clusters),alpha=cxr),color='white')+
scale_fill_manual(breaks=dum$clusters,values=as.character(dum$cl))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",plot.background=element_blank())+
coord_equal()+
scale_x_continuous(expand=c(0,0),breaks=seq(xlm1,xlm2,1),labels=seq(xlm1,xlm2,1),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(ylm1,ylm2,1),labels=seq(ylm1,ylm2,1),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(ylm1,ylm2),xlim=c(xlm1,xlm2))+
    xlab('')+
    ylab('')

dev.off()













############################################################

###DECORRELATION

#####################################################################################
#CALCULATES CORRELATION BETWEEN ALL CHL TIMESERIES
cr<-cor(q,use='pairwise.complete.obs')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
names(cr.t)<-c('cell.x','cell.y','r')
#ADDS GEOGRAPHIC COORDINATES
crds<-unique(subset(dat,select=c('cell','lon','lat')))
names(crds)[1]<-'cell.x'
cr.t<-merge(cr.t,crds,by=c('cell.x'),all.x=TRUE,all.y=FALSE)
names(crds)[1]<-'cell.y'
cr.t<-merge(cr.t,crds,by=c('cell.y'),all.x=TRUE,all.y=FALSE)
#CALCULATES DISTANCES
cr.t$id<-seq(1,dim(cr.t)[1],1)
f<-function(d){
    d$dist<-deg.dist(d$lon.x,d$lat.x,d$lon.y,d$lat.y)
    d$cell.x<-as.character(d$cell.x)
    d$cell.y<-as.character(d$cell.y)
    return(d)
}
cr<-dlply(cr.t,.(id),.fun=f,.progress='text')
gc()
cr<-rbindlist(cr,use.names=TRUE)








#####################################################################################
#FUNCTION TAKES CORRELATIONS BETWEEN ALL TS AND ESTIMATES DISTANCE DECORRELATION SCALE
dt<-unique(subset(dat,select=c('cell','lon','lat')))
dt$id<-seq(1,dim(dt)[1],1)
#dt<-subset(dt,!(cell %in% c("-65.854_45.063", "-65.896_43.688","-65.688_45.063")))

#m<-subset(datcr,decor==-Inf)
d<-subset(dt,cell=="-67.396_44.521")


f<-function(d){
#GET ALL INSTANCES INVOLVING SPECIFIC CELL FROM CORRELATION/DISTANCE DB
    print(unique(d$cell))
    dst<-500
    a<-subset(cr,dist<=dst & (cell.x %in% unique(d$cell) | cell.y %in% unique(d$cell)))
    a<-subset(a,is.na(r)==FALSE)
if(length(unique(na.omit(a$r)))>=5 & max(a$dist)>=100){
a<-rbind.fill(a,data.frame(dist=rep(0,20),r=rep(1,20)))
brks<-seq(min(a$dist)-.001,max(a$dist)+.001,length.out=100)
df<-diff(brks)[1]
lbls<-c(brks+df)[1:length(brks)-1]
a$dbin<-round(as.numeric(as.character(cut(a$dist,breaks=brks,labels=lbls),digits=0)))

dm<-data.frame(dbin=sort(unique(a$dbin)),
               r=tapply(a$r,a$dbin,function(x) mean(x,na.rm=TRUE)))
plot(dm$dbin,dm$r,pch=16,ylim=c(-.5,1),xlim=c(0,500),las=1,xlab='Distance',ylab='Correlation',axes=F,col=alpha('black',.5))
axis(1,at=seq(100,dst,length.out=5))
axis(2,at=seq(-.5,1,.2),las=1,labels=FALSE,tick=TRUE)
abline(h=0,lty=3,lwd=.5)
mod<-gam(r ~ s(dbin,k=6),data=a,gamma=1.4)
s<-summary(mod)
x<-seq(0,max(a$dbin),length.out=10000)
p<-predict(mod,newdata=data.frame(dbin=x),se.fit=TRUE)
lines(x,p$fit)
mtext(as.character(round(as.numeric(d$lat),digits=1)),side=3,line=-1,cex=.75)
dm<-data.frame(x=x,
               p=round(p$fit,digits=2),
               se=round(p$se.fit,digits=4))
dc<-subset(dm,p==.37)

if(dim(dc)[1]>0){
    dcall<-subset(dc,x==min(dc$x))
points(dcall$x,0,pch=18,col='black',cex=3)
} else { dcall<-data.frame(x=-Inf,
                           p=NA,
                           se=NA)
}
d$decor<-round(dcall$x,digits=0)
d$decor.se<-round(dcall$se,digits=4)
d$decor.int<-s$p.table[1,1]

#LINEAR SLOPE NEAR ORIGIN - TELLS IF POSSIBLE TO DERIVE DECORRELATION SCALE
modl<-lm(r ~ dbin,data=subset(a,dist<150))
d$lslope<-modl$coef[2]

return(d)
} else NULL
}

setwd(figsdir)
pdf('decorr_chl_phen_spera.pdf')
par(mfrow=c(5,4),mar=c(2,1,1,1))
phendatcr<-ddply(dt,.(id),.fun=f,.progress='text')
dev.off()












library(rgdal)
par(lwd=.01)
ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }


data('wrld_simpl',package='maptools')
plg<-crop(wrld_simpl,extent(-180,180,-80,80),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
world.df<-fortify(plg)

brmn<-'+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'
mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
brmn<-'+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'
mc<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast<-readOGR('/misc/scratch/dboyce/shapefiles/naturalearthdata_ne_10m_land_poly',layer='ne_10m_land')#works for Alaska - most others don't
#coast<-readOGR('/misc/scratch/dboyce/chl_phenology/data/naturalearthdata_ne_110m_land_poly',layer='ne_110m_land')#works for Alaska - most others don't
fcoast<-fortify(coast)
coast.mc<-crop(coast,extent(-180,180,-90,90),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fcoast.mc<-fortify(coast.mc)
coast.mw<-spTransform(coast.mc, CRS=(brmn))
fcoast.mw<-fortify(coast.mw)
#coast.mw<-spTransform(coast.mc, CRS=(brmn))



adat<-subset(phendatcr,select=c('decor','lon','lat','cell'))
ttl<-'ED'
mn<-0
mx<-350
dg<-1

#FUNCTION TO RETURN MAP
nm<-names(adat)[1]
names(adat)[1]<-'y'
adat$y<-ifelse(adat$y>mx,mx,adat$y)
adat$y<-ifelse(adat$y<mn,mn,adat$y)
adat$y<-round(adat$y,digits=2)
a<-adat

n<-21
brks<-seq((min(adat$y,na.rm=TRUE)-.01),max(adat$y,na.rm=TRUE)+.01,length.out=n)
brks2<-round(seq((min(adat$y,na.rm=TRUE)-.01),max(adat$y,na.rm=TRUE)+.01,length.out=n),digits=dg)
a$ycat<-cut(a$y,breaks=brks)
#b$ycat<-cut(b$y,breaks=brks)
lbls<-sort(unique(a$ycat))
lbls2<-sort(unique(cut(a$y,breaks=brks2)))
cls<-matlab.like(length(lbls))
#cls<-colorRampPalette(c('magenta4','blue3','green','yellow','red3'))
#cls<-(cls(length(lbls)))

return(
ggplot()+
#geom_raster(data=a, aes(x=lon, y=lat,fill=ycat),interpolate=TRUE)+
geom_tile(data=a, aes(x=lon, y=lat,fill=ycat))+
scale_fill_manual(breaks=as.character(lbls),values=cls,labels=lbls2,na.value="transparent",guide=guide_legend(title=paste(ttl)))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)
+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6),plot.background = element_rect(fill = 'green', colour = 'red'))+
scale_x_continuous(expand=c(0,0),breaks=c(-67,-63),labels=as.character(c(-67,-63)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=c(42,46),labels=as.character(c(42,46)),limits=NA)+
#scale_x_continuous(expand=c(0,0),breaks=sort(unique(mlims$lone)),labels=as.character(sort(unique(mlims$lond))),limits=NA)+
#scale_y_continuous(expand=c(0,0),breaks=sort(mlims$late),labels=as.character(mlims$latd),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(min(adat$lat),max(adat$lat)),xlim=c(min(adat$lon),max(adat$lon)))+
    xlab('')+
    ylab('')
)


library(ggmap)
library(rgdal)
nafo<-readOGR('/scratch/dboyce/shapefiles/nafo','Divisions')
nafo<-spTransform(nafo,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fourx<-subset(nafo,ZONE=='4X')
fourx2<-crop(fourx,extent(-180,180,40.9,47))

aut<-data.frame(lon=c(-65.25,-66,-66.25,-66.3),
                lat=c(43.5,43.2,43.9,43.2))
shl<-data.frame(loc=c('Trinity Ledge','German Bank','Lurcher Shoal','Scotts Bay'),
                lon=c(-66.3,-66.4,-66.47,-64.9),
                lat=c(44,43.3,43.87,45.2))
shl<-subset(shl,loc!='Lurcher Shoal')

pngMAP_df<- get_map(location = c(lon = -66.2, lat = 43.5), source = "google", zoom = 7,color='color',maptype='satellite',crop=TRUE)
p<-ggmap(pngMAP_df)

setwd(figsdir)
pdf('spera_spawning_map.pdf',height=7,width=7)
p+
geom_polygon(aes(x=long,y=lat,col='red'),data=fourx2,col=NA,fill='white',alpha=.3,lwd=.2)+
geom_point(aes(x=lon,y=lat,col='red'),data=shl,col='gold3',fill='gold3',alpha=1,size=5,shape=15)+
geom_label(aes(x=lon,y=lat,label=loc),data=shl,nudge_x=c(-.7,-.7,-.6))
dev.off()

setwd(figsdir)
pdf('chl_phen_decor_map_spera.pdf',height=10,width=8)
p+
geom_tile(data=a, aes(x=lon, y=lat,fill=ycat),color='gray90',size=.00001)+
    geom_point(aes(x=-66.25,y=43.5),color='red',size=2)+
    scale_fill_manual(breaks=as.character(lbls),values=cls,labels=lbls2,na.value="transparent",guide=guide_legend(title=paste(ttl)))+
scale_alpha(guide = 'none')+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6),plot.background = element_rect(fill = 'green', colour = 'red'))
dev.off()

deg.dist(-68,43,-69,43)
deg.dist(-68,43,-68.25,43)

+
scale_x_continuous(expand=c(0,0),breaks=c(-67,-63),labels=as.character(c(-67,-63)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=c(42,46),labels=as.character(c(42,46)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(min(adat$lat),max(adat$lat)),xlim=c(min(adat$lon),max(adat$lon)))+
    xlab('')+
    ylab('')













#############################


###########################


map('world',xlim=c(-68,-63),ylim=c(42,46))
pts<-unique(subset(dat,select=c('lon','lat')))
points(pts$lon,pts$lat,pch=16,cex=.5,col='red')

pts<-unique(subset(dd,select=c('clond','clatd')))
points(pts$clond,pts$clatd,pch=16,cex=.5,col='green')


#PUTS TIMESERIES IN PROPER FORMAT TO CALCULATE CORRELATIONS
fun3<-function(d){
    nar<-unique(d$cell)
    dt<-subset(d,select=c('tday','chl'))
    dt<-dt[order(dt$tday),]
    names(dt)[2]<-as.character(nar)
    return(dt)
}
dat2<-dlply(dat,.(cell),.fun=fun3,.progress='text')
q<-Reduce(function(x, y) merge(x, y, by=c('tday'),all=TRUE), dat2)#COMBINE
q<-q[,colSums(is.na(q)) != nrow(q)]#REMOVES COLUMNS THAT ARE ALL MISSING
rownames(q)<-q$tday
q<-q[,-1]



#####################################################################################
#CALCULATES CORRELATION BETWEEN ALL CHL TIMESERIES
cr<-cor(q,use='pairwise.complete.obs')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
names(cr.t)<-c('cell.x','cell.y','r')
#ADDS GEOGRAPHIC COORDINATES
crds<-unique(subset(dat,select=c('cell','lon','lat')))
names(crds)[1]<-'cell.x'
cr.t<-merge(cr.t,crds,by=c('cell.x'),all.x=TRUE,all.y=FALSE)
names(crds)[1]<-'cell.y'
cr.t<-merge(cr.t,crds,by=c('cell.y'),all.x=TRUE,all.y=FALSE)
#CALCULATES DISTANCES
cr.t$id<-seq(1,dim(cr.t)[1],1)
f<-function(d){
    d$dist<-deg.dist(d$lon.x,d$lat.x,d$lon.y,d$lat.y)
    d$cell.x<-as.character(d$cell.x)
    d$cell.y<-as.character(d$cell.y)
    return(d)
}
cr<-dlply(cr.t,.(id),.fun=f,.progress='text')
gc()
cr<-rbindlist(cr,use.names=TRUE)



setwd(datadir)
write.csv(cr,'cormatrix_chl_timeseries_meris.csv',row.names=FALSE)






#####################################################################################
#FUNCTION TAKES CORRELATIONS BETWEEN ALL TS AND ESTIMATES DISTANCE DECORRELATION SCALE

dt<-unique(subset(dat,select=c('cell','lon','lat')))
dt$id<-seq(1,dim(dt)[1],1)
dt<-subset(dt,!(cell %in% c("-66.188_44.438", "-65.896_43.688")))

#m<-subset(datcr,decor==-Inf)

d<-subset(dt,cell=="-66.188_44.438")


f<-function(d){
#GET ALL INSTANCES INVOLVING SPECIFIC CELL FROM CORRELATION/DISTANCE DB
    print(unique(d$cell))
    dst<-500
    a<-subset(cr,dist<=dst & (cell.x %in% unique(d$cell) | cell.y %in% unique(d$cell)))
    a<-subset(a,is.na(r)==FALSE)
if(length(unique(na.omit(a$r)))>=5 & max(a$dist)>=100){
a<-rbind.fill(a,data.frame(dist=rep(0,20),r=rep(1,20)))
brks<-seq(min(a$dist)-.001,max(a$dist)+.001,length.out=100)
df<-diff(brks)[1]
lbls<-c(brks+df)[1:length(brks)-1]
a$dbin<-round(as.numeric(as.character(cut(a$dist,breaks=brks,labels=lbls),digits=0)))

dm<-data.frame(dbin=sort(unique(a$dbin)),
               r=tapply(a$r,a$dbin,function(x) mean(x,na.rm=TRUE)))
plot(dm$dbin,dm$r,pch=16,ylim=c(-.5,1),xlim=c(0,500),las=1,xlab='Distance',ylab='Correlation',axes=F,col=alpha('black',.5))
axis(1,at=seq(100,dst,length.out=5))
axis(2,at=seq(-.5,1,.2),las=1,labels=FALSE,tick=TRUE)
abline(h=0,lty=3,lwd=.5)
mod<-gam(r ~ s(dbin,k=6),data=a,gamma=.8)
s<-summary(mod)
x<-seq(0,max(a$dbin),length.out=10000)
p<-predict(mod,newdata=data.frame(dbin=x),se.fit=TRUE)
lines(x,p$fit)
mtext(as.character(round(as.numeric(d$lat),digits=1)),side=3,line=-1,cex=.75)
dm<-data.frame(x=x,
               p=round(p$fit,digits=2),
               se=round(p$se.fit,digits=4))
dc<-subset(dm,p==.37)

if(dim(dc)[1]>0){
    dcall<-subset(dc,x==min(dc$x))
points(dcall$x,0,pch=18,col='black',cex=3)
} else { dcall<-data.frame(x=-Inf,
                           p=NA,
                           se=NA)
}
d$decor<-round(dcall$x,digits=0)
d$decor.se<-round(dcall$se,digits=4)
d$decor.int<-s$p.table[1,1]

#LINEAR SLOPE NEAR ORIGIN - TELLS IF POSSIBLE TO DERIVE DECORRELATION SCALE
modl<-lm(r ~ dbin,data=subset(a,dist<2000))
d$lslope<-modl$coef[2]

#EXPONENTIAL DECAY MODEL
#GET A PRIORI ESTIMATE OF PARAMETER
dm<-data.frame(dbin=sort(unique(a$dbin)),
               r=tapply(a$r,a$dbin,mean))
dm$rl<-log(dm$r+1.01)
md<-lm(rl~dbin,data=dm)
mds<-summary(md)
k.strt<-mds$coef[2,1]
lwr<-c(k.strt-0.01)
upr<-c(k.strt+0.001)
mod <- nls(r ~ 1*exp(-k * dbin), a, start=c(k=k.strt),control=nls.control(maxiter=500),lower=lwr,upper=upr,algorithm='port')
k<-summary(mod)$coef[1,1]
x<-seq(0,max(a$dbin),length.out=1000)
p<-predict(mod,newdata=data.frame(dbin=x))
lines(x,p,col='red')

dm<-data.frame(x=x,
               p=round(p,digits=2))
dc<-subset(dm,p==.37)

if(dim(dc)[1]>0){
dcall<-subset(dc,x==min(dc$x))
points(dcall$x,0,pch=18,col='red3',cex=3)
} else { dcall<-data.frame(x=-Inf,
                           p=NA)
}
d$k<-k#SLOPE: IF POSITIVE THEN DECORR IS BELOW DETECTION LIMIT
d$decor.ed<-round(dcall$x,digits=0)
return(d)
} else NULL
}

setwd(figsdir)
pdf('decorr_chl_spera.pdf')
par(mfrow=c(5,4),mar=c(2,1,1,1))
datcr<-ddply(dt,.(id),.fun=f,.progress='text')
dev.off()


datcr<-subset(datcr,decor!=-Inf)

mdat<-map_data('world')












##################################################################



#################################################################

centerfun<-function(d,var){
    #DETERMINE DAY OF MAX
    names(d)[1]<-'y'
    mod<-gam(y~s(day,bs='cs',k=5),data=d,gamma=1.4)
    pdat<-data.frame(xx=seq(min(d$day),max(d$day),1))
    pdat$p<-predict(mod,newdata=data.frame(day=pdat$xx))
    mday<-subset(pdat,p==max(pdat$p))$xx[1]
    #CENTER OF TIMING OF MAX
    shft<-(365/2)-mday
    d$day2<-d$day+shft
    d$day2<-ifelse(d$day2>365,d$day2-365,d$day2)
    d$day2<-ifelse(d$day2<0,d$day2+365,d$day2)
    d$day2<-ceiling(d$day2)
    d$var<-var
    dout<-na.omit(subset(d,select=c('year','month','day','day2','y','var','lon','lat')))
    names(dout)[5]<-var
    return(dout)
}
a1<-centerfun(subset(path,select=c('wind','day','year','month','lon','lat')),'wind')
#DON'T NEED TO CENTER SST BECAUSE ALREADY PEAKS IN FALL (DAY~230)
a2<-subset(path,select=c('sst','year','month','day','lon','lat'))
a2$day2<-a2$day
a2$var<-'sst'


phenfun<-function(d){
    print(unique(d$year))
    names(d)[1]<-'y'

if(dim(d)[1]>=10 & length(unique(d$month))>=8){
#FIT MODEL - THE WORKHORSE
mod<-gam(y~s(day2,bs='cc',k=5) + s(lon,lat,k=4),data=d,select=TRUE)
s<-summary(mod)

#FORMATS DATA TO PREDICT AT
    datout<-expand.grid(day2=seq(1,365,1))

datout$lon<-median(d$lon)
datout$lat<-median(d$lat)

#PREDICT
pred<-predict(mod,newdata=datout,na.action=na.pass,se.fit=TRUE)
datout$p<-pred$fit
datout$se<-pred$se.fit

datout$year<-unique(d$year)
datout$lat<-median(d$lat)
datout$lon<-median(d$lon)
datout$aic<-AIC(mod)
datout$dev<-s$dev.expl
datout$r2<-s$r.sq
datout$n<-dim(d)[1]

#GETS DATA AVAILABILITY FOR ID2
datout$nmonths<-length(unique(d$month))
datout$span<-max(d$day,na.rm=TRUE)-min(d$day,na.rm=TRUE)

return(datout)

}else NULL
}

ztemp<-ddply(subset(a2,select=c('sst','lon','lat','day2','month','year')),.(year),.fun=phenfun,.parallel=FALSE,.progress='text')
ztemp$var<-'sst'
zwnd<-ddply(subset(a1,select=c('wind','lon','lat','day2','month','year')),.(year),.fun=phenfun,.parallel=FALSE,.progress='text')
zwnd$var<-'wnd'


f<-function(d){
    lines(d$day,d$p)
}
plot(0,0,xlim=c(0,365),ylim=c(-0,15))
zz<-dlply(zwnd,.(year),.fun=f)
plot(0,0,xlim=c(0,365),ylim=c(5,10))
zz<-dlply(zwnd,.(year),.fun=f)


phenfun2<-function(d){
mx<-subset(d,p==max(d$p))
dout<-data.frame(year=unique(d$year))
dout$daypeak<-mx$day2[1]

if(unique(d$var)=='sst'){
    durpeak<-subset(d,p>=10)
} else {
    durpeak<-subset(d,p>=8)
}
dout$init<-subset(durpeak,day2==min(durpeak$day2))$day2[1]
dout$term<-subset(durpeak,day2==max(durpeak$day2))$day2[1]
dout$dur<-dout$term-dout$init
dout$amp<-max(d$p)-min(d$p)
dout$mx<-max(d$p)
dout$mn<-min(d$p)
return(dout)
}
vtemp<-ddply(ztemp,.(year),.fun=phenfun2)
vwnd<-ddply(zwnd,.(year),.fun=phenfun2)



setwd(datadir)
write.csv(vtemp,'vtemp.csv',row.names=FALSE)
write.csv(vwnd,'vwind.csv',row.names=FALSE)

plot(vtemp$year,vtemp$daypeak,pch=15)
plot(vwnd$year,vwnd$daypeak,pch=15)
plot(vwnd,pch=16)



