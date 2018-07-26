.libPaths()
#assign(".lib.loc", "C:/Rlibrary", envir = environment(.libPaths))
install.packages('cluster',dependencies=TRUE)
install.packages('rcompanion',dependencies=TRUE)
library(AICcmodavg)
library(rcompanion)
library(tidyr)
library(ggpubr)
library(pheatmap)
library(ellipse)
library(RColorBrewer)
library(DataCombine)
library(xts)
require(dlm)
library(maptools)
library(akima)
library(segmented)
library(splines)
library(strucchange)
library(data.table)
library(reshape2)
library(gplots)
library(cluster)
library(vegan)
library(ggplot2)
library(forecast)
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
library(rgdal)
library(car)
library(factoextra)
library(tidyr)
library(moments)
library(rgeos)

datadir<-'N://cluster_2017//scratch//spera//data//stagingdat'
codedir<-'N:/cluster_2017/scratch/spera/code'
figsdir<-"C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures"
setwd(codedir)
source('helper_functions.r')

setwd(datadir)
load("BoF_larval_herring_counts_spera.RData")

data<-subset(data,bathy<0 & spec_stage %in% c('L') & month>=8 & gearprocdesc=='Sawtooth Oblique')
data$stdno<-(data$totno/data$volm3)*abs(data$bathy)#NUMBER PER M2 - INTEGRATED OVER DEPTH
data$stdno<-ifelse(is.na(data$stdno)==TRUE,0,data$stdno)

#TAKE ONLY HERRING
dat<-subset(data,comname=='herring(atlantic)')

#BINS DATA TO 1/12 GRID RES
crds<-eacoordsfun(1/12)

#FUNCTION TO BIN DATA TO SPECIFIED RESOLUTION
dat<-cbind(dat,binfunction(subset(dat,select=c('lon','lat')),1/12))

#MERGE CELL COORDINATES TO DATA
dat<-merge(dat,crds,by=c('newarea5'),all.x=TRUE,all.y=FALSE)


mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast<-readOGR('N://data/shapefiles/naturalearthdata_ne_10m_land_poly',layer='ne_10m_land')#works for Alfcoast<-fortify(coast)
coast.mc<-crop(coast,extent(-180,180,-90,90),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fcoast.mc<-fortify(coast.mc)

a<-subset(dat,month>8 & month<12 & survey=='Bay of Fundy Herring')
a<-subset(a,lat>42.6)
a<-subset(a,!(lat<43 & lon< -66.5))
a<-subset(a,!(lat<43.5 & lon> -64.5))
a$lon<-round(a$lon,digits=1)
a$lat<-round(a$lat,digits=1)

find_hull<-function(df) df[chull(df$lon,df$lat),]
hulls<-find_hull(a)
tr<-as.matrix(subset(hulls,select=c('lon','lat')))
p<-Polygon(tr)
pdum<-Polygons(list(p),'poly')
plglrv<-SpatialPolygons(list(pdum),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plglrv<-erase(plglrv,coast.mc)

aa<-unique(subset(a,select=c('station','lon','lat')))
aa$lon<-round(aa$lon,digits=1)
aa$lat<-round(aa$lat,digits=1)
aa<-unique(aa)
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$lon,y=aa$lat)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=.03, byid=TRUE)
plglrv1 <- gUnaryUnion(buf1)
plglrv1<-erase(plglrv1,coast.mc)

setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/code/analysis/Herring_state_2018/SI/SI_data')
save(plglrv1,file='plglrv1.RData')
save(plglrv,file='plglrv.RData')




datt<-subset(dat,month>=9)

aa<-unique(subset(datt,select=c('lon','lat')))
find_hull<-function(df) df[chull(df$lon,df$lat),]
hulls<-find_hull(aa)
tr<-as.matrix(subset(hulls,select=c('lon','lat')))
p<-Polygon(tr)
pdum<-Polygons(list(p),'poly')
plgall<-SpatialPolygons(list(pdum),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#GETS MEAN NUMBER PER CELL
f<-function(d){
    return(data.frame(meanno=mean(d$stdno)))
}
tdat<-ddply(dat,.(clond,clatd,month,spec_stage),.fun=f,.progress='text')
tdat$lmeanno<-log10(tdat$mean+1)

d<-subset(tdat,select=c('lmeanno','clond','clatd'))
ttl<-'larv'
mn<-0
mx<-2.5
dg<-2
interp<-TRUE
K<-75

names(d)[1]<-'y'
d<-subset(d,is.na(y)==FALSE)
mod2<-gam(y~s(clond,clatd,k=K,bs='ts'),data=d,gamma=.25)
#LON/LAT TO INTERPOLATE AT
clonep<-seq(-70,-60,length.out=750)
clatep<-seq(42,46,length.out=750)
#PREDICTION DATA
pdat2<-expand.grid(clond=clonep,clatd=clatep)
pdat2$p<-predict(mod2,newdata=pdat2,type='response')
adat<-pdat2

#FUNCTION TO RETURN MAP
adat$p<-ifelse(adat$p<mn,mn,adat$p)
adat$p<-round(adat$p,digits=2)

#RETAIN ONLY PREDICTIONS OVER SPATIAL DOMAIN OF DATA
crds2<-SpatialPoints(data.frame(clond=adat$clond,clatd=adat$clatd),proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
adat$ov<-over(crds2,plgall)
adat<-subset(adat,ov==1)

setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/code/analysis/Herring_state_2018/SI/SI_data')
save(adat,file='adat.RData')


###################################################
#READ IN AND FORMAT WEBDROGUE PARTICLE TRACKING
setwd('N:/cluster_2017/scratch/spera/data/assorteddata')
sim1<-read.csv('german_bank_track.csv',header=FALSE,col.names=c('date','time','lon','lat','depth'))
sim1<-subset(sim1,is.na(lon)==FALSE & lon!='Top')
sim1$lon<-as.numeric(as.character(sim1$lon))
sim1$lat<-as.numeric(as.character(sim1$lat))
sim1<-subset(sim1,is.na(lat)==FALSE)

sim2<-read.csv('scotts_bay_track.csv',header=FALSE,col.names=c('date','time','lon','lat','depth'))
sim2<-subset(sim2,is.na(lon)==FALSE & lon!='Top')
sim2$lon<-as.numeric(as.character(sim2$lon))
sim2$lat<-as.numeric(as.character(sim2$lat))
sim2<-subset(sim2,is.na(lat)==FALSE)

sim1$id<-1
sim1$lon<-round(sim1$lon,digits=2)
sim1$lat<-round(sim1$lat,digits=2)
sim2$id<-1
sim2$lon<-round(sim2$lon,digits=2)
sim2$lat<-round(sim2$lat,digits=2)
f<-function(d){
return(data.frame(y=sum(d$id)))
}
sim1<-ddply(sim1,.(lon,lat),.fun=f,.progress='text')
sim2<-ddply(sim2,.(lon,lat),.fun=f,.progress='text')


setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/code/analysis/Herring_state_2018/SI/SI_data')
save(sim1,file='sim1.RData')
save(sim2,file='sim2.RData')


##############################################################

##############################################################
datadirout<-'N://cluster_2017//scratch//spera//data//dataoutput'
datadir1<-'N://cluster_2017//scratch//spera//data//stagingdat'
datadir<-"N:/cluster_2017/scratch/spera/data/finaldat_v2"
figsdir<-"C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures"
#figsdir<-"C://Users//copepod//Documents//aalldocuments//literature//research//active//SPERA//Figures"

mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast<-readOGR('N://data/shapefiles/naturalearthdata_ne_10m_land_poly',layer='ne_10m_land')#works for Alfcoast<-fortify(coast)
coast.mc<-crop(coast,extent(-180,180,-90,90),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/ESS_trophic_control/data')
plg<-readShapePoly('polygons_ecnasap.shp')#COMMAND TO READ BACK IN
plg<-subset(plg,region=='NS')

setwd(datadir)
rvwreadfun<-function(){
  setwd(datadir)
  rvw<-read.csv("herring_weights_RV_survey_spera_allar.csv",header=TRUE)
  rvw<-subset(rvw,month %in% c(6,7,8) & lat<46)
  rvw$strat<-as.character(rvw$strat)
  rvw<-subset(rvw,!(strat %in% c("5Z1","5Z2","5Z9",'493','494','','558','559','440','441','442','445','446','447','443','444','559','447','449','448','450','496','451')))
  return(rvw)
}
rvw<-rvwreadfun()
plg$stratum<-as.numeric(as.character(plg$stratum))
plgj<-subset(plg,stratum %in% c('493','494'))
plga<-subset(plg,!(stratum %in% c('493','494')))
plg2<-subset(plga,stratum>=480 & stratum<=495)

plga<-subset(plg,stratum %in% unique(rvw$strat))
plga<-subset(plga,!(stratum %in% c("5Z1","5Z2","5Z9",'493','494','','558','559','440','441','442','445','446','447','443','444','559','447','449','448','450','496','451')))
plga <- gUnaryUnion(plga)


#READ IN SHAPEFILE TO DEFINE AREA OF INTEREST
setwd('N:/cluster_2017/scratch/spera/data/shapefile_70perc')
plg70<-readShapePoly('plgdf70.shp')#COMMAND TO READ BACK IN
proj4string(plg70)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
buf1<-plgj
proj4string(buf1)<-CRS(mcrt)
buf1 <- gBuffer(buf1, width=.1, byid=TRUE)
plgj2 <- gUnaryUnion(buf1)
plgj2<-erase(plgj2,coast.mc)


xlm<-c(-68,-58.8)
ylm<-c(41,46.5)
atl<-getNOAA.bathy(lon1=xlm[1],lon2=xlm[2],lat1=ylm[1],lat2=ylm[2], resolution=4)

plgall<-subset(plg,stratum %in% unique(rvw$strat))
plgall <- gUnaryUnion(plgall)

nfo<-readOGR('N:/data/shapefiles/NAFO',layer='NAFOsubdiv_WGS84')#works for Alaska - most
nfo2<-subset(nfo,ET_ID %in% c('4wd','4we','4wf','4wg','4wh','4wj','4wk','4wl','4wm','4ww','4xs','4xl','4xm','4xn','4xo','4xp','4xq','4xr','4xs','4xx'))

setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/code/analysis/Herring_state_2018/SI/SI_data')
save(plga,file='plga.RData')
save(plgj2,file='plgj2.RData')
save(atl,file='atl.RData')


setwd("N:/cluster_2017/scratch/spera/data/finaldat_v2")
load('SPERA_andata_new.RData')
dats<-subset(data,select=c('year','herlrv.len','her.len.rv','her.ssbc','her.waa','her.rec1','her.prod','her.ajrat.rv','her.georng','her.totwgt.rv','her.metai.rv','herjuv.metai.rv','her.szpe.rv','her.cf.rv','herjuv.fmass.rv','her.fmass.rv','her.spcv','her.spnug','her.spvar'))

#EXPONENTIAL TRANSFORM TO 'NORMALIZE'
f<-function(x){transformTukey(x,plotit=FALSE)}
dats<-data.frame(cbind(subset(dats,select='year'),apply(dats[,2:dim(dats)[2]],2,f)))

f<-function(x){scale(x,center=TRUE,scale=TRUE)}
dats<-data.frame(cbind(subset(dats,select='year'),apply(dats[,2:dim(dats)[2]],2,f)))

library(tidyr)
dats<-dats %>% gather(var,y,-year)
dats<-subset(dats,year>=1965)

setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/code/analysis/Herring_state_2018/SI/SI_data')
save(dats,file='dats.RData')








###################################################################

setwd(datadir)
load('SPERA_andata_new.RData')
#load('SPERA_andata_spawnar.RData')
#load('SPERA_andata.RData')
nms<-names(data)
nms<-nms[grepl('her',nms)==TRUE]
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('\\.cat',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-subset(df,year>=1965)
df<-df[,!(names(df) %in% c('herjuv.prey.bof','her.jf.bof','her.preytot.bof','her.land','her.expr','herjuv.spvar','herjuv.spcv','herjuv.sprng','herjuv.spnug','her.szdiv.rv','her.dvm2.rv','herjuv.dvm2.rv','her.dur2.rv','herjuv.dur2.rv','her.durng.rv','herjuv.durng.rv','her.state','her.ssb','her.land.pct1','her.land.spdiv','her.land.sprich','herlrv.dep.bof','herlrv.spcv','herlrv.spnug','herlrv.sprng','herlrv.spvar','herjuv.dumx.rv'))]
df$her.dvm.rv<-ifelse(df$her.dvm.rv==Inf,NA,df$her.dvm.rv)
######

#CONVERT TO SHORT FORM, EXAMINE NORMALITY, TRANSFORM USING BEST EXPONENTIAL TRANSFORM
dfs<-df %>% gather(var,y,-year)

#TRANSFORM VARIABLES TO OPTIMIZE NORMALITY
f<-function(d){
if(unique(d$var)=='her.spnug'){
    NULL
    } else { d$y<-transformTukey(d$y,plotit=FALSE)
         }
  return(d)
}
dfs<-ddply(dfs,.(var),.fun=f)


f<-function(d){
  d$y<-(d$y-mean(d$y,na.rm=TRUE))/sd(d$y,na.rm=TRUE)
  return(d)
}
dfs<-ddply(dfs,.(var),.fun=f)

dfl<-spread(data=dfs,key=var, value=y)


setwd(figsdir)
dm<-read.csv('dm.csv',header=TRUE)
nms<-subset(dm,select=c('var','lbl.short'))
names(nms)<-c('var','lbl')
dff2<-dfs
dff2<-merge(dff2,nms,by=c('var'),all.x=TRUE,all.y=FALSE)
dff2<-unique(subset(dff2,select=c('var','lbl')))
dff2$lbl<-ifelse(dff2$var=='her.dvm.rv','Herring diurnal migration',as.character(dff2$lbl))
dff2<-spread(data=dff2,key=var, value=lbl)

#ALL METRICS
d<-subset(dfs,is.na(y)==FALSE)
mod<-gamm(y~as.factor(year),data=d,random=list(var=~1))
pdat<-data.frame(year=sort(unique(d$year)))
p<-predict(mod$gam,newdata=pdat,se.fit=TRUE)
pdat$her.state.all<-p$fit
pdat$her.state.all.se<-p$se.fit

#FINAL INDEX
d<-subset(dfs,var %in% c('herlrv.len','her.len.rv','her.cf.rv','herjuv.metai.rv','herjuv.fmass.rv','her.rec1','her.prod','her.metai.rv','her.spcv','her.spnug','her.spvar','her.ajrat.rv','her.fmass.rv','her.szpe.rv','her.ssbc','her.waa') & is.na(y)==FALSE)
mod<-gamm(y~as.factor(year),data=d,random=list(var=~1))
pdat0<-data.frame(year=sort(unique(d$year)))
p<-predict(mod$gam,newdata=pdat0,se.fit=TRUE)
pdat0$her.state.fin<-p$fit
pdat0$her.state.fin.se<-p$se.fit
pdat<-merge(pdat,pdat0,by=c('year'),all.x=TRUE,all.y=FALSE)


#INCLUDING INCREASING HEALTH INDICES
d<-subset(dfs,var %in% c('herlrv.len','her.len.rv','her.cf.rv','herjuv.metai.rv','herjuv.fmass.rv','her.rec1','her.prod','her.metai.rv','her.spcv','her.spnug','her.spvar','her.ajrat.rv','her.fmass.rv','her.szpe.rv','her.ssbc','her.waa','herjuv.totwgt.rv','her.totwgt.rv','her.georng') & is.na(y)==FALSE)
mod<-gamm(y~as.factor(year),data=d,random=list(var=~1))
pdat0<-data.frame(year=sort(unique(d$year)))
p<-predict(mod$gam,newdata=pdat0,se.fit=TRUE)
pdat0$her.state.19ind<-p$fit
pdat0$her.state.19ind.se<-p$se.fit
pdat<-merge(pdat,pdat0,by=c('year'),all.x=TRUE,all.y=FALSE)

#RV DATA OMITTED
d<-subset(dfs,var %in% c('herlrv.len','her.cf.rv','her.rec1','her.prod','her.ssbc','her.waa') & is.na(y)==FALSE)
mod<-gamm(y~as.factor(year),data=d,random=list(var=~1))
pdat0<-data.frame(year=sort(unique(d$year)))
p<-predict(mod$gam,newdata=pdat0,se.fit=TRUE)
pdat0$her.state.norv<-p$fit
pdat0$her.state.norv.se<-p$se.fit
pdat<-merge(pdat,pdat0,by=c('year'),all.x=TRUE,all.y=FALSE)

#ASSESSMENT DATA OMITTED
d<-subset(dfs,var %in% c('herlrv.len','her.len.rv','her.cf.rv','herjuv.metai.rv','herjuv.fmass.rv','her.metai.rv','her.spcv','her.spnug','her.spvar','her.ajrat.rv','her.fmass.rv','her.szpe.rv','her.waa') & is.na(y)==FALSE)
mod<-gamm(y~as.factor(year),data=d,random=list(var=~1))
pdat0<-data.frame(year=sort(unique(d$year)))
p<-predict(mod$gam,newdata=pdat0,se.fit=TRUE)
pdat0$her.state.noass<-p$fit
pdat0$her.state.noass.se<-p$se.fit
pdat<-merge(pdat,pdat0,by=c('year'),all.x=TRUE,all.y=FALSE)

a<-subset(pdat,select=c('her.state.all','her.state.fin','her.state.norv','her.state.noass'))
plot(a,pch=16)

save(pdat,file='C:\\Users\\sailfish\\Documents\\aalldocuments\\literature\\research\\active\\SPERA\\code\\analysis\\Herring_state_2018\\SI\\SI_data\\pdat.RData')




########################################################################

setwd("N:/cluster_2017/scratch/spera/data/finaldat_v2")
load('SPERA_andata_new.RData')
nms<-names(data)
nms<-nms[grepl('her',nms)==TRUE]
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('\\.cat',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-subset(df,year>=1965)
df<-df[,!(names(df) %in% c('herjuv.prey.bof','her.jf.bof','her.preytot.bof','her.land','her.expr','herjuv.spvar','herjuv.spcv','herjuv.sprng','herjuv.spnug','her.szdiv.rv','her.dvm2.rv','herjuv.dvm2.rv','her.dur2.rv','herjuv.dur2.rv','her.durng.rv','herjuv.durng.rv','her.state','her.ssb','her.land.pct1','her.land.spdiv','her.land.sprich','herlrv.dep.bof','herlrv.spcv','herlrv.spnug','herlrv.sprng','herlrv.spvar','herjuv.dumx.rv'))]
df$her.dvm.rv<-ifelse(df$her.dvm.rv==Inf,NA,df$her.dvm.rv)

df<-df[,!(names(df) %in% c('her.eggprod','her.tbio','her.totno.rv','herjuv.totno.rv'))]

#CONVERT TO SHORT FORM, EXAMINE NORMALITY, TRANSFORM USING BEST EXPONENTIAL TRANSFORM
dfs<-df %>% gather(var,y,-year)


#TRANSFORM VARIABLES TO OPTIMIZE NORMALITY
f<-function(d){
if(unique(d$var)=='her.spnug'){
    NULL
    } else { d$y<-transformTukey(d$y,plotit=FALSE)
         }
  return(d)
}
dfs<-ddply(dfs,.(var),.fun=f)

f<-function(d){
  d$y<-(d$y-mean(d$y,na.rm=TRUE))/sd(d$y,na.rm=TRUE)
  return(d)
}
dfs<-ddply(dfs,.(var),.fun=f)

dfl<-spread(data=dfs,key=var, value=y)

a<-subset(dfl,select=c('year','her.ssbc','her.fmass.rv'))
a<-slide(a,Var='her.ssbc',slideBy=6,NewVar='her.ssbc.t6')
a<-slide(a,Var='her.ssbc',slideBy=5,NewVar='her.ssbc.t5')
aa<-na.omit(subset(a,select=c('year','her.ssbc.t6','her.fmass.rv')))
mod<-lm(her.ssbc.t6~her.fmass.rv,data=a)

mod<-gls(her.ssbc.t6~her.fmass.rv,data=aa,method='ML',correlation=corCAR1(form = ~year))
mod2<-gls(her.ssbc.t6~her.fmass.rv,data=aa,method='ML')
pdat<-a
pdat<-subset(pdat,is.na(her.fmass.rv)==FALSE)
p<-predictSE(mod2,newdata=data.frame(her.fmass.rv=pdat$her.fmass.rv),se.fit=TRUE)
pdat$pmod<-p$fit
pdat$pmod.se<-p$se.fit
#pdat$pmod2<-predict(mod2,newdata=data.frame(her.fmass.rv=pdat$her.fmass.rv))
pdat2<-data.frame(year=seq(min(pdat$year)+6,max(pdat$year)+6,1),
 pmod.t6=slide(pdat,Var='pmod',slideBy=-6,NewVar='pmod.t6')$pmod.t6,
 pmod.se.t6=slide(pdat,Var='pmod.se',slideBy=-6,NewVar='pmod.se.t6')$pmod.se.t6)

pdat<-subset(pdat,select=c('year','her.ssbc'))
pdat2<-merge(pdat2,pdat,by=c('year'),all=TRUE)
save(pdat2,file='C:\\Users\\sailfish\\Documents\\aalldocuments\\literature\\research\\active\\SPERA\\code\\analysis\\Herring_state_2018\\SI\\SI_data\\pdat2.RData')
