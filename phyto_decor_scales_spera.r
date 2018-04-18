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


#datadir<-'/scratch/dboyce/chl_phenology/data'
datadir<-'/scratch/dboyce/spera/data'
figsdir<-'/scratch/dboyce/spera/figs'
codedir<-'/scratch/dboyce/spera/code'

setwd(codedir)
source('helper_functions.r')

setwd(datadir)
load('phytoplankton_biomass_all_spera.RData')


#nlat<-180*12
#nlon<-360*12

###################################################################
#GENERATES EQUAL-AREA GRID COORDINATES IN DEGREES AND PROJECTED UNITS
crds<-eacoordsfun(1/12)

#FUNCTION TO BIN DATA TO SPECIFIED RESOLUTION
dat<-cbind(data,binfunction(subset(data,select=c('lon','lat')),1/12))

#MERGE CELL COORDINATES TO DATA
data<-merge(dat,crds,by=c('newarea5'),all.x=TRUE,all.y=FALSE)

#AVERAGE VALUES PER CELL
avefun<-function(d){
    return(data.frame(chl=mean(d$chl,na.rm=TRUE),
                      clond=unique(d$clond),
                      clatd=unique(d$clatd),
                      cell=unique(d$newarea5),
                      day=unique(d$day),
                      year=unique(d$year),
                      month=unique(d$month),
                      db=unique(d$db),
                      bathy=mean(d$bathy,na.rm=TRUE),
                      dist=mean(d$dist,na.rm=TRUE),
                      tday=unique(d$tday),
                      mday=unique(d$mday),
                      cl=unique(d$cl),
                      tbin3=unique(d$tbin3),
                      tbin5=unique(d$tbin5)))
}
l<-dlply(data,.(db,newarea5,year,day),.fun=avefun,.progress='text')
data<-rbindlist(l)

save(data,file='phytoplankton_biomass_all_spera_gridded.RData')
load('phytoplankton_biomass_all_spera_gridded.RData')



data$tday<-(data$tday-min(data$tday))+1
data$id<-gsub(' ','',paste(data$cell,'_',data$db))






#########################################################################

#########################################################################

####PHENOLOGY
dat<-subset(data,year>=2000)

f1<-function(d0){

f<-function(d){
day<-unique(d$day)
d2<-subset(d0,select=c('day'))#GET ALL DAYS IN ID2 OBJECT
d2$day2<-d2$day-(day)
d2$day2<-ifelse(d2$day2<=0,d2$day2+(365),d2$day2)
    return(data.frame(day=day,
                   span=diff(range(min(d2$day2),max(d2$day2)))))
}
z<-ddply(d0,.(day),.fun=f)
z<-subset(z,span==max(z$span))
return(subset(z,day==min(z$day)))
}
dayshiftcell<-dlply(unique(subset(dat,select=c('id','day'))),.(id),.fun=f1,.progress='text',.parallel=FALSE)
dayshiftcell<-rbindlist(dayshiftcell)
dayshiftcell$cell<-sort(unique(dat$cell))


f2<-function(d){
d2<-subset(dayshiftcell,cell==unique(d$cell))
day<-d2$day#NUMER OF DAYS TO SHIFT
span<-d2$span
d$day2<-d$day-(day)
d$day2<-ifelse(d$day2<=0,d$day2+365,d$day2)
d$celloffset<-day
d$span<-span
return(d)
}
dat<-dlply(dat,.(cell),.fun=f2,.progress='text',.parallel=FALSE)
dat<-rbindlist(dat)



##############################################
#USE CALCULATIONS ABOVE TO ADJUST DAY VALUES FOR EACH CELL AND DB TO MAXIMIZE PHENOLOGY SPAN

d<-subset(dat,cell=="-64.042_45.375")

phenfun<-function(d){
#print(unique(d$cell))
options(warn=-1)

if(length(unique(d$year))>1 & dim(d)[1]>=10 & length(unique(d$month))>=8){
#FIT MODEL - THE WORKHORSE
mod<-try(gam(log10(chl+.1)~s(day2,bs='cc',k=6,id=1) + as.factor(year),data=d,select=TRUE,gamma=1.4))
s<-summary(mod)

#FORMATS DATA TO PREDICT AT
datout<-data.frame(day2=seq(1,365,1),
                   year=sort(unique(d$year),decreasing=TRUE)[1],
                   tday=sort(unique(d$tday))[1])

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

    if(length(unique(datout$p))<=2){nunique<-1
    } else nunique<-2
    
datout$cell<-unique(d$cell)
datout$db<-unique(d$db)
datout$newarea5<-unique(d$newarea5)
datout$lat<-unique(d$lat)
datout$lon<-unique(d$lon)
datout$aic<-AIC(mod)
datout$dev<-s$dev.expl
datout$r2<-s$r.sq
datout$celloffset<-unique(d$celloffset)
datout$nunique<-nunique
    
#BACK-TRANSFORMS DAYS
shiftfun<-function(d){
d$day<-d$day2+((unique(d$celloffset)))
d$day<-ifelse(d$day>365,d$day-365,d$day)
return(d)
}
datout<-shiftfun(datout)
    

#GETS DATA AVAILABILITY FOR ID2
datout$n.id<-dim(d)[1]
datout$nmonths.id<-length(unique(d$month))
datout$span.id<-max(d$day2,na.rm=TRUE)-min(d$day2,na.rm=TRUE)

return(datout)

}else NULL
}

z<-dlply(dat,.(cell),.fun=phenfun,.parallel=FALSE,.progress='text')
z<-rbindlist(z)
z<-subset(z,nunique==2)#OMIT INSTANCES WHERE PHENOLOGY IS FLAT


#ADD MONTHS
date<-strptime(gsub(' ','',paste(2000,'-',ceiling(z$day))),"%Y-%j")
z$month<-month(date)

dum<-data.frame(cell=sort(unique(dat$cell)),
                nmonth=tapply(dat$month,dat$cell,function(x) length(unique(x))))
dum<-subset(dum,nmonth>=11)
z<-subset(z,cell %in% unique(dum$cell))

#WRITE TO FILE
setwd(datadir)
write.csv(z,'phenology_spera_swfs.csv',row.names=FALSE)
gc()





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
pngMAP_df<- get_map(location = c(lon = -66.5, lat = 43.5), source = "google", zoom = 7,color='color',maptype='satellite',crop=FALSE)
p<-ggmap(pngMAP_df)


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






















