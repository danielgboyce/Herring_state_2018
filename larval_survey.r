library(rgdal)
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

datadir<-'N://cluster_2017//scratch//spera//data//stagingdat'
codedir<-'N:/cluster_2017/scratch/spera/code'
figsdir<-"C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures"
setwd(codedir)
source('helper_functions.r')

setwd(datadir)
load("BoF_larval_herring_counts_spera.RData")
load('BoF_larval_herring_lengths_spera.RData')


data<-subset(data,bathy<0 & spec_stage %in% c('L') & month>=8 & gearprocdesc=='Sawtooth Oblique')
data$stdno<-(data$totno/data$volm3)*abs(data$bathy)#NUMBER PER M2 - INTEGRATED OVER DEPTH
data$stdno<-ifelse(is.na(data$stdno)==TRUE,0,data$stdno)

data2<-subset(data2, spec_stage %in% c('L') & month>=8 & gearprocdesc=='Sawtooth Oblique')


#BINS INTO 3 YEAR INTERVALS
brks<-seq(1886,2015,3); lbls<-seq(min(brks)+1.5,max(brks)-1.5,3)
data$tbin3<-cut(data$year,breaks=brks,labels=lbls)

#BINS INTO 5 YEAR INTERVALS
brks<-seq(1885,2015,5);lbls<-seq(min(brks)+2.5,max(brks)-2.5,5)
data$tbin5<-cut(data$year,breaks=brks,labels=lbls)

#BINS INTO 10 YEAR INTERVALS
brks<-seq(1885,2015,10);lbls<-seq(min(brks)+5,max(brks)-5,10)
data$tbin10<-cut(data$year,breaks=brks,labels=lbls)
#data<-data[,-1]

#TAKE ONLY HERRING
dat<-subset(data,comname=='herring(atlantic)')
dat2<-subset(data2,comname=='herring(atlantic)')

#BINS DATA TO 1/12 GRID RES
crds<-eacoordsfun(1/12)

#FUNCTION TO BIN DATA TO SPECIFIED RESOLUTION
dat<-cbind(dat,binfunction(subset(dat,select=c('lon','lat')),1/12))
dat2<-cbind(dat2,binfunction(subset(dat2,select=c('lon','lat')),1/12))

#MERGE CELL COORDINATES TO DATA
dat<-merge(dat,crds,by=c('newarea5'),all.x=TRUE,all.y=FALSE)
dat2<-merge(dat2,crds,by=c('newarea5'),all.x=TRUE,all.y=FALSE)

f<-function(d){
    return(data.frame(len=weighted.mean(d$length,w=d$clen)))
}
ot<-ddply(dat2,.(clond,clatd,yday),.fun=f,.progress='text')
ot2<-subset(ot,clond> -65.5 & clatd<45)

plot(ot$clond,ot$clatd,pch=16)
plot(dat2$lon,dat2$lat,pch=16)

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
plg <- gUnaryUnion(buf1)
plg<-erase(plg,coast.mc)

library(marmap)
xlm<-c(-68,-58.8)
ylm<-c(41,46.5)
atl<-getNOAA.bathy(lon1=xlm[1],lon2=xlm[2],lat1=ylm[1],lat2=ylm[2], resolution=4)

setwd(figsdir)
pdf('sampling_spatial_domains_lrv.pdf.pdf',height=4,width=7)
par(mfrow=c(1,3),mar=c(0,0,0,0),oma=c(4,1,4,1))
cl3<-'green3'
xlm<-c(-68.2,-58.5)
ylm<-c(41.5,46.5)
lw<-.5
mpcl<-'gray80'
map('worldHires',col=mpcl,fill=TRUE,border=NA,xlim=xlm,ylim=ylm)
plot(atl, deep=-200, shallow=-20, step=200,add=TRUE,col=c(mpcl,'gray60'),lty=c(1,1),lwd=lw)
plot(plglrv,add=TRUE,col=alpha(cl3,.4),border=cl3,lwd=.001,xlim=xlm,ylim=ylm,lwd=1)
plot(plg,add=TRUE,col=alpha(cl3,.7),xlim=xlm,ylim=ylm,border=NA)
box()

map('worldHires',col=mpcl,fill=TRUE,border=NA,xlim=xlm,ylim=ylm)
plot(atl, deep=-200, shallow=-20, step=200,add=TRUE,col=c(mpcl,'gray60'),lty=c(1,1),lwd=lw)
plot(plglrv,add=TRUE,col=alpha(cl3,.4),border=cl3,lwd=.001,xlim=xlm,ylim=ylm,lwd=1)
plot(plg,add=TRUE,col=alpha(cl3,.7),xlim=xlm,ylim=ylm,border=NA)
box()
dev.off()

#GETS MEAN NUMBER PER CELL
f<-function(d){
    return(data.frame(meanno=mean(d$stdno)))
}
tdat<-ddply(dat,.(clond,clatd,month,spec_stage),.fun=f,.progress='text')
tdat$lmeanno<-log10(tdat$mean+1)



mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast<-readOGR('N://data/shapefiles/naturalearthdata_ne_10m_land_poly',layer='ne_10m_land')#works for Alfcoast<-fortify(coast)
coast.mc<-crop(coast,extent(-180,180,-90,90),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fcoast.mc<-fortify(coast.mc)


ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }


adat<-subset(tdat,month==6,select=c('lmeanno','clond','clatd'))
ttl<-'ED'
mn<-0
mx<-2.5
dg<-2

#FUNCTION TO RETURN MAP
pltfun<-function(adat,ttl,mx,dg){
if(dim(adat)[1]>25){
mn<-0
nm<-names(adat)[1]
names(adat)[1]<-'y'
adat$y<-ifelse(adat$y>mx,mx,adat$y)
adat$y<-ifelse(adat$y<mn,mn,adat$y)
adat$y<-round(adat$y,digits=2)
adat<-rbind.fill(data.frame(y=c(seq(0,mx,length.out=1000))),adat)
a<-adat

n<-21
brks<-seq((min(adat$y,na.rm=TRUE)-.01),max(adat$y,na.rm=TRUE)+.01,length.out=n)
brks2<-round(seq((min(adat$y,na.rm=TRUE)-.01),max(adat$y,na.rm=TRUE)+.01,length.out=n),digits=dg)
a$ycat<-cut(a$y,breaks=brks)
lbls<-sort(unique(a$ycat))
lbls2<-sort(unique(cut(a$y,breaks=brks2)))
cls<-matlab.like(length(lbls))
#cls<-colorRampPalette(c('magenta4','blue3','green','yellow','red3'))
#cls<-(cls(length(lbls)))

return(
ggplot()+
geom_tile(data=a, aes(x=clond, y=clatd,fill=ycat),col=NA,size=.0001)+
scale_fill_manual(breaks=as.character(lbls),values=cls,labels=lbls2,na.value="transparent",guide=guide_legend(title=''))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.95,.25),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.06, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-70,-59,1),labels=as.character(seq(-70,-59,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(42,46,1),labels=as.character(seq(42,46,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(41,46),xlim=c(-69,-59.5))+
    xlab('')+
    ylab('')+
       labs(title = '',
       subtitle = ttl,
       caption = '')

)
} else NULL
}
p3<-pltfun(subset(tdat,spec_stage=='L' &month==3,select=c('lmeanno','clond','clatd')),'Month=3',3,1)
p4<-pltfun(subset(tdat,spec_stage=='L' &month==4,select=c('lmeanno','clond','clatd')),'Month=4',3,1)
p8<-pltfun(subset(tdat,spec_stage=='L' &month==8,select=c('lmeanno','clond','clatd')),'Month=8',3,1)
p10<-pltfun(subset(tdat,spec_stage=='L' &month==10,select=c('lmeanno','clond','clatd')),'Month=10',3,1)
p11<-pltfun(subset(tdat,spec_stage=='L' &month==11,select=c('lmeanno','clond','clatd')),'Month=11',3,1)
p0<-pltfun(subset(tdat,spec_stage=='L',select=c('lmeanno','clond','clatd')),'Total log no',3,1)
p00<-pltfun(subset(tdat,spec_stage=='L',select=c('meanno','clond','clatd')),'Total No',3,1)

setwd(figsdir)
pdf('herlarv_abund_map_spera.pdf',height=18,width=13)
grid.arrange(p3,p4,p8,p10,p11,p0,p00,ncol=2)
dev.off()







#GETS MEAN NUMBER PER CELL
f<-function(d){
    return(data.frame(meanlnth=weighted.mean(d$length,w=d$clen)))
}
tdat<-ddply(dat2,.(clond,clatd,month,spec_stage),.fun=f,.progress='text')



p3<-pltfun(subset(tdat,spec_stage=='L' &month==3,select=c('meanlnth','clond','clatd')),'Month=3',35,1)
p4<-pltfun(subset(tdat,spec_stage=='L' &month==4,select=c('meanlnth','clond','clatd')),'Month=4',35,1)
p8<-pltfun(subset(tdat,spec_stage=='L' &month==8,select=c('meanlnth','clond','clatd')),'Month=8',35,1)
p10<-pltfun(subset(tdat,spec_stage=='L' &month==10,select=c('meanlnth','clond','clatd')),'Month=10',35,1)
p11<-pltfun(subset(tdat,spec_stage=='L' &month==11,select=c('meanlnth','clond','clatd')),'Month=11',35,1)
p0<-pltfun(subset(tdat,spec_stage=='L', select=c('meanlnth','clond','clatd')),'Total log no',35,1)

setwd(figsdir)
pdf('herlarv_length_map_spera.pdf',height=18,width=13)
grid.arrange(p3,p4,p8,p10,p11,p0,ncol=2)
dev.off()


p3<-pltfun(subset(tdat,spec_stage=='J' &month==3,select=c('meanlnth','clond','clatd')),'Month=3',130,1)
p4<-pltfun(subset(tdat,spec_stage=='J' &month==4,select=c('meanlnth','clond','clatd')),'Month=4',130,1)
p11<-pltfun(subset(tdat,spec_stage=='J' &month==11,select=c('meanlnth','clond','clatd')),'Month=11',130,1)
p0<-pltfun(subset(tdat,spec_stage=='J', select=c('meanlnth','clond','clatd')),'Total log no',130,1)

setwd(figsdir)
pdf('herlarv_juv_length_map_spera.pdf',height=10,width=13)
grid.arrange(p3,p4,p11,p0,ncol=2)
dev.off()








###########################################
#POLYGON OVER ALL POINTS
#SPAWNING ON GERMAN BANK AUG-SEPT, EARLIER ON LURCHER
datt<-subset(dat,month>=9)
datt2<-subset(dat,month<9)

aa<-unique(subset(datt,select=c('lon','lat')))
find_hull<-function(df) df[chull(df$lon,df$lat),]
hulls<-find_hull(aa)
tr<-as.matrix(subset(hulls,select=c('lon','lat')))
p<-Polygon(tr)
pdum<-Polygons(list(p),'poly')
plgall<-SpatialPolygons(list(pdum),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#GETS MEAN NUMBER PER CELL
f<-function(d){    return(data.frame(meanno=mean(d$stdno)))}
tdat<-ddply(subset(datt,spec_stage=='L'),.(clond,clatd),.fun=f,.progress='text')
tdat$lmeanno<-log10(tdat$mean+1)

p1<-pltfun(subset(tdat,select=c('lmeanno','clond','clatd')),'Total no',2.75,1)
p1

p1+
geom_point(aes(lon,lat),data=aa,col='green')+
geom_polygon(aes(long,lat),data=plgall,alpha=.5,col='gold4')

z<-sdum-coast.mc









d<-subset(tdat,select=c('lmeanno','clond','clatd'))
ttl<-'larv'
mn<-0
mx<-2.5
dg<-2
interp<-TRUE
K<-80

#pfunint<-function(d,mn,mx,ttl,interp,K){
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
#adat$p<-ifelse(adat$p>mx,mx,adat$p)
adat$p<-ifelse(adat$p<mn,mn,adat$p)
adat$p<-round(adat$p,digits=2)

#RETAIN ONLY PREDICTIONS OVER SPATIAL DOMAIN OF DATA
crds2<-SpatialPoints(data.frame(clond=adat$clond,clatd=adat$clatd),proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
adat$ov<-over(crds2,plgall)
adat<-subset(adat,ov==1)

a<-adat
#SOME NONSENSE TO GET THE COLOUR SCALE RIGHT
aa<-data.frame(p=seq((min(adat$p,na.rm=TRUE)),max(abs(adat$p),na.rm=TRUE),length.out=100))
a<-rbind.fill(a,aa)

n<-20
brks<-seq(min(a$p,na.rm=TRUE)-.01,max(abs(a$p),na.rm=TRUE)+.01,length.out=n)
brks2<-round(seq(min(a$p,na.rm=TRUE)+.01,max(abs(a$p),na.rm=TRUE)+.01,length.out=n),digits=2)
a$ycat<-cut(a$p,breaks=brks)
lbls<-sort(unique(a$ycat))
lbls2<-sort(unique(cut(a$p,breaks=brks2)))
#cls<-colorRampPalette(c('magenta4','blue3','green','yellow','red3'))
cls<-matlab.like(length(lbls))
#cls<-colorRampPalette(brewer.pal(name='RdYlBu',8))
#cls<-rev(cls(length(lbls)))

p1int<-(
ggplot()+
geom_tile(data=a, aes(x=clond, y=clatd,fill=ycat))+
scale_fill_manual(breaks=as.character(lbls),values=cls,labels=lbls2,na.value="transparent",guide=guide_legend(title=paste(ttl)))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-69,-62,1),labels=as.character(seq(-69,-62,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(42,46,1),labels=as.character(seq(42,46,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(42,max(adat$clatd)),xlim=c(-69,-60))+
    xlab('')+
    ylab('')
)

p1int


p1int

setwd(figsdir)
pdf('herlarv_abund_map_interp_spera3.pdf',width=6,height=8)
grid.arrange(p1int,p1int,ncol=1)
dev.off()



#CREATE POLYGONS WHERE HERRING LARVAL ABUNDANCES ARE ABOVE UPPER 90 OR 95 PERCENTILE
library(rgeos)
qtl95<-quantile(subset(a)$p,.95)
aa<-subset(a,p>qtl95 & is.na(clond)==FALSE & clatd>42.8& clond>=-67 & clond< -65)
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$clond,y=aa$clatd)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=.01, byid=TRUE)
plg95 <- gUnaryUnion(buf1)
#buf <-  gBuffer(buf2, width=0)
plg95<-erase(plg95,coast.mc)


qtl85<-quantile(subset(a)$p,.85)
aa<-subset(a,p>qtl85 & is.na(clond)==FALSE & clond< -65.9 & clatd>42.75& clond>=-67.2)
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$clond,y=aa$clatd)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=.01, byid=TRUE)
plg85 <- gUnaryUnion(buf1)
plg85<-erase(plg85,coast.mc)

qtl80<-quantile(subset(a)$p,.80)
aa<-subset(a,p>qtl80 & is.na(clond)==FALSE & clond< -65.9 & clatd>42.75& clond>=-67.2 & !(clatd>44.5 & clond> -66))
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$clond,y=aa$clatd)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=.01, byid=TRUE)
plg80 <- gUnaryUnion(buf1)
plg80<-erase(plg80,coast.mc)


qtl75<-quantile(subset(a)$p,.75)
aa<-subset(a,p>qtl75 & is.na(clond)==FALSE & clond< -65.7 & clatd>42.75& clond>=-67.2 & !(clatd>44.5 & clond> -66.1))
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$clond,y=aa$clatd)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=.01, byid=TRUE)
plg75 <- gUnaryUnion(buf1)
plg75<-erase(plg75,coast.mc)
p175<-p1int+
    geom_polygon(aes(x=long,y=lat,color='green'),data=plg75,alpha=.75,col=NA)


qtl70<-quantile(subset(a)$p,.70)
aa<-subset(a,p>qtl70 & is.na(clond)==FALSE & clond< -65.7 & clatd>42.75& clond>=-67.2 & !(clatd>44.5 & clond> -66.1))
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$clond,y=aa$clatd)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=.02, byid=TRUE)
plg70 <- gUnaryUnion(buf1)
plg70<-erase(plg70,coast.mc)
p170<-p1int+
    geom_polygon(aes(x=long,y=lat,color='green'),data=plg70,alpha=.70,col=NA)
p170


#THIS VERSION HAS A LARGER BUFFER AROUND IT
qtl70<-quantile(subset(a)$p,.70)
aa<-subset(a,p>qtl70 & is.na(clond)==FALSE & clond< -65.75 & clatd>42.75& clond>=-68 & !(clatd>44.5 & clond> -65) & !(clatd<=43.1 & clond< -66.75))
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$clond,y=aa$clatd)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=.5, byid=TRUE)
plg70 <- gUnaryUnion(buf1)
plg70<-erase(plg70,coast.mc)
p170<-p1int+
    geom_polygon(aes(x=long,y=lat,color='green'),data=plg70,alpha=.70,col=NA)
p170




qtl70<-quantile(subset(a)$p,.75)
aa<-subset(a,p>qtl70 & is.na(clond)==FALSE & clond< -65.7 & clatd>42.75& clond>=-67.2 & !(clatd>44.5 & clond> -66.1))
#aa<-subset(a,p>qtl70 & is.na(clond)==FALSE)
pts0<-SpatialPoints(as.matrix(data.frame(x=aa$clond,y=aa$clatd)),CRS(mcrt))
buf1 <- gBuffer(pts0, width=0, byid=TRUE)
plg70 <- gUnaryUnion(buf1)

buf1<-spTransform(buf1, CRS( "+init=epsg:3347"))
buf1<-gBuffer(buf1,byid=TRUE,width=1)
buf1<-spTransform(buf1, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#plg70<-erase(plg70,coast.mc)
p170<-p1int+
    geom_polygon(aes(x=long,y=lat,color='green'),data=plg70,alpha=.75,col=NA)

p170


p175<-p1int+
    geom_polygon(aes(x=long,y=lat,color='green'),data=plg75,alpha=.75,col=NA)
p185<-p1int+
    geom_polygon(aes(x=long,y=lat),data=plg85,alpha=.75,col=NA)
p180<-p1int+
    geom_polygon(aes(x=long,y=lat),data=plg80,alpha=.75,col=NA)

library(sp)
library(maptools)
setwd('/scratch/dboyce/spera/data/shapefile_80perc')
lu<-data.frame(percent=rep('eighty',length(plg80)))
plgdf80 <- SpatialPolygonsDataFrame(plg80, lu)
proj4string(plgdf80)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writePolyShape(plgdf80,'plgdf80.shp')


setwd('/scratch/dboyce/spera/data/shapefile_75perc')
lu<-data.frame(percent=rep('seventyfive',length(plg75)))
plgdf75 <- SpatialPolygonsDataFrame(plg75, lu)
proj4string(plgdf75)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writePolyShape(plgdf75,'plgdf75.shp')


setwd('N:/cluster_2017/scratch/spera/data/shapefile_70perc')
lu<-data.frame(percent=rep('seventy',length(plg70)))
plgdf70 <- SpatialPolygonsDataFrame(plg70, lu)
proj4string(plgdf70)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writePolyShape(plgdf70,'plgdf70.shp')





setwd(figsdir)
pdf('herlarv_abund_interp_map_spera.pdf',width=6,height=10)
grid.arrange(p185,p180,p175,ncol=1)
dev.off()

setwd(datadir)
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



setwd(figsdir)
pdf('herlarv__webd_map_spera.pdf',width=6,height=5)
ggplot()+
geom_point(aes(lon,lat),data=sim1,col='green',alpha=.15,size=.5)+
geom_point(aes(lon,lat),data=sim2,col='purple',alpha=.15,size=.5)+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6),plot.background = element_rect(fill = 'green', colour = 'red'))+
scale_x_continuous(expand=c(0,0),breaks=seq(-70,-62,1),labels=as.character(seq(-70,-62,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(43,46,1),labels=as.character(seq(43,46,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(42.75,max(adat$clatd)),xlim=c(-68,-64))+
    xlab('')+
    ylab('')

ggplot()+
geom_polygon(aes(x=long,y=lat,color='red'),data=plg75,alpha=.75,col=NA,fill='darkred')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6),plot.background = element_rect(fill = 'green', colour = 'red'))+
scale_x_continuous(expand=c(0,0),breaks=seq(-70,-62,1),labels=as.character(seq(-70,-62,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(43,46,1),labels=as.character(seq(43,46,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(42.75,max(adat$clatd)),xlim=c(-68,-64))+
    xlab('')+
    ylab('')
dev.off()


setwd(figsdir)
pdf('herlarv_abund_interp_webd_75pct_map_spera.pdf',width=6,height=5)
p175+
geom_point(aes(lon,lat),data=sim1,col='green',alpha=.15,size=.5)+
geom_point(aes(lon,lat),data=sim2,col='purple',alpha=.15,size=.5)
dev.off()








d<-subset(datt,tbin5==1972.5)
K<-50
mx<-15
d<-subset(datt,spec_stage=='L' & tbin5==tb[1])
K<-60
mx<-200
plgn<-plg75

###MAKES MAP OF PREDICTED PATTERN OF LARVAL HERRING ABUNDANCE FOR 5 YEAR INTERVALS
mapfun<-function(d,K,mx,plgn){

f<-function(dd){
    return(data.frame(meanno=mean(dd$stdno)))
}
dd<-ddply(d,.(clond,clatd),.fun=f,.progress='text')
dd$lmeanno<-log10(dd$mean+1)

dd<-subset(dd,select=c('lmeanno','clond','clatd'))
ttl<-as.character(unique(d$tbin5))
dg<-1

crd<-unique(subset(d,select=c('clond','clatd')))

    names(dd)[1]<-'y'
    dd<-subset(dd,is.na(y)==FALSE)
    mod2<-gam(y~s(clond,clatd,k=K,bs='ts'),data=dd,gamma=.25)
#    mod2<-gam(stdno~s(clond,clatd,k=K,bs='ts'),data=d,gamma=.25,family=nb)
    #LON/LAT TO INTERPOLATE AT
    clonep<-seq(-70,-62,length.out=600)
    clatep<-seq(42,46,length.out=600)
    #PREDICTION DATA
    pdat2<-expand.grid(clond=clonep,clatd=clatep)
    pdat2$p<-predict(mod2,newdata=pdat2,type='response')
    adat<-pdat2

crds2<-SpatialPoints(data.frame(clond=adat$clond,clatd=adat$clatd),proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

adat$ov<-over(crds2,plgall)
adat<-subset(adat,ov==1)

adat$p<-adat$p*100
print(summary(adat$p))
#FUNCTION TO RETURN MAP
adat$p<-ifelse(adat$p>mx,mx,adat$p)
adat$p<-ifelse(adat$p<0,0,adat$p)
adat$p<-round(adat$p,digits=dg)
a<-adat
#SOME NONSENSE TO GET THE COLOUR SCALE RIGHT
aa<-data.frame(p=seq((min(adat$p,na.rm=TRUE)),max(abs(adat$p),na.rm=TRUE),length.out=100))
a<-rbind.fill(a,aa)

n<-15
brks<-seq(min(a$p,na.rm=TRUE)-.01,max(abs(a$p),na.rm=TRUE)+.01,length.out=n)
brks2<-round(seq(min(a$p,na.rm=TRUE)+.01,max(abs(a$p),na.rm=TRUE)+.01,length.out=n),digits=2)
a$ycat<-cut(a$p,breaks=brks)
lbls<-sort(unique(a$ycat))
lbls2<-sort(unique(cut(a$p,breaks=brks2)))
#cls<-colorRampPalette(c('magenta4','blue3','green','yellow','red3'))
cls<-matlab.like(length(lbls))

return(
ggplot()+
geom_tile(data=a, aes(x=clond, y=clatd,fill=ycat))+
geom_tile(data=crd, aes(x=clond, y=clatd,alpha=0,fill=NA,na.value=NA,color='black'),color=NA)+
geom_polygon(aes(x=long,y=lat),data=plgn,alpha=.5,col=NA)+
scale_fill_manual(breaks=as.character(lbls),values=cls,labels=lbls2,na.value="transparent",guide=guide_legend(title=paste(ttl)))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.1, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-70,-64,1),labels=as.character(seq(-70,-64,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(43,46,1),labels=as.character(seq(43,46,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(42.75,max(adat$clatd)),xlim=c(-68,-64))+
    xlab('')+
    ylab('')
)

}


tb<-sort(unique(datt$tbin5))
p1<-mapfun(subset(datt,spec_stage=='L' & tbin5==tb[1]),60,200,plg75)
p2<-mapfun(subset(datt,spec_stage=='L' &tbin5==tb[2]),60,200,plg75)
p3<-mapfun(subset(datt,spec_stage=='L' &tbin5==tb[3]),60,200,plg75)
p4<-mapfun(subset(datt,spec_stage=='L' &tbin5==tb[4]),60,200,plg75)
p5<-mapfun(subset(datt,spec_stage=='L' &tbin5==tb[5]),60,200,plg75)
p6<-mapfun(subset(datt,spec_stage=='L' &tbin5==tb[6]),60,200,plg75)
p7<-mapfun(subset(datt,spec_stage=='L' &tbin5==tb[7]),60,200,plg85)

setwd(figsdir)
pdf('herlarv_abund_interp_map_timeseries_85overlay_spera.pdf',width=15,height=10)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,ncol=3)
dev.off()


p1<-mapfun(subset(datt,tbin5==tb[1]),60,200,plg80)
p2<-mapfun(subset(datt,tbin5==tb[2]),60,200,plg80)
p3<-mapfun(subset(datt,tbin5==tb[3]),60,200,plg80)
p4<-mapfun(subset(datt,tbin5==tb[4]),60,200,plg80)
p5<-mapfun(subset(datt,tbin5==tb[5]),60,200,plg80)
p6<-mapfun(subset(datt,tbin5==tb[6]),60,200,plg80)
p7<-mapfun(subset(datt,tbin5==tb[7]),60,200,plg80)

setwd(figsdir)
pdf('herlarv_abund_interp_map_timeseries_80overlay_spera.pdf',width=15,height=10)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,ncol=3)
dev.off()








######################################################
#SUMMARY PLOTS OF HERRING ASSESSMENTS
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


plot(her$year,her$sumf)
plot(her$year,her$mnf)
par(mfrow=c(2,2))
plot(her$year,her$s,type='l')
plot(her$year,her$ssb,type='l')
plot(her$year,her$r,type='l')
plot(her$year,her$larv,type='l')

d<-subset(her,select=c('r','year'))
cl<-'red'


library(segmented)
library(RColorBrewer)
cor(her$ssb,her$f,use='pairwise.complete.obs')
plot(diff(her$ssb),diff(her$f))
cor(diff(her$ssb),diff(her$f),use='pairwise.complete.obs')
plot(her$year,rescale(her$ssb,newrange=c(0,100)),pch=16,col='black')
points(her$year,rescale(her$f,newrange=c(0,100)),pch=16,col='red')

plot(her$year,rescale(her$ssb,newrange=c(0,100)),pch=16,col='black',type='l')
points(her$year,rescale(her$f,newrange=c(0,100)),pch=16,col='red',type='l')
dm<-na.omit(subset(her,select=c('f','ssb')))
ccf(as.ts(dm$f),as.ts(dm$ssb))

dm<-na.omit(subset(her,select=c('s','ssb')))
ccf(as.ts(dm$s),as.ts(dm$ssb),ylim=c(-1,1))

plot(her$ssb,her$f)
plot(her$ssb,log(her$f))
hist(her$f,breaks=100)

setwd(figsdir)
pdf('herring_stock_breakpointtrends.pdf',height=9,width=8)
par(mfrow=c(3,2),mar=c(2,4,1,1))

d<-subset(her,select=c('r','year'))
ylabb<-'Tons [x1000]'

segfun<-function(d,ylabb){
cl<-'black'
d<-na.omit(d)
nm<-names(d)[1]
names(d)[1]<-'y'
if(ylabb=='Tons [x1000]'){d$y<-d$y/1000
} else NULL
par(mfrow=c(2,2))
AIC(lin.mod,lin.mod2,segmod)

lin.mod <- lm(y~year,data=d)
lin.mod2 <- lm(y~year + I(year^2) + I(year^3),data=d)
s<-summary(lin.mod)
segmod <- segmented(lin.mod, seg.Z = ~year,method='ML')
psi<-segmod$psi[2]
s2<-summary(segmod)
summary(segmod)
print(segmod$psi)
out<-data.frame(slope(segmod)$year)
pdat<-data.frame(xx=seq(min(d$year),max(d$year),length.out=100))
plot(d$year,d$y,las=1,xlab='Year',ylab=ylabb,pch=15)

if((AIC(lin.mod)-AIC(segmod))>0){
p<-predict(segmod,newdata=data.frame(year=pdat$xx),se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat$upr<-pdat$p+(1.96*pdat$se)
pdat$lwr<-pdat$p-(1.96*pdat$se)
abline(v=psi,lty=3)
text(psi,min(d$y),gsub(' ','',paste('(',round(psi,digits=0),')')),pos=2)
} else {
p<-predict(lin.mod,newdata=data.frame(year=pdat$xx),se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat$upr<-pdat$p+(1.96*pdat$se)
pdat$lwr<-pdat$p-(1.96*pdat$se)
}
polygon(c(pdat$xx,pdat$xx[length(pdat$xx):1]),c(pdat$upr,pdat$lwr[length(pdat$lwr):1]),col=alpha('lightgray',.75),border=NA)
lines(pdat$xx,pdat$p,col=cl)
points(d$year,d$y,pch=15,col=cl)
legend('topleft',paste(c(paste('r2 (lin) =',round(s$adj.r.squared,digits=2)),paste('r2 (seg) =',round(s2$adj.r.squared,digits=2)))),bty='n')
out$psi<-segmod$psi[2]
out$var<-nm
return(out)
}
s.ssb<-segfun(subset(her,select=c('ssb','year')),'Tons [x1000]')
s.larv<-segfun(subset(her,select=c('larv','year')),'Larvae [n m3]')
s.f<-segfun(subset(her,select=c('f','year')),'Average F')
s.s<-segfun(subset(her,select=c('s','year')),'Larval survivorship')
s.r<-segfun(subset(her,select=c('r','year')),'Recruitment')



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


setwd(figsdir)
pdf('her_ssb_timetrend_timeline.pdf',height=5,width=10)
a<-unique(subset(s.wgt,segtrue==1,select=c('psi','var')))
a<-data.frame(psi=mean(a$psi),
              var='weight')
psidat<-unique(subset(rbind(s.ssb,s.larv,s.f),select=c('psi','var')))
psidat<-rbind(psidat,a)

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
rect(1965,-100,1988,7e+06,col=alpha('green3',.1),border=NA)
rect(1988,-100,1995,7e+06,col=alpha('gold',.1),border=NA)
rect(1995,-100,2006,7e+06,col=alpha('blue',.1),border=NA)
cl<-'red3'
rug(psidat$psi,lwd=2,col=cl)
text(psidat$psi[1],0,'SSB collapse starts',srt=40,adj=0,col=cl,cex=.75)
text(psidat$psi[2],0,'Herring larvae decline starts',srt=40,adj=0,col=cl,cex=.75)
text(psidat$psi[3],0,'F increase starts',srt=40,adj=0,col=cl,cex=.75)
text(psidat$psi[4],0,'Herring weights start to decline',srt=40,adj=0,col=cl,cex=.75)
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
