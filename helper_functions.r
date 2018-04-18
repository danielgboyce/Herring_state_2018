#GENERATES EQUAL-AREA GRID COORDINATES IN DEGREES AND PROJECTED UNITS
eacoordsfun<-function(stp){
#SET UP FULL GLOBAL COORDINATES
crds<-expand.grid(lon=seq(-180+(stp/2),180-(stp/2),by=stp),lat=seq(-90+(stp/2),90-(stp/2),by=stp))

#BINS LON INTO INTERVALS
brks<-seq(-180,180,by=stp)#VALUES TO DEFINE GRID EDGES
lbls<-seq(min(brks)+(stp/2),max(brks)-(stp/2),by=stp)#LABELS FOR LON GRID
crds$lonref<-round(as.numeric(as.character(cut(crds$lon,breaks=brks,labels=lbls))),digits=4)

#BINS LAT INTO INTERVALS
brks<-seq(-90,90,by=stp)#VALUES TO DEFINE GRID EDGES
lbls<-seq(min(brks)+(stp/2),max(brks)-(stp/2),by=stp)#LABELS FOR LON GRID
crds$latref<-round(as.numeric(as.character(cut(crds$lat,breaks=brks,labels=lbls))),digits=4)
crds$newarea5<-gsub(' ','',paste(crds$lonref,'_',crds$latref))
crds<-subset(crds,select=c('lon','lat','newarea5'))
names(crds)[1:2]<-c('clond','clatd')

#COMBINE ALL AND OUTPUT
return(crds)
}



###################################################################
#FUNCTION TO BIN DATA TO SPECIFIED RESOLUTION
binfunction<-function(d,stp){

if(dim(d)[1]>=1){
names(d)<-c('lond','latd')

#NEED TO PAD LON/LAT TO BE GLOBAL PRIOR TO PROJECTING
dum<-data.frame(lond=c(-180,-180,180,180),
                latd=c(-90,90,-90,90),
                dum=c(1,1,1,1))
d<-rbind.fill(d,dum)

#BINS LON INTO INTERVALS
brks<-seq(min(d$lon),max(d$lon),by=stp)#VALUES TO DEFINE GRID EDGES
lbls<-seq(min(brks)+(stp/2),max(brks)-(stp/2),by=stp)#LABELS FOR LON GRID
d$lonref<-round(as.numeric(as.character(cut(d$lon,breaks=brks,labels=lbls))),digits=4)

#BINS LAT INTO INTERVALS
brks<-seq(min(d$lat),max(d$lat),by=stp)#VALUES TO DEFINE GRID EDGES
lbls<-seq(min(brks)+(stp/2),max(brks)-(stp/2),by=stp)#LABELS FOR LON GRID
d$latref<-round(as.numeric(as.character(cut(d$lat,breaks=brks,labels=lbls))),digits=4)
d$newarea5<-gsub(' ','',paste(d$lonref,'_',d$latref))
d<-subset(d,is.na(dum)==TRUE,select=c('newarea5'))#REMOVE PADDED LAT/LON
return(d)
} else { return(subset(data.frame(lone=-9999,
                         late=-9999),lone>0))
}
}    
