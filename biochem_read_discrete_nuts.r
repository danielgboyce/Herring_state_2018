library(plyr)
library(fossil)
library(maps)
library(mapdata)
library(lubridate)
library(scales)
#NOTES:
#TIMES: STRANGE FORMAT: SOME ARE 0,1,2,... OTHERS 24 HOUR FORMAT

figsdir<-'N:/data/plankton_BIOCHEM/figures'
datadir<-'N:/data/plankton_BIOCHEM/data'

#DATA FORMATTING: ORGANIZED IN TO MISSIONS (TRIPS OR STUDIES) AND EVENTS
#'PARAMETER.NAME' GIVES WHAT WAS MEASURED, 'HEADER.START.DEPTH' GIVES DEPTH AT WHICH IT WAS MEASURED


#READS IN COORDINATES FOR ALL SECTION STATIONS
setwd(datadir)
azmp<-read.csv('azmp_sections_stations.csv',header=FALSE,skip=3,col.names=c('id','program',"section.code", "section.name","section.name.fra", "station.code","station.name","station.name.fra", "lat","lon", "depth"))
azmp$lon<-azmp$lon*-1
azmp<-subset(azmp,id<=161)

datad<-read.csv('dan_discrete.csv',header=TRUE)#DISCRETE CHL
#DATA FORMATTING: ORGANIZED IN TO MISSIONS (TRIPS OR STUDIES) AND EVENTS
#'PARAMETER.NAME' GIVES WHAT WAS MEASURED, 'HEADER.START.DEPTH' GIVES DEPTH AT WHICH IT WAS MEASURED
names(datad)<-tolower(names(datad))
unique(datad$parameter_name)
gc()

datad<-subset(datad,header_start_depth<=20)

#REMOVE INSTANCES WHERE START AND END DEPTH DIFFER
datad$deltad<-datad$header_start_depth-datad$header_end_depth
#datad<-subset(datad,abs(deltad)==0)

#FUNCTION FORMATS HEADER START DATE TO GET YEARS, MONTHS, DAYS, ETC..
datefun<-function(d){
names(d)<-tolower(gsub('_','.',names(d)))
d$mission.sdate<-strptime(d$mission.start,"%d-%b-%Y")
d$mission.edate<-strptime(d$mission.end,"%d-%b-%Y")
d$event.sdate<-strptime(d$event.start,"%d-%b-%Y")
d$mission.sjul<-as.numeric(julian(d$mission.sdate,origin=as.POSIXct('1900-01-01',tz='GMT')))
d$mission.ejul<-as.numeric(julian(d$mission.edate,origin=as.POSIXct('1900-01-01',tz='GMT')))
d$event.sjul<-as.numeric(julian(d$event.sdate,origin=as.POSIXct('1900-01-01',tz='GMT')))

d$date<-strptime(d$header.start,"%d-%b-%Y")
d$year<-as.numeric(substr(d$date,1,4))
d$month<-as.numeric(substr(d$date,6,7))
d$dom<-as.numeric(substr(d$date,9,10))
d$djul<-as.numeric(julian(d$date,origin=as.POSIXct('1900-01-01',tz='GMT')))
d$quarter<-quarters(d$date,abbreviate=TRUE)
d$doy<-yday(d$date)

d$lon<-d$header.start.lon
d$lat<-d$header.start.lat

d$oc<-ifelse(d$lat <= 65 & d$lat > 30 & d$lon > -110 & d$lon <= 60  | d$lat <= 30 & d$lat > 20 &
        d$lon > -100 & d$lon <= 30 | d$lat <= 20 & d$lat > 16 & d$lon > -98 & d$lon <= 20
        | d$lat <= 10 & d$lat > -60 & d$lon > -70 & d$lon <= 20 | d$lat <= 20 & d$lat > 9 & d$lon > -83 & d$lon <= 20,'atlantic','pacific')
d$oc<-ifelse(d$lat<50 & d$lat>40 & d$lon>27 & d$lon<45,'black',d$oc)
d$oc<-ifelse(d$lat<50 & d$lat>30 & d$lon>42 & d$lon<65,'caspian', d$oc)
d$oc<-ifelse(d$lat<48 & d$lat>30 & d$lon> 0 & d$lon<25 | d$lat<40 & d$lat>30 & d$lon>= 25 & d$lon<40 |  d$lat<40 & d$lat>30 & d$lon> -5 & d$lon<=0,'med',d$oc)
d$oc<-ifelse(d$lat <= 30 & d$lat > -60 & d$lon > 20 & d$lon <= 100 | d$lat <=0 & d$lat > -60 & d$lon >100 & d$lon <= 110 | d$lat <=-10 & d$lat > -60 & d$lon > 110 & d$lon <= 140 | d$lat <= - 30 & d$lat > -60 & d$lon > 140 & d$lon <= 150,'indian',d$oc)
d$oc<-ifelse(d$lat >= 65,'arctic',d$oc)
d$oc<-ifelse(d$lat <= -60,'southern',d$oc)
d$oc<-ifelse(d$lon>=105 & d$lon<120 & d$lat<5 & d$lat>-7,'pacific',d$oc)
d$oc<-ifelse(d$lon>=135 & d$lon<160 & d$lat<5 & d$lat>-20,'pacific',d$oc)
d$oc<-ifelse(d$lon>=100 & d$lon<=133 & d$lat< -9 & d$lat>-20,'indian',d$oc)
d$oc<-ifelse(d$lon>=109.9 & d$lon<=110.5 & d$lat< -7 & d$lat>-9,'indian',d$oc)
d$oc<-ifelse(d$lon>=99 & d$lon<=100 & d$lat< 12 & d$lat>10,'pacific',d$oc)
d$oc<-ifelse(d$lon>=20 & d$lon<=30 & d$lat< 41 & d$lat>38,'med',d$oc)

#adds grid cells
d$lon<-ifelse(d$lon== -180, -179.99,d$lon)
d$lon<-ifelse(d$lon== 180, 179.99,d$lon)
latref10 = floor((d$lat - (90 + 10))/10*-1); lonref10 = floor((d$lon - (180 + 10))/10*-1); d$newarea10 = latref10 + (lonref10/100)#10 degree
latref5 = floor((d$lat - (90 + 5))/5*-1); lonref5 = floor((d$lon - (180 + 5))/5*-1); d$newarea5 = latref5 + (lonref5/100)#10 degree
latref = floor((d$lat - (90 + 1))/1*-1);lonref = floor((d$lon - (180 + 1))/1*-1);d$newarea1<- latref + (lonref/1000)
latref = floor((d$lat - (90 + 2))/2*-1);lonref = floor((d$lon - (180 + 2))/2*-1);d$newarea2<- latref + (lonref/1000)

return(d)
}
datad<-datefun(datad)
datad<-subset(datad,lon> -68 & lon< -64 & lat<45 & lat>42 & oc %in% c('atlantic'))
#datad$date<-as.Date(datad$date)




library(mapdata)
map('worldHires',fill=TRUE,col='gray',xlim=c(-70,-60),ylim=c(42,45))
map.axes()
points(datad$lon,datad$lat,pch=16,col='red')




names(datad)<-ifelse(names(datad) %in% c('parameter.name'),'var',names(datad))
names(datad)<-ifelse(names(datad) %in% c('data.value'),'value',names(datad))
names(datad)<-ifelse(names(datad) %in% c('header_start_depth'),'sdepth',names(datad))
names(datad)<-ifelse(names(datad) %in% c('header_end_depth'),'edepth',names(datad))

datad<-subset(datad,select=c('mission.seq','event.seq','data.type.desc','var','value','data.quality','year','month','dom','djul','doy','lon','lat'))
datad<-subset(datad,var %in% c('Phosphate','Silicate','Nitrate'))

datad<-subset(datad,!(data.quality %in% c('QC has been performed; element appears to be doubtful.','QC has been performed; element appears to be erroneous.','QC has been performed; element appears to be inconsistent with other elements.')))


nit<-subset(datad,var=='Nitrate')
sil<-subset(datad,var=='Silicate')
phos<-subset(datad,var=='Phosphate')

f<-function(d){return(value=mean(d$value))}
nit<-ddply(nit,.(year,month,dom,djul,doy,lon,lat),.fun=f)
sil<-ddply(sil,.(year,month,dom,djul,doy,lon,lat),.fun=f)
phos<-ddply(phos,.(year,month,dom,djul,doy,lon,lat),.fun=f)
names(nit)[8]<-'value'
names(sil)[8]<-'value'
names(phos)[8]<-'value'

sil<-subset(sil,value<15)
phos<-subset(phos,value<2.5)

setwd('N:/data/Spera/final')
write.csv(nit,'biochem_nitrate_spera.csv',row.names=FALSE)
write.csv(sil,'biochem_silicate_spera.csv',row.names=FALSE)
write.csv(phos,'biochem_phosphate_spera.csv',row.names=FALSE)

par(mfrow=c(2,2))
hist(nit$value,breaks=100)
hist(sil$value,breaks=100)
hist(phos$value,breaks=100)
