library(maps)
library(plyr)
library(lubridate)
library(rgdal)
library(rgeos)

#SPECIFY DIRECTORIES
datadir<-'N:/data/chl_phenology/data'
#figsdir<-'C:/Users/sailfish/Documents/aalldocuments/literature/postdoc_2013/chl_phenology/figures'
setwd(datadir)

#READ IN NAFO SHAPEFILE TO DEFINE SPATIAL BOUNDARIES
nafo<-readOGR('N:/data/dynamic_trophic_control/Assessments/data/fishing_shapefiles/nafo','Divisions')#works for Alaska - most others don't
nafo<-spTransform(nafo,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fourx<-subset(nafo,ZONE=='4X')
fourwx<-subset(nafo,ZONE%in% c('4X','4W'))

#######################################################################

olay<-function(dat,plg){
dat<-subset(dat,lat>41 & lat<46 & lon>-70 & lon< -58)
crds<-SpatialPoints(data.frame(lon=dat$lon,lat=dat$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gint<-over(crds,plg)
dat$dum<-gint[,1]
dat<-subset(dat,is.na(dum)==FALSE)
return(dat)
}


#N ATLANTIC CPR DATA
setwd(datadir)
data<-read.csv('CPR.2015.csv',header=TRUE)
names(data)<-tolower(names(data))

data<-subset(data,taxonname!='')
names(data)<-c('id','date','time','ampm','lat','lon','chl','taxonid','taxon','aphiaid','counts')
data<-subset(data,aphiaid!=0 & is.na(aphiaid)==FALSE)
#REPLACE ERRONEOUS APHIAIDS
data$aphiaid<-ifelse(data$aphiaid==156652,162936,data$aphiaid)
data$aphiaid<-ifelse(data$aphiaid==341581,109901,data$aphiaid)
data$aphiaid<-ifelse(data$aphiaid==367334,104473,data$aphiaid)

#ADDS TAXONOMIC HEIRARCHY - EXTRACTED PREVIOUSLY FROM WORMS
setwd('N:/data/iDiv_data/CPR')
taxolist<-read.csv('CPR_NWAtlantic_worms_taxolist.csv',header=TRUE)
taxolist<-subset(taxolist,select=c("aphia.id", "rank", "valid_name","kingdom","phylum","order","family","genus"))
names(taxolist)[3]<-'scientificname'
names(taxolist)[1]<-'aphiaid'

data<-merge(data,taxolist,by=c('aphiaid'),all.x=TRUE,all.y=FALSE)
data<-subset(data,rank!="")

#data$dir<-ifelse((regexpr("N",data$Latitude)>0)==TRUE, 'North',NA)
data$date<-strptime(data$date,"%d/%m/%Y")
data$year<-as.numeric(substr(data$date,1,4))
data$month<-as.numeric(substr(data$date,6,7))
data$dom<-as.numeric(substr(data$date,9,10))
data$djul<-as.numeric(julian(data$date,origin=as.POSIXct('1900-01-01',tz='GMT')))
data$quarter<-quarters(data$date,abbreviate=TRUE)
data$doy<-yday(data$date)
data<-subset(data,select=c("aphiaid","id","time", "ampm", "lat", "lon", "chl", "taxonid","taxon",  "counts", "rank", "scientificname", "kingdom","phylum", "order", "family", "genus","year", "month",'dom','djul','quarter','doy'))


a<-subset(data,phylum=='Cnidaria')
a<-subset(a,lon< -50 & lat>41)
dm<-data.frame(year=sort(unique(a$year)),
               n=tapply(a$counts,a$year,sum))
plot(dm$year,dm$n,pch=15)
plot(dm$year,log(dm$n+1),pch=15)


plot(data$lon,data$lat,pch='.')
map('world',add=TRUE,col='gray',fill=TRUE)

library(scales)
aa<-subset(a,counts>0)
plot(aa$lon,aa$lat,pch=16,col=alpha('darkred',.3),cex=2)
map('world',add=TRUE,col='gray',fill=TRUE)






####NOAA CPR
###############################################
#THIS GETS TAXONOMIC HEIRARCHY - ONLY NEED TO RUN ONCE
install.packages("XML", dependencies = TRUE)
download.file("http://www.omegahat.org/Prerelease/XMLSchema_0.8-0.tar.gz", "XMLSchema")
install.packages("XMLSchema", type="source", repos = NULL)
download.file("http://www.omegahat.org/Prerelease/SSOAP_0.91-0.tar.gz", "SSOAP")
install.packages("SSOAP", type="source", repos = NULL)
library(SSOAP)
w = processWSDL("http://www.marinespecies.org/aphia.php?p=soap&wsdl=1")
iface = genSOAPClientInterface(, w)
AphiaID = iface@functions$getAphiaID("Solea solea",1,('http://www.marinespecies.org/aphia.php?p=soap'))
print(AphiaID)
#should output '[0] 127160'


#GET APHIA ID FROM SPECIES NAME
AphiaMatch <- function(x) {
result<-NULL
for (i in 1:length(x)) {
AphiaRecord <- iface@functions$getAphiaID(x[i],1,('http://www.marinespecies.org/aphia.php?p=soap'))
result<-c(result, AphiaRecord)
}
return(result)
}

#FUNCTION TAKES WORMS APHIA ID AND RETREIVES TAXO INFO
getFullRecord <- function(x) {
result<-NULL
for (i in 1:length(x)) {
AphiaRecord <- iface@functions$getAphiaRecordByID(x[i],('http://www.marinespecies.org/aphia.php?p=soap'))
slotnames <- slotNames(AphiaRecord)
slotlist <- data.frame(rbind(1:length(slotnames)))
names(slotlist) <- slotnames
for(y in slotnames) {
#R cannot handle a slot name "class"
if (y == "CLASS") {slotlist[1,y] <- '(empty)'}
else {slotlist[1, y] <- slot(AphiaRecord,  y)}
}
result<-rbind(result, slotlist)
}
return(result)
}


#CPR DATA
setwd('N:/data/iDiv_data/CPR/CPR_GoM')
cprdat2<-read.csv('NOAACPRDATA.csv',header=TRUE,skip=6,col.names=c('cruise','sample','year','month','mday','hour','minute','lat','lon','chl','marmap.taxcode','abundance','sciname'))
cprdat2$lon<-cprdat2$lon*-1
cprdat2$date<-strptime(gsub(' ','',paste(cprdat2$mday,'/',cprdat2$month,'/',cprdat2$year)),"%d/%m/%Y")
cprdat2$day<-yday(cprdat2$date)
cprdat2$month<-as.numeric(substr(cprdat2$date,6,7))
cprdat2$dom<-as.numeric(substr(cprdat2$date,9,10))
cprdat2$djul<-as.numeric(julian(cprdat2$date,origin=as.POSIXct('1900-01-01',tz='GMT')))
cprdat2$quarter<-quarters(cprdat2$date,abbreviate=TRUE)
cprdat2$doy<-yday(cprdat2$date)

cprdat2<-subset(cprdat2,select=c('lon','lat','year','day','month','dom','djul','quarter','doy','marmap.taxcode','abundance','sciname'))
cprdat2$origsciname<-cprdat2$sciname

speclist<-unique(subset(cprdat2,select=c('origsciname','sciname')))
speclist$sciname<-gsub(" SPP.", "", speclist$sciname)
speclist$sciname<-gsub("Rhizosolenia hebetata semispina", "Rhizosolenia hebetata", speclist$sciname)
speclist$sciname<-gsub("Chaetoceros \\(Hyalochaetae) ", "Chaetoceros", speclist$sciname)
speclist$sciname<-gsub("Skeletonima costatum", "Skeletonema costatum", speclist$sciname)
speclist$sciname<-gsub("Rizosolenia", "Rhizosolenia", speclist$sciname)
speclist$sciname<-gsub("Rhizosolenia imricata shrobsolei", "Rhizosolenia imbricata", speclist$sciname)
speclist$sciname<-gsub("Proboscia  alata f. inermis", "Proboscia alata", speclist$sciname)
speclist$sciname<-gsub("Proboscia  alata", "Proboscia", speclist$sciname)
speclist$sciname<-gsub("RHIZOSOLENIA HEBETATA HIEMALUS", "RHIZOSOLENIA HEBETATA", speclist$sciname)
speclist$sciname<-gsub("Gonyaulacales \\(not Ceratium", "Gonyaulacales", speclist$sciname)
speclist$sciname<-gsub("RHIZOSOLENIA HEBETATA HIEMALUS", "RHIZOSOLENIA HEBETATA", speclist$sciname)
speclist$sciname<-gsub("Ceratium macroceros gallicum", "Tripos macroceros", speclist$sciname)
speclist$sciname<-gsub("SkeletonimaMelosira moniloformis", "Melosira moniliformis", speclist$sciname)
speclist$sciname<-gsub("CERATIUM VULTUR SUMATRANUM", "CERATIUM VULTUR", speclist$sciname)
speclist$sciname<-gsub("Ceratium evarcuatum", "Tripos euarcuatus", speclist$sciname)
speclist$sciname<-gsub("Melosira moniloformis", "Melosira moniliformis", speclist$sciname)

#RETREIVE APHIA ID FROM WORMS
speclist$aphiaid<-AphiaMatch(speclist$sciname)
speclist<-na.omit(speclist)

#ADD APHIA IDS THAT WORMS CAN'T FIND
speclist$aphiaid<-ifelse(speclist$sciname=='Ceratium tripos',109982,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Ceratium longipes',841259,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Ceratium macroceros',841260,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Rhizosolenia',149069,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Rhizosolenia alata f. indica',248181,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Navicula',149142,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Nitzschia',119270,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Ceratium massiliense',109968,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Hemiaulus',163248,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Ceratium trichoceros',109981,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Ceratium declinatum',841202,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Diplosalis',109516,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='Rhizosolenia acuminata',196805,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='CERATIUM GRACILE',841245,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$origsciname=='Ceratium evarcuatum',841210,speclist$aphiaid)

#RETRIEVE FULL TAXONOMIC HEIRARCHY FROM WORMS USING APHIA IDS
speclist<-subset(speclist,aphiaid>0)
aa<-getFullRecord(speclist$aphiaid)
taxolist<-cbind(speclist,aa)

taxolist<-subset(taxolist,select=c('origsciname','aphiaid', "rank", "valid_name","kingdom","phylum","order","family","genus"))
names(taxolist)[1]<-'sciname'
names(taxolist)[4]<-'scientificname'

cprdat3<-merge(cprdat2,taxolist,by=c('sciname'),all.x=TRUE,all.y=FALSE)
cprdat3<-subset(cprdat3,is.na(kingdom)==FALSE)

names(cprdat3)<-ifelse(names(cprdat3)%in% c('abundance'),'counts',names(cprdat3))
names(data)<-ifelse(names(data)%in% c('taxon'),'origsciname',names(data))
cprdat3$origin<-rep('NOAA',dim(cprdat3)[1])
data$origin<-rep('SAHFOS',dim(data)[1])

data2<-rbind.fill(data,cprdat3)


data2<-subset(data2,select=c("aphiaid","lat", "lon","counts", "rank", "scientificname", "kingdom","phylum", "order", "family", "genus","year", "month",'dom','djul','quarter','doy','origin','origsciname'))


cprdat<-olay(data2,fourwx)#DIVISION 4WX OVERLAY


library(maps)
map('world',xlim=c(-70,-60),ylim=c(40,48))
points(cprdat$lon,cprdat$lat,pch=16,cex=.5,col='red')
points(cprdat.ess$lon,cprdat.ess$lat,pch=16,cex=.5,col='red')

cprphyto<-subset(cprdat,kingdom %in% c('Chromista','Plantae'))
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
write.csv(cprphyto,'phyto_cpr_counts_spera.csv',row.names=FALSE)

cprzoop<-subset(cprdat,kingdom %in% c('Animalia'))
setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
write.csv(cprzoop,'zoop_cpr_counts_spera.csv',row.names=FALSE)


setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
phyt<-read.csv('phyto_cpr_counts_spera.csv',header=TRUE)
phyt$tax<-'phyto'
phyt<-subset(phyt,counts>0)

zoop<-read.csv('zoop_cpr_counts_spera.csv',header=TRUE)
zoop$tax<-'zoop'
zoop<-subset(zoop,counts>0)





#CALCULATE RICHNESS AND SHANNON H
#NEED TO REMOVE 0'S OR CAN'T CALCULATE SHANNON INDEX
richfun<-function(d){
d<-na.omit(d)
#RICHNESS
dout<-data.frame(richness=length(unique(d$scientificname)))

#RICHNESS WITH 1% CUTOFF LOCALLY
tcells<-sum(d$counts)
d$prp<-(d$counts/tcells)*100
d2<-subset(d,prp>=1)
dout$richness1<-length(unique(d2$scientificname))

#SHANNON I
f2<-function(d3){ return(data.frame(p=sum(d3$counts)/tcells)) }
sdat<-ddply(subset(d,select=c('scientificname','counts')),.(scientificname),.fun=f2)
sdat$lp<-log(sdat$p)
sdat$pp<--1*(sdat$p*sdat$lp)
sdat<-subset(sdat,is.na(pp)==FALSE)

dout$shannon<-sum(sdat$pp,na.rm=TRUE)
dout$pe<-sum(sdat$pp)/(log(length(unique(d$scientificname))))
dout$totalabundance<-sum(d$counts)
dout$nspecies<-length(unique(d$scientificname))

#COPEPODS
if(unique(d$tax)=='zoop'){
a<-subset(d,order %in%  c('Calanoida', 'Cyclopoida','Harpacticoida','Poecilostomatoida'))
dout$copepods<-sum(a$counts,na.rm=TRUE)
} else NULL
return(dout)
}
zrich<-ddply(zoop,.(year,month,dom,doy,djul,origin,lon,lat),.fun=richfun,.progress='text')
names(zrich)[9:15]<-gsub(' ','',paste(names(zrich)[9:15],'.z'))

prich<-ddply(phyt,.(year,month,dom,doy,djul,origin,lon,lat),.fun=richfun,.progress='text')
names(prich)[9:14]<-gsub(' ','',paste(names(prich)[9:14],'.p'))

cprrich<-merge(zrich,prich,by=c('year','month','dom','doy','djul','origin','lon','lat'),all=TRUE)

setwd('N:/cluster_2017/scratch/spera/data/stagingdat')
write.csv(cprrich,'plank_cpr_richness_spera.csv',row.names=FALSE)







#####
