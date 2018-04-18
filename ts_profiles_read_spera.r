library(plyr)
library(data.table)
setwd('N:/data/TS_profiles_maritimes/CTD_profiles_DFOhydrograph')
fls<-c('clim_8464_j.txt','clim_8462_j.txt','clim_8461_j.txt')

l<-list()
for(i in 1:length(fls)){
print(fls[i])
d<-read.table(paste(as.character(fls[i])),header=TRUE,sep=',')
d<-subset(d,SIGMAT<9000,select=c('STN_ID','DEPTH','PRESSURE','CRUISE_DATE','LONGITUDE','LATITUDE','TEMPERATURE','SALINITY','SIGMAT'))
l[[i]]<-d
gc()
}
dat1<-rbindlist(l)
names(dat1)<-tolower(names(dat1))

library(lubridate)
date<-strptime(dat1$cruise_date,"%d/%m/%Y")
dat1$year<-as.numeric(substr(date,1,4))
dat1$month<-as.numeric(substr(date,6,7))
dat1$dom<-as.numeric(substr(date,9,10))
dat1$doy<-yday(date)
dat1$cruise_date<-NULL

names(dat1)[1]<-'stationid'
names(dat1)[4]<-'lon'
names(dat1)[5]<-'lat'
dat1<-subset(dat1,lat>40 & lat< 46 & lon< -63 & lon > -68)

save(dat1,file='DFO_hydrographic_TS_profiles.RData')
setwd('N:/data/TS_profiles_maritimes/CTD_profiles_DFOhydrograph')
#load('DFO_hydrographic_TS_profiles.RData')



setwd('N:/data/TS_profiles_maritimes/RPettipas')
fls<-c('Data_2008-2015.csv','Data_2014.csv','Data_2015.csv','Data_2016.csv')

l<-list()
for(i in 1:length(fls)){
print(fls[i])
if(fls[i]=='Data_2014.csv'){
d<-read.csv(paste(as.character(fls[i])),header=TRUE)
}else {
d<-read.csv(paste(as.character(fls[i])),header=FALSE,col.names=c('MissionID','Latitude','Longitude','Year','Month','Day','Hour','Minute','Pressure','Temperature','Salinity','SigmaT','StationID'))
}
d<-subset(d,SigmaT> -1)
l[[i]]<-d
gc()
}
dat2<-rbindlist(l)
names(dat2)<-tolower(names(dat2))
dat2<-unique(dat2)
dat2<-subset(dat2,lat>40 & lat< 46 & lon< -63 & lon > -68)

library(lubridate)
date<-strptime(gsub(' ','',paste(dat2$day,'/',dat2$month,'/',dat2$year)),"%d/%m/%Y")
dat2$dom<-dat2$day
dat2$doy<-yday(date)

dat2<-subset(dat2,select=c('stationid','longitude','latitude','year','month','dom','doy','pressure','temperature','salinity','sigmat'))
names(dat2)<-c('stationid','lon','lat','year','month','dom','doy','pressure','temperature','salinity','sigmat')


x<-sin(dat2$lat/57.29578)^2
g<-9.780318*(1+((5.2788*10^-3)+(2.36*10^-5)*x)*x)+(1.092*10^-6)*dat2$pressure
dat2$depth<-(((((-1.82*10^-15) * dat2$pressure+ (2.279*10^-10)) * dat2$pressure- (2.2512*10^-5)) * dat2$pressure+ 9.72659) * dat2$pressure) / g

save(dat2,file='DFO_hydrographic_TS_profiles_RPettipas.RData')


#COMBINE ALL DATA
data<-rbind.fill(dat1,dat2)
setwd('N:/data/TS_profiles_maritimes/RPettipas')
save(data,file='DFO_hydrographic_TS_profiles_1914_2016_BoF.RData')

data$id<-gsub(' ','',paste(data$stationid,'_',data$year,'_',data$doy))

#ENSURE THAT PROFILES HAVE READINGS AT SURFACE (<5M) AND DEPTH (50M)
dum<-data.frame(id=sort(unique(data$id)),
            n=tapply(data$sigmat,data$id,function(x) length(unique(x))),
            mxd=tapply(data$depth,data$id,max),
            mnd=tapply(data$depth,data$id,min))
dum<-subset(dum,n>=10 & mxd>=25 & mnd<=5)
data<-subset(data,id %in% dum$id & salinity<9000)

library(akima)
d<-subset(data,id==sort(unique(data$id))[1])

stratfun<-function(d){
as<-aspline(d$depth,d$sigmat,xout=c(1,25,50))
    at<-aspline(d$depth,d$temperature,xout=c(1,25,50))
    asl<-aspline(d$depth,d$salinity,xout=c(1,25,50))
return(data.frame(id=unique(d$id),
                year=unique(d$year),
                month=unique(d$year),
                doy=unique(d$doy),
                dom=unique(d$dom),
                lon=unique(d$lon),
                lat=unique(d$lat),
                maxdepth=max(d$depth),
                mindepth=min(d$depth),
                s25=(as$y[2]-as$y[1])/(as$x[2]-as$x[1]),
                s50=(as$y[3]-as$y[1])/(as$x[3]-as$x[1]),
                temp1=at$y[1],
                temp25=at$y[2],
                temp50=at$y[3],
                sal1=asl$y[1],
                sal25=asl$y[2],
                sal50=asl$y[3]))
}
out<-dlply(data,.(id),.fun=stratfun,.progress='text')
out<-rbindlist(out)

out$s50<-ifelse(out$maxdepth<45,NA,out$s50)
out$temp50<-ifelse(out$maxdepth<45,NA,out$temp50)
out$sal50<-ifelse(out$maxdepth<45,NA,out$sal50)
