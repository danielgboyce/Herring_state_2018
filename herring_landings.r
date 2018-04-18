library(tidyr)
library(rgdal)
library(ggplot2)
library(gganimate)
library(lubridate)
library(plyr)
datadir<-'N:/data/Spera/herring/landings'
figsdir<-"C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures"
figsdir<-"C://Users//copepod//Documents//aalldocuments//literature//research//active//SPERA//Figures"

#COMBINES HERRING LANDINGS BY COUNTY FROM
#1) HISTORICAL EXTRACTIONS (1900-1963)
#2) COMMLAND DB (1967-2001)
#3) MARFIS DB (2000-ONWARD)

NEED TO FIGURE OUT
- DO MARFIS DATA CONTAIN ALL COUNTIES (NS, AND NB)
- WHAT ARE UNITS IN COMLAND
- DO I HAVE TOO MANY OBS IN MARFIS? DB PERSON TOLD ME TO BE CAREFUL ABOUT 'MON_DOC_ID'


#######################################################
#READ IN AND FORMAT MARFIS LANDINGS DATA (2000 ONWARD)
#UNITS ARE IN KGS
setwd(datadir)
pro.spec<-read.csv('pro_spec_info.csv',header=TRUE)
names(pro.spec)<-tolower(names(pro.spec))
dt<-strptime(pro.spec$landed_date,"%y-%m-%d")
pro.spec$year<-as.numeric(substr(dt,1,4))
pro.spec<-subset(pro.spec,species_code==200)

lond<-as.numeric(substr(pro.spec$longitude,1,2))
lonm<-round(as.numeric(substr(pro.spec$longitude,3,4))/60,digits=2)*100
lons<-as.numeric(substr(pro.spec$longitude,6,6))
lons<-ifelse(is.na(lons)==TRUE,0,lons)
pro.spec$lon<-as.numeric(gsub(' ','',paste(lond,'.',lonm,lons)))*-1


latd<-as.numeric(substr(pro.spec$latitude,1,2))
latm<-round(as.numeric(substr(pro.spec$latitude,3,4))/60,digits=2)*100
lats<-as.numeric(substr(pro.spec$latitude,6,6))
lats<-ifelse(is.na(lats)==TRUE,0,lats)
pro.spec$lat<-as.numeric(gsub(' ','',paste(latd,'.',latm,lats)))

setwd(figsdir)
pdf('marfis_herring_landings.pdf')
map('world',xlim=c(-70,-45),ylim=c(40,50),col='black')
points(pro.spec$lon,pro.spec$lat,pch=16,col=alpha('darkred',.1),cex=.5)
map.axes()
map('world',xlim=c(-70,-45),ylim=c(40,50),col='black',add=TRUE)
dev.off()
library(maps)

#TABLE TO LINK COMMUNITIES TO COUNTIES
comm<-read.csv('communities.csv',header=TRUE)
names(comm)<-tolower(names(comm))

#TABLE CONTAINS COUNTIES
cnty<-read.csv('districts_counties_provinces.csv',header=TRUE)
names(cnty)<-tolower(names(cnty))
names(cnty)[1]<-'district_id'

ps<-join(pro.spec,comm,by=c('community_code'),type='left')
ps<-join(ps,cnty,by=c('district_id'),type='left')
ps<-subset(ps,is.na(county_name)==FALSE)


f<-function(d){
    return(data.frame(lnd=sum(d$rpt_weight_kgs,na.rm=TRUE)))
}
marf<-ddply(ps,.(year,county_name),.fun=f,.progress='text')


##########################################################
#PLOTS SPATIAL DISTRIBUTION OF LANDINGS BY COUNTY
f<-function(d){
map('world',xlim=c(-70,-45),ylim=c(40,50),col='black')
points(d$lon,d$lat,pch=16,col=alpha('green',1),cex=.75)
map.axes()
map('world',xlim=c(-70,-45),ylim=c(40,50),col='black',add=TRUE)
legend('bottomright',legend=unique(d$county_name),bty='n')
}
setwd(figsdir)
pdf('marfis_herring_landings.pdf',height=12,width=10)
par(mfrow=c(4,3),mar=c(1,1,1,1))
z<-dlply(ps,.(county_name),.fun=f,.progress='text')
dev.off()
library(maps)


############################################################
#READ IN AND FORMAT COMLAND LANDINGS DATA (2000 ONWARD)
cnty<-read.csv('districts_counties_provinces.csv',header=TRUE)
names(cnty)<-tolower(names(cnty))
names(cnty)[1]<-'dist.id'
load('herland.RData')
names(herland)<-tolower(names(herland))
dat<-herland
names(dat)[11]<-'dist.id'
dat<-subset(dat,select=c('dist.id','land_date','live_wt','trip_num','unit_code'))
dat<-subset(dat,unit_code!='K' | is.na(unit_code)==TRUE)
dat$year<-as.numeric(substr(dat$land_date,1,4))
dat$month<-as.numeric(substr(dat$land_date,6,7))
dat$mday<-as.numeric(substr(dat$land_date,9,10))
dat<-merge(dat,cnty,by=c('dist.id'),all.x=TRUE,all.y=FALSE)
dat$county_name<-tolower(dat$county_name)
dat<-subset(dat,is.na(county_name)==FALSE&county_name!='unknown ns')

f<-function(d){
    d2<-(data.frame(year=sort(unique(d$year)),
                    lnd=tapply(d$live_wt,d$year,function(x) sum(x,na.rm=TRUE))))
    db<-data.frame(year=seq(1967,2001,1))
    d2<-merge(d2,db,by=c('year'),all=TRUE)
return(d2)
    }
cland<-ddply(dat,.(county_name),.fun=f,.progress='text')
cland$lnd<-ifelse(is.na(cland$lnd)==TRUE,0,cland$lnd)




##############################################################

##############################################################
#HISTORICAL EXTRACTED DATA
#UNITS IN 000'S LBS
setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/data/herring')
hdat<-read.csv('herring_landings_historical_bycounty.csv',header=TRUE)
names(hdat)[28]<-'unit'
hdats<-hdat %>% gather(county_name,y,charlotte:kings)
hdats$y<-ifelse(hdats$year<1947,hdats$y/10,hdats$y)#CONVERT HUNDREDWEIGHTS TO 1000'S OF LBS
hdats$y<-(hdats$y*1000)*0.45359237#CONVERT TO KGS

f<-function(d){
    return(data.frame(lnd=sum(d$y,na.rm=TRUE)))
}
histl<-ddply(hdats,.(year,county_name),.fun=f,.progress='text')


data<-rbind.fill(cland,histl,marf)
data$county_name<-tolower(data$county_name)
data<-subset(data,!(county_name %in% c('unknown','unknown ns')))
data<-subset(data,!(county_name %in% c('inverness','victoria','antigonish','pictou','cape breton','kent','gloucester','restigouche','northumberland')))




############################################################
# DIVERSITY, EVENNESS, RICHNESS OFLANDINGS ACROSS COUNTIES
dfun<-function(d){

dout<-data.frame(richness=length(unique(d$county_name)))
#RICHNESS WITH 1% CUTOFF LOCALLY
tland<-sum(d$lnd)
d$prp<-(d$lnd/tland)*100
d2<-subset(d,prp>=1)

#SHANNON I
f2<-function(d3){ return(data.frame(p=sum(d3$lnd)/tland)) }
sdat<-ddply(subset(d,select=c('county_name','lnd')),.(county_name),.fun=f2)
sdat$lp<-log(sdat$p)
sdat$pp<--1*(sdat$p*sdat$lp)
sdat<-subset(sdat,is.na(pp)==FALSE)

dout$shannon<-sum(sdat$pp,na.rm=TRUE)
dout$pe<-sum(sdat$pp)/(log(length(unique(d$county_name))))
dout$totalabundance<-sum(d$lnd)
dout$nspecies<-length(unique(d$county_name))
return(dout)
}
ldiv<-ddply(subset(data,lnd>0),.(year),.fun=dfun,.progress='text')
ldiv<-subset(ldiv,is.na(pe)==FALSE)


#PROPORTION OF ALL CATCH LANDED IN TOP COUNTY
dfun2<-function(d){
dout<-data.frame(richness=length(unique(d$county_name)))
tland<-sum(d$lnd)

#SHANNON I
f2<-function(d3){ return(data.frame(p=round(sum(d3$lnd)/tland,digits=2))) }
sdat<-ddply(subset(d,select=c('county_name','lnd')),.(county_name),.fun=f2)
sdat<-subset(sdat,p %in% sort(unique(sdat$p),decreasing=TRUE)[1])
sdat$p<-sdat$p*100
return(sdat)
}
ldiv2<-ddply(subset(data,lnd>0 & year<2018),.(year),.fun=dfun2,.progress='text')
dm<-data.frame(county_name=sort(unique(ldiv2$county_name)),
   cl=c('firebrick3','gold3','yellow','pink','dodgerblue3','lightblue','forestgreen'))
ldiv2<-merge(ldiv2,dm,by=c('county_name'),all.x=TRUE,all.y=FALSE)


#PROPORTION OF ALL CATCH LANDED IN TOP COUNTY
dfun3<-function(d){
dout<-data.frame(richness=length(unique(d$county_name)))
tland<-sum(d$lnd)

#SHANNON I
f2<-function(d3){ return(data.frame(p=round(sum(d3$lnd)/tland,digits=2))) }
sdat<-ddply(subset(d,select=c('county_name','lnd')),.(county_name),.fun=f2)
sdat$p<-sdat$p*100
return(sdat)
}
ldiv3<-ddply(subset(data,lnd>0 & year<2018),.(year),.fun=dfun3,.progress='text')

f<-function(d){
    plot(d$year,d$p,pch=15,las=1)
    legend('top',unique(d$county_name),bty='n')
}
par(mfrow=c(5,5),mar=c(1,1,1,1))
z<-dlply(ldiv3,.(county_name),.fun=f)
median(c(87,56,27,35,44,34,26,82,79,37,68,88,61,63,51,42))
mean(c(87,56,27,35,44,34,26,82,79,37,68,88,61,63,51,42))

f<-function(d){
    return(data.frame(lnd=sum(d$lnd)))
    }
dm<-ddply(data,.(county_name,year),.fun=f,.progress='text')


setwd(figsdir)
pdf('landings_halifaxcounty.pdf',height=7,width=8)
par(mfrow=c(2,2))
a<-subset(ldiv3,county_name=='halifax')
a$cl<-ifelse(a$year<1967,'firebrick3','dodgerblue3')
a$cl<-ifelse(a$year>2000,'forestgreen',a$cl)
plot(a$year,a$p,col=alpha(as.character(a$cl),.5),pch=16,las=1,xlab='Year',ylab='Proportion of all landings',main='Halifax')
points(a$year,a$p,pch=1,lwd=.001)
legend('top',legend=c('Bureau of statistics','ComLand','Marfis'),col=c('firebrick3','dodgerblue3','forestgreen'),bty='n',pch=15)

b<-subset(ldiv3,county_name=='queens')
b$cl<-ifelse(b$year<1967,'firebrick3','dodgerblue3')
b$cl<-ifelse(b$year>2000,'forestgreen',b$cl)
plot(b$year,b$p,col=alpha(as.character(b$cl),.5),pch=16,las=1,xlab='Year',ylab='Proportion of all landings',main='Queens')
points(b$year,b$p,pch=1,lwd=.001)
legend('top',legend=c('Bureau of statistics','ComLand','Marfis'),col=c('firebrick3','dodgerblue3','forestgreen'),bty='n',pch=15)
dev.off()


library(akima)
library(mgcv)
library(RColorBrewer)

pfun<-function(d,ylm,ylb){
names(d)[1:2]<-c('year','y')
#ylm<-c(min(d$y),max(d$y))
pdat<-data.frame(year=seq(min(d$year),max(d$year),length.out=1000))
mod<-gam(y~s(year,k=20),data=d)
pdat$y<-aspline(d$year,d$y,xout=pdat$year)$y
pdat$y2<-predict(mod,newdata=pdat,se.fit=FALSE)
#plot(pdat$year,pdat$y,type='l',las=1,xlab='Year')
plot(0,0,las=1,xlab='Year',ylim=ylm,xaxt='n',xlim=c(1915,2020),ylab=ylb,col='white')
axis(1,seq(1910,2020,5),labels=FALSE)
axis(1,seq(1910,2020,10),labels=TRUE)
points(d$year,d$y,pch=16,col=alpha('gold3',.9),cex=1.5)
points(d$year,d$y,pch=1,col='gray50',cex=1.5,lwd=.01)
lines(pdat$year,pdat$y2,col=alpha('gray30',.75),lwd=2,lty=2)
}
setwd(figsdir)
pdf('landings_bycounty_tsplots.pdf',height=8,width=11)
par(mfrow=c(2,2),mar=c(4,4,1,1))
pfun(subset(ldiv,select=c('year','richness')),c(7,19),'Number of counties reporting landings')
pfun(subset(ldiv,select=c('year','shannon')),c(1,2.5),'Diversity of landings across counties')
pfun(subset(ldiv2,year<2018,select=c('year','p')),c(10,70),'Proportion of all landings accounted for by top county')

plot(ldiv2$year,ldiv2$p,col=alpha(as.character(ldiv2$cl),.75),pch=16,cex=1.5,las=1,ylim=c(10,70),xlab='Year',ylab='Proportion of all landings accounted for by top county',xaxt='n')
axis(1,seq(1910,2020,5),labels=FALSE)
axis(1,seq(1910,2020,10),labels=TRUE)
points(ldiv2$year,ldiv2$p,pch=1,lwd=.01,cex=1.5)
legend('topleft',legend=unique(dm$county_name),col=unique(as.character(dm$cl)),pch=15,bty='n')
dev.off()




plot(ldiv2$year,ldiv2$p)


dat2$id<-dat2$county_name


setwd('C:/Users/sailfish/Downloads')
shp<-readOGR('C:/Users/sailfish/Downloads',layer='gcd_000b11a_e')#works for Alaska - most
a<-subset(shp,PRNAME %in% c('Nova Scotia / Nouvelle-Écosse','New Brunswick / Nouveau-Brunswick'))
a$CDNAME<-tolower(a$CDNAME)

a<-subset(a,!(CDNAME %in% c('restigouche','madawaska','gloucester','northumberland','victoria','carleton','york','sunbury')))
a2<-subset(a,!(CDNAME %in% c('restigouche','madawaska','gloucester','northumberland','victoria','carleton','york','sunbury','kent','madawaska','carleton','york','sunbury','pictou','inverness','antigonish')))
a2<-subset(a2,!(PRNAME=='New Brunswick / Nouveau-Brunswick' & CDNAME=='queens'))
a2<-subset(a2,!(PRNAME=='New Brunswick / Nouveau-Brunswick' & CDNAME=='kings'))

setwd(figsdir)
pdf('counties.pdf',height=10,width=10)
plot(a,lwd=.001)
plot(a2,col='firebrick3',add=TRUE,border='black',lwd=.001)
cdat <- data.frame(getSpPPolygonsLabptSlots(a2))
cdat$cnty<-a2$CDNAME
text(cdat$X1,cdat$X2,labels=cdat$cnty,col='white',cex=.75)
dev.off()


cdat <- data.frame(getSpPPolygonsLabptSlots(a))
cdat$cnty<-a$CDNAME
text(cdat$X1,cdat$X2,labels=cdat$cnty,col='black',cex=.75)

#am<-fortify(a,county='id')
a$id<-a$CDNAME
am<-fortify(a,region='id')
am<-subset(am,id%in% dat2$id)
mydat<-merge(am,dat2,by=c('id'))

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_blank(),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position="right",
                         plot.title = element_text(size=16)))

mydat$lnd<-ifelse(mydat$lnd==0,NA,mydat$lnd)
pdf('prac.pdf')
ggplot(mydat, aes(x = long, y = lat, group=group, fill=lnd)) +
  theme_opts +
  coord_equal() +
  geom_polygon(data=a,aes(long,lat,group=group),fill=NA,colour='black') +
  #  geom_polygon(color = 'grey') +
  labs(title = "Herring landings between 1967 and 2015",
       subtitle = "Static Map",
       caption = 'Data source: Statistics Canada') +
  scale_fill_distiller(palette='Spectral')
dev.off()

mydat$lnd2<-log10(mydat$lnd+1)
hist(mydat$lnd)

p<-ggplot(mydat, aes(x = long, y = lat, group = group, fill=lnd2,frame=year)) +
  theme_opts +
  coord_equal() +
  geom_polygon(color = 'grey',size=.0001) +
#  geom_polygon(aes(long,lat,group=group),fill=NA,colour='black',data=a) +
  labs(title = "Herring landings between 1967 and 2001",
       subtitle = "Dynamic Map",
       caption = 'Data source: Statistics Canada') +
  scale_fill_distiller(palette='Spectral')

  gganimate(p,'output.gif')
  gganimate(p,'output.log.gif')
  gganimate(p,'output.mp4')
#devtools::install_github('dgrtwo/gganimate')
#devtools::install_github('yihui/animation')
library(animation)

unique(a$CDNAME)
plot(a,col='white',border=NA)
map('world',add=TRUE,fill=TRUE,col='gray',border=NA)
plot(a,add=TRUE,lwd=.1,col='cornflowerblue',border='darkblue')
plot(a)
text(a,labels=as.character(a$CDNAME))

library(googleVis)
unique(dat$unit_code)



