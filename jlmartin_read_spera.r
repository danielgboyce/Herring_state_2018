library(lattice)
library(mgcv)
library(plyr)
library(maps)
datadir<-'N:/data/iDiv_data/JMARTIN/data'
setwd(datadir)

data<-read.csv('jlmartin_mergedplanktondata.csv',header=TRUE)
data$lcount<-log10(data$count+1)
data$sid<-gsub(' ','',paste(data$stn,'_',data$date,'_',data$depth))
stns<-c('wolves','brandy_cove','deadmans_harbour','lime_kiln_bay')#STATIONS THAT SPAN 1987-2011
data<-subset(data,stn %in% stns)


jlphyto<-subset(data,kingdom %in% c('Plantae','Chromista'))
jlzoop<-subset(data,kingdom %in% c('Animalia','Protozoa'))


#CALCULATE RICHNESS, EVENNESS
#NEED TO REMOVE 0'S OR CAN'T CALCULATE SHANNON INDEX
richfun<-function(d){
dout<-unique(subset(d,select=c('stn','lon','lat','station','depth','year','month','day','djul','quarter','doy','bathy')))
dout$scount<-sum(d$count)
dout$lscount<-log10(dout$scount+1)

d<-subset(d,count>0)
#RICHNESS
dout$richness<-length(unique(d$scientificname))
#RICHNESS WITH 1% CUTOFF LOCALLY
tcells<-sum(d$count)
d$prp<-(d$count/tcells)*100
d2<-subset(d,prp>=1)
dout$richness1<-length(unique(d2$scientificname))

if(dim(d)[1]>1){
#SHANNON I
f2<-function(d3){ return(data.frame(p=sum(d3$count)/tcells)) }
sdat<-ddply(subset(d,select=c('scientificname','count')),.(scientificname),.fun=f2)
sdat$lp<-log10(sdat$p)
sdat$pp<--1*(sdat$p*sdat$lp)
sdat<-subset(sdat,is.na(pp)==FALSE)
dout$shannon<-sum(sdat$pp,na.rm=TRUE)
dout$pe<-sum(sdat$pp)/(log(length(unique(d$scientificname))))
dout$nspecies<-length(unique(d$scientificname))
} else {
dout$shannon<-0
dout$pe<-NA
dout$nspecies<-0
}
return(dout)
}
jlphyto2<-ddply(jlphyto,.(sid),.fun=richfun,.progress='text')
jlzoop2<-ddply(jlzoop,.(sid),.fun=richfun,.progress='text')


#EXAMIEN CHANGES IN DEPTH OF OCCUPANCY - NONE
depstn<-c('deadmans_harbour','lime_kiln_bay','wolves')

f<-function(d){
print(unique(d$stn))
print(summary(d$depth))
names(d)<-ifelse(names(d) %in% c('scount'),'y',names(d))
d$depth<-d$depth+1
 mod<-gam(depth~s(day,k=5) +s(year),weights=d$y,data=d,family=Gamma('log'))
pdat<-data.frame(year=seq(min(d$year),max(d$year),length.out=100),
                 day=250)
    p<-predict(mod,newdata=pdat,type='link')
    pdat$p<-p
return(pdat)
}
ddat<-ddply(subset(jlzoop2,stn %in% depstn),.(stn),.fun=f)


#ESTIMATE ANNUAL VALUES IN Z/P COMMUNITY COMPOSITION AVERAGED ACROSS DIFFERENT STATIONS USING MIXED MODEL WITH RANDOM INTERCEPT
f<-function(d,nm,nm2){
#print(unique(d$stn))
#names(d)<-ifelse(names(d) %in% c('scount'),'y',names(d))
d$y<-log10(d$scount+1)

mod<-gamm(y~s(day,k=5) +as.factor(year) + s(depth,k=4),data=d,random=list(stn=~1))
mods<-gamm(shannon~s(day,k=5) +as.factor(year) + s(depth,k=4),data=d,random=list(stn=~1))
modr<-gamm(richness~s(day,k=5) +as.factor(year) + s(depth,k=4),data=d,random=list(stn=~1))
pdat<-data.frame(year=sort(unique(d$year)),
                 day=250,
                 depth=1)
    p<-predict(mod$gam,newdata=pdat,se.fit=TRUE)
    ps<-predict(mods$gam,newdata=pdat,se.fit=TRUE)
    pr<-predict(modr$gam,newdata=pdat,se.fit=TRUE)
    pdat$mn<-p$fit
    pdat$div<-ps$fit
    pdat$rich<-pr$fit
    pdat$mn.se<-p$se.fit
    pdat$div.se<-ps$se.fit
    pdat$rich.se<-pr$se.fit

#FOR EVENNESS: LOTS OF MISSINGS SO NEED ALTERNATE APPROACH
d2<-subset(d,is.na(pe)==FALSE)
pdat2<-data.frame(year=sort(unique(d2$year)),
                  day=250,
                  depth=1)
if(dim(d2)[1]>20){
modp<-gamm(pe~s(day,k=5) +as.factor(year) + s(depth,k=4),data=d2,random=list(stn=~1))
    pp<-predict(modp$gam,newdata=pdat2,se.fit=TRUE)
    pdat2$pe<-pp$fit
    pdat2$pe.se<-pp$se.fit
} else {
    pdat2$pe<-NA
    pdat2$pe.se<-NA
}

pdat<-merge(pdat,pdat2,by=c('year','day','depth'),all=TRUE)
names(pdat)[4:11]<-gsub(' ','',paste(nm,'.',names(pdat)[4:11],'.',nm2))
return(pdat)

}
tdatz<-f(jlzoop2,'zp','sabs')
tdatp<-f(jlphyto2,'ph','sabs')
tdat<-merge(tdatz,tdatp,by=c('year','day','depth'),all=TRUE)
tdat<-tdat[,-2]
tdat<-tdat[,-2]





###############################################################
#ESTIMATES TRENDS FROM STATIONS INDIVIDUALLY TO VERIFY THAT MIXED MODEL IS APPROPRIATE
f<-function(d){
print(unique(d$stn))
#names(d)<-ifelse(names(d) %in% c('scount'),'y',names(d))
d$y<-log10(d$scount+1)

#d$y<-(d$y-mean(d$y))/sd(d$y)
if(length(unique(d$depth))>3){
mod<-gam(y~s(day,k=5) +as.factor(year) +s(depth,k=4),data=d)
mods<-gam(shannon~s(day,k=5) +as.factor(year) +s(depth,k=4),data=d)
modr<-gam(richness~s(day,k=5) +as.factor(year),data=d)
} else {
mod<-gam(y~s(day,k=5) +as.factor(year),data=d)
modr<-gam(richness~s(day,k=5) +as.factor(year),data=d)
mods<-gam(shannon~s(day,k=5) +as.factor(year),data=d)
}
pdat<-data.frame(year=sort(unique(d$year)),
                     day=250,
                     depth=1)
    p<-predict(mod,newdata=pdat,se.fit=TRUE)
    ps<-predict(mods,newdata=pdat,se.fit=TRUE)
    pr<-predict(modr,newdata=pdat,se.fit=TRUE)
    pdat$scount<-p$fit
    pdat$shannon<-ps$fit
    pdat$rich<-pr$fit
    pdat$scount.se<-p$se.fit
    pdat$shannon.se<-ps$se.fit
    pdat$rich.se<-pr$se.fit

#FOR EVENNESS: LOTS OF MISSINGS SO NEED ALTERNATE APPROACH
d2<-subset(d,is.na(pe)==FALSE)
pdat2<-data.frame(year=sort(unique(d2$year)),
                     day=250,
                     depth=1)
if(length(unique(d2$depth))>3 & dim(d2)[1]>20){
modp<-gam(pe~s(day,k=5) +as.factor(year) +s(depth,k=4),data=d2)
    pp<-predict(modp,newdata=pdat2,se.fit=TRUE)
    pdat2$pe<-pp$fit
    pdat2$pe.se<-pp$se.fit
} else if(length(unique(d2$depth))<=3 & dim(d2)[1]>20){
modp<-gam(pe~s(day,k=5) +as.factor(year),data=d)
    pp<-predict(modp,newdata=pdat2,se.fit=TRUE)
    pdat2$pe<-pp$fit
    pdat2$pe.se<-pp$se.fit
} else {
    pdat2$pe<-NA
    pdat2$pe.se<-NA
}

pdat<-merge(pdat,pdat2,by=c('year','day','depth'),all=TRUE)
return(pdat)
}
#tdatz<-ddply(jlzoop2,.(stn),.fun=f)
#tdatp<-ddply(jlphyto2,.(stn),.fun=f)

xyplot(scount~year|stn,data=tdatp,pch=15)
xyplot(shannon~year|stn,data=tdatp,pch=15)
xyplot(pe~year|stn,data=tdatp,pch=15)
xyplot(rich~year|stn,data=tdatp,pch=15)










###############################################################

#CALCULATES ANNUAL TURNOVER RATE BETWEEN SUCCESSIVE TIME POINTS
ff<-function(d){
dt<-sort(unique(d$year))

l2<-list()
   for(j in 1:(length(dt)-1)){
       dt1<-subset(d,year==dt[j])
       dt2<-subset(d,year==dt[j+1])
     total<-length(unique(c(dt1$scientificname,dt2$scientificname)))
       #GETS SPECIES AT FIRST TIME THAT ARE NOT PRESENT IN SECOND
       emmigrants<-dim(subset(dt1,!(scientificname %in% dt2$scientificname)))[1]
       #GETS SPECIES AT FIRST TIME THAT ARE NOT PRESENT IN SECOND
       immigrants<-dim(subset(dt2,!(scientificname %in% dt1$scientificname)))[1]
       to<-(immigrants+emmigrants)/total
       yrs<-as.numeric(dt[j+1]-dt[j])

       l2[[j]]<-data.frame(year=dt[j+1],
                           to=to/yrs)
}
return(data.frame(do.call('rbind',l2)))
}
p.to<-ddply(subset(jlphyto,count>0 & depth<80),.(stn),.fun=ff,.progress='text')
z.to<-ddply(subset(jlzoop,count>0 & depth<80),.(stn),.fun=ff,.progress='text')

xyplot(turnover~year|stn,data=p.to,pch=15,type='b')
xyplot(turnover~year|stn,data=z.to,pch=15,type='b')


f<-function(d,nm,nm2){
mod<-gamm(to~as.factor(year),data=d,random=list(stn=~1))
pdat<-data.frame(year=sort(unique(d$year)))
p<-predict(mod$gam,newdata=pdat,se.fit=TRUE)
pdat$to<-p$fit
pdat$to.se<-p$se.fit
names(pdat)[2:3]<-gsub(' ','',paste(nm,'.',names(pdat)[2:3],'.',nm2))
return(pdat)
}
phytoto<-f(p.to,'ph','sabs')
zoopto<-f(z.to,'zp','sabs')

tdat2<-merge(phytoto,zoopto,by=c('year'),all=TRUE)

plank_standrews<-merge(tdat,tdat2,by=c('year'),all=TRUE)
setwd("N://cluster_2017//scratch//spera//data//finaldat_v2")
save(plank_standrews,file='plank_standrews.RData')


print(co
sort(unique(data$order))
save(








unique(jlphyto$stn)
f<-function(d){
t<-tapply(d$month,d$year,function(x) length(unique(x)))
return(data.frame(minyear=min(d$year),
                  maxyear=max(d$year),
                  maxdepth=max(d$depth),
                  nyr=length(unique(d$year)),
                  nmonths=mean(t)))
}
dum<-ddply(jlphyto,.(stn),.fun=f)


f<-function(d){
dout<-(data.frame(year=sort(unique(d$year)),
                  n=tapply(d$count,d$year,function(x) sum(x,na.rm=TRUE))))
plot(dout$year,dout$n,pch=15)
return(dout)
}
par(mfrow=c(4,4),mar=c(1,1,1,1))
dat<-ddply(jlphyto,.(stn),.fun=f)

lbls<-unique(subset(jlphyto,select=c('lon','lat','stn')))
map('world',xlim=c(-67.1,-66),ylim=c(44.5,45.1))
map('world',xlim=c(-67.5,-66),ylim=c(44,45.5))
points(jlphyto$lon,jlphyto$lat,pch=15,col='red')
text(lbls$lon,lbls$lat,labels=as.character(lbls$stn))
map.axes()

a1<-subset(dat,stn=='brandy_cove')
a2<-subset(dat,stn=='deadmans_harbour')
a3<-subset(dat,stn=='lime_kiln_bay')
a4<-subset(dat,stn=='wolves')
a<-merge(a1,a2,by=c('year'),all=TRUE)
a<-merge(a,a3,by=c('year'),all=TRUE)
a<-merge(a,a4,by=c('year'),all=TRUE)

cor(subset(a,is.na(stn.x)==FALSE,select=c('n.x','n.y','

setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/data/zoop/comm_comp/jlmartin')
write.csv(jlzoop,'zoop_jlmartin_counts_spera.csv',row.names=FALSE)

setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/data/phyto/comm_comp/jlmartin')
write.csv(jlphyto,'phyto_jlmartin_counts_spera.csv',row.names=FALSE)
