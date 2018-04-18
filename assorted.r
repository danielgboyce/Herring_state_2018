library(RColorBrewer)
library(segmented)
library(splines)
library(strucchange)
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
figsdir<-'C://Users//copepod//Documents//aalldocuments//literature//research//active//SPERA//Figures'
figsdir<-'C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures'

#LON/LAT FOR GERMAN BANK
gblon<--66.3
gblat<-43.3


#################################################################
#VARIOUS EXPLORATORY
#################################################################


#EXAMINE SIZE/AGE SPECIFICITY OF GEAR TYPES
setwd(datadir)
herc<-read.csv('her_size_bygeartype_2006asses.csv',header=TRUE,na.strings=c('-','- ',' - '))
herc$variable<-gsub(' ','',paste(as.character(herc$variable)))
vr<-names(herc)[2:12]
l<-list()
for(i in 1:length(vr)){
d<-subset(herc,select=c(paste(vr[i]),'variable','Total','year','fishery','season'))
names(d)[1]<-'y'
d$age<-i
l[[i]]<-d
}
herc<-data.frame(do.call('rbind',l))

a<-subset(herc,variable=='%numbers')
b<-subset(herc,variable=='%catchwt.')
f<-function(d,cl){
    mod<-gam(y~s(age,k=10),data=d,gamma=.75)
    pdat<-data.frame(age=seq(1,11,length.out=100))
    p<-predict(mod,newdata=pdat)
    pdat$p<-p
    lines(pdat$age,pdat$p,col=cl,lwd=2)
    mx<-subset(pdat,p==max(pdat$p))
    points(mx$age,-.015,pch=17,col=cl,cex=1.5)
}
setwd(figsdir)
pdf('herring_catch_byage_geartype.pdf',height=8,width=7)
par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(0,0,xlim=c(1,11),ylim=c(0,.4),las=1,col='white',xaxt='n',xlab='Age',ylab='Proportion of total catch (numbers)')
legend('right',c('Purse seine','Gillnet','Weir'),col=c('firebrick4','dodgerblue4','green4'),lwd=2,bty='n')
f(subset(a,fishery=='ps'),'firebrick4')
f(subset(a,fishery=='gill'),'dodgerblue4')
f(subset(a,fishery %in% c('wier','wiers')),'green')
axis(1,seq(1,11,1))

plot(0,0,xlim=c(1,11),ylim=c(0,.4),las=1,col='white',xaxt='n',xlab='Age',ylab='Proportion of total catch (weight)')
legend('right',c('Purse seine','Gillnet','Weir'),col=c('firebrick4','dodgerblue4','green4'),lwd=2,bty='n')
f(subset(b,fishery=='ps'),'firebrick4')
f(subset(b,fishery=='gill'),'dodgerblue4')
f(subset(b,fishery=='wier'),'green')
axis(1,seq(1,11,1))
dev.off()


#######################################
a<-read.csv('her_fishing_distribution_2006assess.csv',header=TRUE,na.strings=c('-',' -',' - ','- '))
f<-function(d){
    return(data.frame(ps=sum(d$X4Xs.Fall...Winter.Purse.Seine,d$X4Xqr.Summer.Purse.Seine,d$X4W.Winter.Purse.Seine,na.rm=TRUE),
                      gill=d$X4X.Summer.Gillnet,
                      wier=sum(d$X4Xr.Nova.Scotia.Weir,d$Non.Stock.4Xs.N.B..Weir...Shutoff,na.rm=TRUE)))
}
aa<-ddply(a,.(Year),.fun=f)
aa$gill.ps<-aa$gill/aa$ps
plot(aa$Year,aa$gill.ps,pch=15,las=1,xlab='Year',ylab='Ratio gillnet:purse')



##############################################
###############################################
setwd(datadir)
load('phytoplankton_biomass_all_spera_spawnar.RData');phyt<-dat
sz<-read.csv("phyto_size_monthly_1997_2007_spera_spawnar.csv",header=TRUE)
fg<-read.csv("phyto_functional_monthly_1997_2010_spera_spawnar.csv",header=TRUE)
#strat<-read.csv("physical_stratification_spera_spawnar.csv",header=TRUE)
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw<-subset(rvw,month %in% c(6,7,8))
#vtemp<-read.csv('vtemp.csv',header=TRUE)
#vwind<-read.csv('vwind.csv',header=TRUE)
#ctdat<-read.csv('ctdat.csv',header=TRUE)
#ctdat<-read.csv('ct_data.csv',header=TRUE)
#names(ctdat)[1]<-'year'
#load("BoF_larval_herring_lengths_spera_spawnar.RData")
larvs<-dat



###########################################################
#AVERAGE DEPTH OF OCCUPANCY FROM RV SURVEY
setwd(datadir)
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw<-subset(rvw,month>=6 & month<=8)
#rvw2<-subset(rvw,is.na(depth)==FALSE & lon< -65.5)
rvw2<-subset(rvw,is.na(depth)==FALSE)
mod<-gam(depth~as.factor(year) + s(lon,lat),weights=rvw2$totno,data=rvw2,gamma=1)
her2<-data.frame(year=sort(unique(rvw2$year)),
                 lon=gblon,
                 lat=gblat)
p<-predict(mod,newdata=her2,se.fit=TRUE)
her2$her.dep<-p$fit
her2$her.dep.se<-p$se.fit
plot(her2$year,her2$her.dep,pch=15)




###########################################################
#EXAMINE THERMAL TOLERANCE RANGE OF HERRING FROM RV SURVERY
rvw<-read.csv("herring_weights_RV_survey_spera_spawnar.csv",header=TRUE)
rvw2<-subset(rvw,is.na(surface_temperature)==FALSE)#ONLY 358 MISSING
rvw2$temp<-round((rvw2$bottom_temperature+rvw2$surface_temperature)/2,digits=1)
rvw2$btemp<-round(rvw2$bottom_temperature,digits=0)
rvw2$stemp<-round(rvw2$surface_temperature,digits=0)
rvw2$occ<-ifelse(rvw2$totno>0,TRUE,FALSE)
save(rvw2,file='rvw2.RData')

stdat<-data.frame(temp=sort(unique(rvw2$temp)),
stemp=tapply(rvw2$surface_temperature,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
btemp=tapply(rvw2$bottom_temperature,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
swgt=tapply(rvw2$sampwgt,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
twgt=tapply(rvw2$totwgt,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
tno=tapply(rvw2$totno,rvw2$temp,function(x) mean(x,na.omit=TRUE)),
n=tapply(rvw2$totno,rvw2$temp,function(x) length(x)))
save(stdat,file='stdat.RData')













#AVERAGE AGE OF CATCH VERSUS LARVAL LENGTH: OLDER FISH YIELD LARGER LARVAE
#WEAK CORRELATION FOR POOLED DATA (.07), BUT STRONGER WHEN EXCLUDING OUTLYING VALUES (YEARS 1987-1994; R=0.67)
a<-merge(hlen,agedat,by=c('year'),all=FALSE)
plot(a$age.ave,a$p,pch=16,las=1,xlab='Average age',ylab='Average larvae size',xlim=c(3,5.5))
text(a$age.ave,a$p,labels=a$year,cex=1.5,adj=0)
aa<-subset(a,!(year %in% c(1988,1989,1990,1991,1992,1993,1994)))
cor(a$age.ave,a$p)
cor(aa$age.ave,aa$p)

a<-merge(hlen,agedat,by=c('year'),all=FALSE)
a<-merge(a,phen4,by=c('year'),all=FALSE)
a<-merge(a,her,by=c('year'),all=TRUE)
par(mfrow=c(2,3))
plot(a$age.ave,a$p,pch=15)
plot(a$mxfall,a$p,pch=15)
plot(a$tim,a$p,pch=15)
plot(a$dur20,a$p,pch=15)
plot(a$p,a$ssb,pch=15)
text(a$p,a$ssb,labels=a$year)
plot(a$p,log(a$ssb),pch=15)

a1<-subset(a,select=c('year','ssb'))
a2<-subset(a,select=c('year','p'))
a2$year<-a2$year+6
aa<-merge(a1,a2,by=c('year'),all=FALSE)
plot(aa$p,aa$ssb,pch=15,las=1,xlab='Average larval length (t-6)',ylab='SSB')
rv<-round(cor(aa$p,aa$ssb,use='pairwise.complete.obs'),digits=2)
legend('topleft',paste('r =',rv),bty='n')

a1<-subset(a,select=c('year','ssb'))
a2<-subset(a,select=c('year','larv'))
a2$year<-a2$year-1
aa<-merge(a1,a2,by=c('year'),all=FALSE)
plot(aa$ssb,aa$larv,pch=15,las=1,ylab='Average larval abundance',xlab='SSB (t-1)')
rv<-round(cor(log10(aa$larv),aa$ssb,use='pairwise.complete.obs'),digits=2)
legend('topleft',paste('r =',rv),bty='n')

#LAG0 R=-.06; LAG1=0.14; LAG2 R=0.25; LAG3 R=0.46; LAG4=0.67; LAG5=0.76; LAG6=0.81; LAG7=0.78; LAG8=0.76


mod<-stepAIC(lm(p~age.ave+mxfall+tim+dur20,data=a))

mod1<-lm(p~age.ave,data=a)
mod2<-lm(p~mxfall,data=a)
mod3<-lm(p~age.ave+mxfall,data=a)
mod4<-lm(p~age.ave*mxfall,data=a)
AIC(mod1,mod2,mod3,mod4)
hist(a$p,breaks=100)
z<-subset(a,year>=1988 & year<=1994)
plot(z$mxfall,z$p)

plot(hlen$year,hlen$p,pch=15,xlim=c(1970,2010),las=1,ylab='Herring larvae length',xlab='Year')
par(new=TRUE)
plot(her$year,her$larv,pch=16,col='red',ylim=c(0,100),xlim=c(1970,2010),axes=FALSE,xlab='',ylab='')
axis(4,seq(0,100,25),labels=TRUE,las=1)
a<-merge(hlen,her,by=c('year'),all=FALSE)
a<-na.omit(subset(a,select=c('year','p','larv')))
r<-cor(log(a$p),log(a$larv),use='pairwise.complete.obs')
par(fig = c(0,1,0,1))
plot(log(a$p),log(a$larv),pch=16,las=1,xlab='Herring larvae length',ylab='Herring larvae abundance')
legend('topright',legend=paste('r =',round(r,digits=2)),bty='n')

a<-merge(hlen,her,by=c('year'),all=FALSE)
a<-na.omit(subset(a,select=c('p','f')))
cor(a$p,a$f,use='pairwise.complete.obs')
ccf(a$p,a$f)
ccf(log(a$p),log(a$f))
plot(a$p,a$f)
plot(log(a$p),log(a$f))

