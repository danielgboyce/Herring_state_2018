library(ggplot2)
library(hexbin)
library(colorRamps)
library(scales)
library(plyr)
setwd('N:/data/Trawl data/DFO_NS_2016/tables')
vessel<-read.csv('GSVESSEL.csv',header=TRUE)
stratum<-read.csv('GSSTRATUM.csv',header=TRUE);names(stratum)[2:6]<-gsub(' ','',paste('strt.',names(stratum)[2:6]))
area2<-read.csv('GSAREA2.csv',header=TRUE)
howobt<-read.csv('GSHOWOBT.csv',header=TRUE)
wind<-read.csv('GSWIND.csv',header=TRUE)
force<-read.csv('GSFORCE.csv',header=TRUE)
curnt<-read.csv('GSCURNT.csv',header=TRUE)
xtype<-read.csv('GSXTYPE.csv',header=TRUE)
aux<-read.csv('GSAUX.csv',header=TRUE)
gear<-read.csv('GSGEAR.csv',header=TRUE)
missions<-read.csv('GSMISSIONS.csv',header=TRUE)
inf<-read.csv('GSINF.csv',header=TRUE)
species<-read.csv('GSSPECIES.csv',header=TRUE)
names(species)[1]<-'SPECIES2'
spec<-read.csv('GSSPEC.csv',header=TRUE)
change_log<-read.csv('GS_CHANGE_LOG.csv',header=TRUE)
maturity<-read.csv('GSMATURITY.csv',header=TRUE)
agemat<-read.csv('GSAGEMAT.csv',header=TRUE)
edge<-read.csv('GSEDGE.csv',header=TRUE)
ager<-read.csv('GSAGER.csv',header=TRUE)
det<-read.csv('GSDET.csv',header=TRUE)
cat<-read.csv('GSCAT.csv',header=TRUE)
#TAXONOMIC HEIRARCHIES
taxo<-read.csv('ecnasap_species_codes_taxo_heirarchies.csv',header=TRUE)
taxo<-subset(taxo,region=='NS')
names(taxo)[2]<-'SPEC'



#SAMPLE DATA
samp<-merge(missions,vessel,by=c('VESEL'),all.x=TRUE,all.y=FALSE)
samp<-merge(inf,samp,by=c('MISSION'),all.x=TRUE,all.y=FALSE)
samp<-merge(samp,stratum,by=c('STRAT'),all.x=TRUE,all.y=FALSE)


#CONVERT LON/LAT FROM DEGREES/MINS TO DECIMAL DEGREES/MINS
ltmin<-as.numeric(substr(samp$SLAT,3,4))/60
lnmin<-as.numeric(substr(samp$SLON,3,4))/60
ltdeg<-as.numeric(substr(samp$SLAT,1,2))
lndeg<-as.numeric(substr(samp$SLON,1,2))
samp$lon<-(lndeg+lnmin)*-1
samp$lat<-(ltdeg+ltmin)

#CATCHES DATA
catch<-merge(cat,species,by.x=c('SPEC'),by.y=c('CODE'),all.x=TRUE,all.y=FALSE)
catch<-merge(catch,taxo,by=c('SPEC'),all.x=TRUE,all.y=FALSE)

##################################################
##################################################
names(det)<-tolower(names(det))
a<-subset(det,spec==60)
par(mfrow=c(2,2))
plot(a$flen,a$fwt)
#a$flen<-ifelse(a$flen>55,((a$flen/10)*1.0866)+0.95632,a$flen)
a$flen<-ifelse(a$mission=='NED2016016',a$flen/10,a$flen)
plot(a$flen,a$fwt)
a$flen<-ifelse(a$mission!='NED2016016',(a$flen*1.0866)+0.95632,a$flen)
plot(a$flen,a$fwt)

a$flen<-ifelse(a$mission=='NED2016016',,a$flen)
plot(a$flen,a$fwt)

plot(a$flen,a$fwt,col=ifelse(a$mission=='NED2016016','red','black'),las=1,xlab='')

pdf('herring_rvsurvey_length_weight.pdf',height=5,width=7)
ggplot(a,aes(flen,fwt))+
    stat_binhex(bins=100,col='white',size=.3)+
  scale_fill_gradientn(colours=matlab.like(100),name = "Frequency",na.value=NA)
dev.off()


#GET ONLY HERRING OBSERVATIONS FROM CATCHES
catch<-subset(catch,comname=='herring(atlantic)' & is.na(TOTNO)==FALSE)

#LENGTHS DATA FOR INDIVIDUALS
lnth<-merge(catch,det,by=c('MISSION','SETNO','SPEC'),all.x=TRUE,all.y=FALSE)

#MERGE HERRING OBS WITH SAMPLING INFO
dat<-merge(samp,catch,by=c('MISSION','SETNO'),all.x=TRUE,all.y=FALSE)
dat<-merge(dat,gear,by=c('GEAR'),all.x=TRUE,all.y=FALSE)
dat$sid<-gsub(' ','',paste(dat$MISSION,'_',dat$SETNO))

ldat<-merge(lnth,samp,by=c('MISSION','SETNO'),all.x=TRUE,all.y=FALSE)
ldat<-merge(ldat,gear,by=c('GEAR'),all.x=TRUE,all.y=FALSE)
ldat$sid<-gsub(' ','',paste(ldat$MISSION,'_',ldat$SETNO))

library(data.table)
#MERGE CATCHES AND LENGTHS WITH SAMPLE INFO
#samp<-data.table(samp,keys=c('MISSION','SETNO'))
#catch<-data.table(catch,keys=c('MISSION','SETNO'))
#lnth<-data.table(lnth,keys=c('MISSION','SETNO'))
#dat<-merge(catch,samp,by=c('MISSION','SETNO'),all.x=TRUE,all.y=FALSE)
#ldat<-merge(lnth,samp,by=c('MISSION','SETNO'),all.x=TRUE,all.y=FALSE)

#MAKE NAMES EASIER TO DEAL WITH
names(dat)<-tolower(names(dat))
names(ldat)<-tolower(names(ldat))

#herr<-subset(dat,comname=='herring(atlantic)')
#herrl<-subset(ldat,comname=='herring(atlantic)')

#GET ALL SAMPLES WHERE HERRING ARE NOT RECORDED
#dm<-subset(dat,comname!='herring(atlantic)' & !(sid %in% unique(herr$sid)))

dat$totno<-ifelse(is.na(dat$totno)==TRUE,0,dat$totno)
dat$totwgt<-ifelse(is.na(dat$totwgt)==TRUE,0,dat$totwgt)
dat<-subset(dat,is.na(lon)==FALSE & is.na(lat)==FALSE)


setwd('N:/data/Spera/final')
write.csv(dat,'herring_weights_RV_survey_spera.csv',row.names=FALSE)
write.csv(ldat,'herring_lengths_RV_survey_spera.csv',row.names=FALSE)

setwd('N:/data/Spera/final')
a<-read.csv('herring_lengths_RV_survey_spera.csv',header=TRUE)



hadd<-subset(dat,comname=='haddock')
haddl<-subset(ldat,comname=='haddock')

#GET ALL SAMPLES WHERE HERRING ARE NOT RECORDED
dm<-subset(dat,comname!='haddock' & !(sid %in% unique(hadd$sid)))
dm$sampwgt<-0
dm$totwgt<-0
dm$totno<-0
dm$calwt<-0
hadd<-rbind.fill(hadd,dm)

setwd('N:/data/Spera/final')
write.csv(hadd,'haddock_weights_RV_survey_spera.csv',row.names=FALSE)
write.csv(haddl,'haddock_lengths_RV_survey_spera.csv',row.names=FALSE)



#VERIFIES THAT ALL 0'S HAVE BEEN ADDED
dum<-dat
length(unique(dum$sid))
length(unique(herr$sid))





library(maps)
map('world',xlim=c(-80,-50),ylim=c(40,50))
points(herr$lon,herr$lat,pch=16,col='blue')
points(dat$lon,dat$lat,pch=16,col='red')
deg<-substr(dat$SLON,1,2)


##############################################
