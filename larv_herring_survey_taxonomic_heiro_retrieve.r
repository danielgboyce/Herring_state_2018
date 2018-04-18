library(plyr)
library(rvest)

#READ IN LIST OF SPECIES AND SPECIES CODES
setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/SPERA/data/fish_larvae/larv_herring_BOFsurv_DAN')
setwd('C:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/fish_larvae/larv_herring_BOFsurv_DAN')
speclist<-read.csv('species_taxonomy_codes.csv',header=TRUE)
speclist<-subset(speclist,region=='NS')

#INSTALL PACKAGES NECESSARY TO RETRIEVE HEIRARCHIES FROM WORMS; NEEDS TO RUN ON <R 3.2
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

#FIX SOME SPECIES NAMES
speclist$sciname2<-speclist$sciname
speclist$sciname2<-tolower(speclist$sciname2)
speclist$sciname2<-gsub(' sp.','',speclist$sciname2)
speclist$sciname2<-gsub(' sp','',speclist$sciname2)
speclist$sciname2<-gsub(' p.','',speclist$sciname2)
speclist$sciname2<-gsub(' s.p.','',speclist$sciname2)
speclist$sciname2<-gsub(' s.c.','',speclist$sciname2)
speclist$sciname2<-gsub(' s.f.','',speclist$sciname2)
speclist$sciname2<-gsub(' s.o.','',speclist$sciname2)
speclist$sciname2<-gsub(' f.','',speclist$sciname2)
speclist$sciname2<-gsub(' c.','',speclist$sciname2)
speclist$sciname2<-gsub(' o.','',speclist$sciname2)
speclist$sciname2<-gsub(' (obsolete)','',speclist$sciname2)
speclist$sciname2<-gsub(' eggs','',speclist$sciname2)
speclist$sciname2<-gsub(' larvae','',speclist$sciname2)

#GET APHIA ID FROM SPECIES NAME
AphiaMatch <- function(x) {
result<-NULL
for (i in 1:length(x)) {
AphiaRecord <- iface@functions$getAphiaID(x[i],1,('http://www.marinespecies.org/aphia.php?p=soap'))
result<-c(result, AphiaRecord)
}
return(result)
}
speclist$aphiaid<-AphiaMatch(speclist$sciname2)

#ADD SOME APHIA IDS THAT WORMS CAN'T FIND
speclist$aphiaid<-ifelse(speclist$sciname=='MICROCALANUS PUSILLUS',157675 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='GNATHIA CERINA',157890 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='GEPHYREA (SIPUNCULA) P.',254787 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='MUGIL CUREMA',159416 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='CALANUS FINMARCHICUS',104464 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='EUALUS PUSIOLUS',107507 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='ERYTHROPS ERYTHROPTHALMA',852968 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='ISCHYROCEROS ANGUIPES',102412 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='PARACALANUS PARVUS',104685 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='GADOIDEI S.O.',125732 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='CENTROBRANCHUS NIGRO-OCELLATUS',126584 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='LIOPSETTA PUTNAMI',305774 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='ENCHELYOPUS/UROPHYCIS SP.',125742 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='SEBASTES FASCIATUS',127252 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='CLUPEIDAE/OSMERIDAE F.',125464 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='MYROPHIS PUNCTATUS',158643 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='clupeidae/osmeridae f.',125464,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='evermanella indica',126339 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='stephauge sp.',100763 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='actiuge sp.',100750 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='cryptotica affinis',140525 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='stephasterias albula',123808 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='lamiria sp.',138125 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='whalefishes, redmouth',125463 ,speclist$aphiaid)
speclist$aphiaid<-ifelse(speclist$sciname=='COPEPODA,NAUPLII',1080 ,speclist$aphiaid)


#SPLIT INTO ENTRIES THAT RUN THROUGH WORMS AND THOSE FOR ITIS
speclist$aphiaid<-ifelse(speclist$aphiaid==-999,NA,speclist$aphiaid)
speclistw<-subset(speclist,is.na(aphiaid)==FALSE)#GET TAXO VIA WORMS
speclisti<-subset(speclist,is.na(aphiaid)==TRUE)#GET TAXO VIA ITIS


#FUNCTION TAKES WORMS APHIA ID AND RETREIVES TAXO INFO
getFullRecord <- function(x) {
result<-NULL
for (i in 1:length(x)) {
print(x[i])
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
aa<-getFullRecord(speclistw$aphiaid)
taxolist.worms<-cbind(speclistw,aa)


#FINAL TAXONOMIC HEIRARCHY FOR WORMS
taxolist.worms<-subset(taxolist.worms,select=c('spec','sciname','comname','region','aphiaid', "rank", "valid_name","kingdom","phylum","order","family","genus",'tsn'))



#########################################################
#THOSE THAT WORMS CAN'T FIND RUN THROUGH ITIS
#########################################################
setwd('N:/data/dynamic_trophic_control/Surveys/ecnasap/species_codes_philgrayson/itis')
itis<-read.csv('ITIS_TSNs.csv',header=F,col.names=c('tsn', 'usage', 'unaccept_reason', 'parent_tsn', 'kingdom_id', 'rank_id', 'complete_name', 'rank_name', 'dir_parent_rank_id', 'req_parent_rank_id', 'completename', 'vernacular_name'),stringsAsFactors=FALSE,fileEncoding="latin1")
itis$complete_name<-tolower(itis$complete_name)
itis$completename<-tolower(itis$completename)
itis$kingdom_id<-tolower(itis$kingdom_id)
itis$rank_name<-tolower(itis$rank_name)
itis<-itis[,1:11]#DROP VERNACULAR - DUPLICATED
itis<-unique(itis)


########################################
#ADDS ITIS SPECIES NAME ACCORIDNG TO TSN
dt<-subset(speclisti,is.na(tsn)==FALSE)
dt<-merge(dt,itis,by=c('tsn'),all.x=TRUE,all.y=FALSE)
dt<-subset(dt,select=c('tsn','spec','sciname','comname','region','aphiaid','kingdom_id','rank_name','completename'))

##############################
#READS IN TAXONOMIC HEIRARCHIES
hr<-read.csv('ITIS_TSN_HIERARCHY_NEW.csv',header=TRUE,stringsAsFactors=FALSE,fileEncoding="latin1")
hr<-hr[,1:20]
hr<-subset(hr,tsn %in% speclisti$tsn)

#FUNCTION ASSEMBLES THE FULL TAXONOMIC HEIRARCHY
heirf<-function(d){
tsn<-d$tsn
print(tsn)
d2<-data.frame(t(d))
names(d2)<-'tsn'

d2$id<-seq(1,dim(d2)[1],1)
d3<-merge(d2,itis,by=c('tsn'),all.x=TRUE,all.y=FALSE)
d3<-d3[order(d3$id),]
d3<-na.omit(subset(d3,select=c('rank_name','completename')))
d4<-data.frame(t(d3))
names(d4)<-d3$rank_name
d4<-d4[2,]
d4$tsn<-tsn
return(d4)
}
out<-dlply(hr,.(tsn),.fun=heirf)
out<-rbind.fill(out)
out<-subset(out,select=c('tsn','kingdom','phylum','class','order','infraorder','family','genus','species','subspecies'))
dt<-merge(dt,out,by=c('tsn'),all.x=TRUE,all.y=FALSE)

#FULL TAXONOMIC HEIRARCHIES FROM ITIS
taxolist.itis<-subset(dt,select=c('sciname','comname','region','aphiaid',"rank_name","completename","kingdom","phylum","order","family","genus",'species','spec','tsn'))
names(taxolist.itis)<-c('sciname','comname','region','aphiaid', "rank", "valid_name","kingdom","phylum","order","family","genus",'species','spec','tsn')

#MERGE THE TWO - SOME ENTRIES STILL CAN'T BE FOUND BUT THEY ARE GARBAGE/MUD, ETC...
taxolist.full<-rbind.fill(taxolist.worms,taxolist.itis)
setwd('C:/Users/copepod/Documents/aalldocuments/literature/research/active/SPERA/data/fish_larvae/larv_herring_BOFsurv_DAN')
write.csv(taxolist.full,'BoF_larval_herring_taxon_heirarchy.csv',row.names=FALSE)






