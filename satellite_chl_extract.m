cd('C:\Users\sailfish\Documents\aalldocuments\literature\research\active\match_mismatch_global\data\fishing')
dats=shaperead('RAM_3.shp')
plot(dats)

cd('C:\Users\sailfish\Documents\aalldocuments\literature\research\complete\2016_synchrony_scales\data')
dats=shaperead('cod_stock_polygons.shp')
plot(dats);
 plot([dats.X],[dats.Y]);
 
 x=out(:,1);
 y=out(:,2);
 xp=dats.X;
 yp=dats.Y;
 
 datss=inpolygon(x,y,xp,yp);
 out(:,4)=datss;
 
 out2=out(out(:,4)==1,:);
scatter(out2(:,1),out2(:,2),'.');
scatter(out(:,1),out(:,2),'.');


dmin = dat(dat(:,2) == minchl,:);

        


file='L3m_20020429__GLOB_4_AV-MER_CHL1_DAY_00.nc';

                lon2=ncread(file,'lon')';
                %lon2=lon2(lon2(:,1)<0,:)';
                lat2=ncread(file,'lat');
                latlen=length(lat2);
                lonlen=length(lon2);
                latout = repmat(lat2,lonlen,1);
                lonout = repmat(lon2,latlen,1);
                lonout = reshape(lonout,length(latout),1);
                lnth=length(lonout); 
                
                 a=flipud(rot90(ncread(file,'CHL1_mean')));
                a(a < 0) = NaN;%replace missing with NaN
                a = reshape(a, lnth,1);
                a(a==1) = NaN;%replace 1's with NaN; some years have weird 1's
                a=reshape(a,length(lonout),1);
                out=[lonout latout a];
                
                
                
%EXTRACT CHL SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('I:\MERIS_chl_daily')

%ncdisp('L3m_20020429__GLOB_4_AV-MER_CHL1_DAY_00.nc');
file='L3m_20020429__GLOB_4_AV-MER_CHL1_DAY_00.nc';

                lon2=ncread(file,'lon')';
                %lon2=lon2(lon2(:,1)<0,:)';
                lon2=lon2(:,2662:2810);
                lat2=ncread(file,'lat');
                lat2=lat2(1010:1200,:);
                latlen=length(lat2);
                lonlen=length(lon2);
                latout = repmat(lat2,lonlen,1);
                lonout = repmat(lon2,latlen,1);
                lonout = reshape(lonout,length(latout),1);
                lnth=length(lonout); 
      
                
subplot(2,2,3)
spy (a); figure(gcf)
aa=flipud(rot90(a));
image(a*1000);

            
cd('I:\MERIS_chl_daily\chl1\4km')
year=2003;
month=4;
day=29;
     

tic
for year=2003:2012;
    
out2= zeros(lnth*365,5);
out2(:,:) = NaN;

for month=1:12;
    
            if (month<=9); 
            eval(['month1 = ''', '0',num2str(month),'''';]);
            else
            eval(['month1 = ''', '',num2str(month),'''';]);
            end;
            
               for day=1:31;
            
               if (day<=9); 
               eval(['day1 = ''', '0',num2str(day),'''';]);
               else
               eval(['day1 = ''', '',num2str(day),'''';]);
               end;          
            
                eval(['file = ''L3m_', num2str(year),month1,day1,'__GLOB_4_AV-MER_CHL1_DAY_00.nc''';]);

                %try to extract chl but inserts NaN if file missing
                try
                  
                %calculates day of the year
                DC=datevec(datestr([year,month,day,0,0,0]));
                DV  = DC(:, 1:3);   % [N x 3] array, no time
                DV2 = DV;
                DV2(:, 2:3) = 0;    % [N x 3], day before 01.Jan
                doy = datenum(DV) - datenum(DV2);


                a=flipud(rot90(ncread(file,'CHL1_mean')));
                a=a(1010:1200,2662:2810);
                a(a < 0) = NaN;%replace missing with NaN
                %a = resizem(a, [180 360]);%resize to 5 degree
                a = reshape(a, lnth,1);
                a(a==1) = NaN;%replace 1's with NaN; some years have weird 1's
                a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                %out = out(out(:,2) > 41 & out(:,2) < 47 & out(:,1) > -68 & out(:,1) < -63,:);%restrict to div 4X
                out(:,4)=year;
                out(:,5)=doy;
                %out = single(out);%converts to single to save memory
    
                index1=((doy-1)*lnth)+1;
                index2=index1+(lnth-1);
                out2(index1:index2,:)=out;
                clear doy;
                
                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 'a' num2str(year) ' = out2;' ]);
    
end;

toc


mdata = [a2002; a2003; a2004; a2005; a2006; a2007; a2009; a2010; a2011; a2012];%2008 contains missing data so not included

scatter(mdata(:,1),mdata(:,2),'.')
scatter(a2002(:,1),a2002(:,2),'.')

cd('N:\data\chl_phenology\data')
csvwrite('meris_chla_1day_4km_spera.csv',mdata);

                
out(any(isnan(out),2),:) = [];
scatter(out(:,1),out(:,2),'.')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODIS
year=2003;
day=1;

cd('I:\modis_chl_daily')

for year=2004:2014;
    
out2= zeros(lnth*365,5);
out2(:,:) = NaN;

        for day=1:365;
            
            if (day<=9); 
            eval(['file = ''A', num2str(year),'00',num2str(day),'.L3m_DAY_CHL_chlor_a_4km''';]);
            elseif (day>9 && day<100)
            eval(['file = ''A', num2str(year),'0',num2str(day),'.L3m_DAY_CHL_chlor_a_4km''';]);
            else
            eval(['file = ''A', num2str(year),num2str(day),'.L3m_DAY_CHL_chlor_a_4km''';]);
            end;

                %try to extract chl but inserts NaN if file missing
                try
                a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                %a(a < 0) = 0;%replace missing with NaN
                % a(isnan(a)) = 0 ;
                a=a(1010:1200,2662:2810);
           
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                a(a < 0) = NaN;%replace missing with NaN
                %a = resizem(a, [180 360]);%resize to 5 degree
                a = reshape(a, lnth,1);
                %a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=day;
                out = single(out);%converts to single to save memory
    
                index1=((day-1)*lnth)+1;
                index2=index1+(lnth-1);
                out2(index1:index2,:)=out;

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 'a' num2str(year) ' = out2;' ]);
    
end;


mdata = [a2003; a2004; a2005; a2006; a2007; a2009; a2010; a2011; a2012; a2013; a2014];%2008 contains missing data so not included
cd('N:\data\chl_phenology\data')
csvwrite('modis_chla_1day_4km_spera.csv',mdata);

scatter(mdata(:,1),mdata(:,2),'.')


scatter(a2003(:,1),a2003(:,2),'.')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('C:\Users\sailfish\Downloads\swfs')
cd('I:\SWFS_chl_daily')
file='S1998001.L3m_DAY_CHL_chlor_a_9km'
A2014365.L3m_DAY_CHL_chlor_a_4km
subplot(2,2,3)
spy (a); figure(gcf)

%can't figure out how to get lon/lat from hdf file (below), so instead read in 9km coords from netcdf file 
ncdisp('N:\Files\Data\chl\remote.sensing\giovanni\seawifs_chl\S1997 (9).nc');
lat2=flipud(ncread('N:\Files\Data\chl\remote.sensing\giovanni\seawifs_chl\S1997 (9).nc','G3fakeDim0'));
lon2=ncread('N:\Files\Data\chl\remote.sensing\giovanni\seawifs_chl\S1997 (9).nc','G3fakeDim1')';
lon2=lon2(:,1330:1405);
lat2=lat2(515:590,:);
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);
lnth=length(lonout); 


for year=1998:2008;
    
out2= zeros(lnth*365,5);
out2(:,:) = NaN;

        for day=1:365;
            
            if (day<=9); 
            eval(['file = ''S', num2str(year),'00',num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            elseif (day>9 && day<100)
            eval(['file = ''S', num2str(year),'0',num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            else
            eval(['file = ''S', num2str(year),num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            end;

                %try to extract chl but inserts NaN if file missing
                try
                a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});
                % a(isnan(a)) = 0 ;
                
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                a(a < 0) = NaN;%replace missing with NaN
                a=a(515:590,1330:1405);
                a = reshape(a, lnth,1);
                %a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=day;
                out = single(out);%converts to single to save memory
    
                index1=((day-1)*lnth)+1;
                index2=index1+(lnth-1);
                out2(index1:index2,:)=out;

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 's' num2str(year) ' = out2;' ]);
    
end;


sdata = [s1998; s1999; s2000; s2001; s2002; s2000; s2003; s2004; s2005; s2006; s2007];%2008 contains missing data so not included
cd('N:\data\chl_phenology\data')
csvwrite('swfs_chla_1day_9km_spera.csv',sdata);

scatter(sdata(:,1),sdata(:,2),'.')
scatter(s1998(:,1),s1998(:,2),'.')








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('I:\CZCZ_chl_daily\mapped')
a=l3m_data;
a(a < 0) = NaN;%replace missing with NaN
a(isnan(a)) = 0 ;
spy (a); figure(gcf)
%a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
   
%can't figure out how to get lon/lat from hdf file (below), so instead read in 9km coords from netcdf file 
lat2=flipud(ncread('N:\Files\Data\chl\remote.sensing\giovanni\seawifs_chl\S1997 (9).nc','G3fakeDim0'));
lon2=ncread('N:\Files\Data\chl\remote.sensing\giovanni\seawifs_chl\S1997 (9).nc','G3fakeDim1')';
lon2=lon2(:,1330:1405);
lat2=lat2(515:590,:);
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);
lnth=length(lonout); 


cd('I:\CZCZ_chl_daily\mapped')
for year=1978:1986;
    
out2= zeros(lnth*365,5);
out2(:,:) = NaN;

        for day=1:365;
            
            if (day<=9); 
            eval(['file = ''C', num2str(year),'00',num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            elseif (day>9 && day<100)
            eval(['file = ''C', num2str(year),'0',num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            else
            eval(['file = ''C', num2str(year),num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            end;

                %try to extract chl but inserts NaN if file missing
                try
                a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});
                a=a(515:590,1330:1405);
                
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                a(a < 0) = NaN;%replace missing with NaN
                
                if (sum(sum(isnan(a)))<9200000)
                a = reshape(a, lnth,1);
                %a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=day;
                out = single(out);%converts to single to save memory
    
                index1=((day-1)*lnth)+1;
                index2=index1+(lnth-1);
                out2(index1:index2,:)=out;
                
                else display('sparse matrix');
                end;

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 'c' num2str(year) ' = out2;' ]);
    
end;


cdata = [c1978; c1979; c1980; c1981; c1982; c1983; c1984; c1985; c1986];
cd('N:\data\chl_phenology\data')
csvwrite('czcs_chla_1day_9km_spera.csv',cdata);

cdata = c1986;
scatter(cdata(:,1),cdata(:,2),'.')









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%extracts alvain 2005 frequency of occurence for individual phytoplankotn
%%functional groups for specified region 4X averaged monthly  from 1998 to 2010
%date is in gregorian days and some months are missing , so needed to
%modify code to deal with this
cd('N:\Files\Data\phytoplankton\physat\monthly_frequency_groups')

year=2008;
month=1;


for year=1997;

%read data in
eval(['file = ''PHYSAT_', num2str(year),'_freq-pro.nc''';]);
eval(['file2 = ''PHYSAT_', num2str(year),'_freq-slc.nc''';]);
eval(['file3 = ''PHYSAT_', num2str(year),'_freq-dia.nc''';]);

lon = getnc(file,'lon')';
lat = flipud(getnc(file,'lat'));
lon=lon(:,1330:1405);
lat=lat(515:590,:);
latlen=length(lat);
lonlen=length(lon);
latout = repmat(lat,lonlen,1);
lonout = repmat(lon,latlen,1);
lonout = reshape(lonout,length(latout),1);
lnth=length(lonout); 

pro = getnc(file,'freq-pro',-1,-1,1);
pro(pro < -9.9990) = NaN;

slc = getnc(file2,'freq-slc',-1,-1,1);
slc(slc < -9.9990) = NaN;

dia = getnc(file3,'freq-dia',-1,-1,1);
dia(dia < -9.9990) = NaN;

dt = getnc(file, 'date');
%montha=datestr(dt, 'mm');


    for i = 1:length(dt);
        
        %month=str2num(montha(i,:));
        %display(month);

        prob=pro(i,:,:);
        prob = squeeze(prob);
        prob=flipud(prob);
        prob=prob(515:590,1330:1405);
        prob = reshape(prob, lnth,1);
        
        diab=dia(i,:,:);
        diab = squeeze(diab);
        diab=flipud(diab);
        diab=diab(515:590,1330:1405);
        diab = reshape(diab, lnth,1);
        
        slcb=slc(i,:,:);
        slcb = squeeze(slcb);
        slcb=flipud(slcb);
        slcb=slcb(515:590,1330:1405);
        slcb = reshape(slcb, lnth,1);
        
        out=[lonout latout prob diab slcb];
        
        month=str2num(datestr(dt(i,:), 'mm'));
        display(month);
        out(:,6)=month;
                
        eval([ 'out' num2str(month) ' = out;' ]);

    end;
       
    outt = [out1;out2;out3;out4;out5;out6;out7;out8;out9;out10;out11;out12];
    outt(any(isnan(outt),2),:) = [];
    outt(:,7)=year;
    eval([ 'out' num2str(year) ' = outt;' ]);

clear slc dia pro;

end;

               
out = [out1997; out1998; out1999; out2000; out2001;out2002;out2003;out2004;out2005;out2006;out2007;out2008;out2009;out2010];


scatter(out(:,1),out(:,2));

csvwrite('alvain_dominance__monthly_9km_1997_2010_spera.csv',out);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%extracts brewin phytoplankton size classes from 1997 to 2007
%References  = 'Brewin, R.J.W., Sathyendranath, S., Hirata, T., Lavender, S., Barciela, R.M. & Hardman- Mountford, N.J. (2010). A three-component model of phytoplankton size class for the Atlantic Ocean. Ecological Modelling, 221, 1472−1483.'

cd('C:\Users\sailfish\Downloads\brewin')
file='SeaWiFS_GLOBAL_19970901_Phytoplankton-Size-Class.nc';
ncdisp(file);

lon2 = getnc(file,'longitude',-1,-1,1);
lat2 = flipud(getnc(file,'latitude',-1,-1,1));
lon2=lon2(1,1330:1405);
lat2=lat2(515:590,1);
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);
lnth=length(lonout); 


for year=1997:2007;

        for month=1:12;
          
             if (month<=9); 
            eval(['month1 = ''', '0',num2str(month),'''';]);
            else
            eval(['month1 = ''', '',num2str(month),'''';]);
            end;
            
%read data in
eval(['file = ''SeaWiFS_GLOBAL_', num2str(year),num2str(month1),'01_Phytoplankton-Size-Class.nc''';]);   

try 
    
pic = flipud(getnc(file,'Pico_percent_Tchl',-1,-1,1));
nano = flipud(getnc(file,'Nano_percent_Tchl',-1,-1,1));
mic = flipud(getnc(file,'Micro_percent_Tchl',-1,-1,1));
pic=pic(515:590,1330:1405);
nano=nano(515:590,1330:1405);
mic=mic(515:590,1330:1405);
pic(pic < -9.9990) = NaN;
nano(nano < -9.9990) = NaN;
mic(mic < -9.99) = NaN;

pic = reshape(pic, lnth,1);
nano = reshape(nano, lnth,1);
mic = reshape(mic, lnth,1);
        
out=[lonout latout pic nano mic];
        
        %month=str2num(datestr(dt(i,:), 'mm'));
        display(month);
        out(:,6)=month;
                
        eval([ 'out' num2str(month) ' = out;' ]);
        
catch ME
    out=[NaN NaN NaN NaN NaN NaN];
        eval([ 'out' num2str(month) ' = out;' ]);
end;

    end;
       
    outt = [out1;out2;out3;out4;out5;out6;out7;out8;out9;out10;out11;out12];
    outt(any(isnan(outt),2),:) = [];
    outt(:,7)=year;
    eval([ 'out' num2str(year) ' = outt;' ]);
end;


out = [out1997; out1998; out1999; out2000; out2001;out2002;out2003;out2004;out2005;out2006;out2007];

scatter(out(:,1),out(:,2));

csvwrite('brewin_phyto_sizegroups_monthly_9km_1997_2007_spera.csv',out);

    








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%import meris pp data extracted from http://www.science.oregonstate.edu/ocean.productivity/standard.product.php

year=2003;
day=1;

ncdisp('N:\Files\Data\chl\remote.sensing\giovanni\seawifs_chl\S1997 (9).nc');
lat2=flipud(ncread('N:\Files\Data\chl\remote.sensing\giovanni\seawifs_chl\S1997 (9).nc','G3fakeDim0'));
lon2=ncread('N:\Files\Data\chl\remote.sensing\giovanni\seawifs_chl\S1997 (9).nc','G3fakeDim1')';
lon2=lon2(:,1330:1405);
lat2=lat2(515:590,:);
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);
lnth=length(lonout); 

cd('C:\Users\sailfish\Downloads\meris_pp')

for year=2002:2016;
    
out2= zeros(lnth*365,5);
out2(:,:) = NaN;

        for day=1:365;
            
            if (day<=9); 
            eval(['file = ''vgpm.', num2str(year),'00',num2str(day),'.hdf''';]);
            elseif (day>9 && day<100)
            eval(['file = ''vgpm.', num2str(year),'0',num2str(day),'.hdf''';]);
            else
            eval(['file = ''vgpm.', num2str(year),num2str(day),'.hdf''';]);
            end;

                %try to extract chl but inserts NaN if file missing
                try
                a = hdfread(file, '/npp', 'Index', {[1  1],[1  1],[2160  4320]});
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});
                % a(isnan(a)) = 0 ;
                %image(a);
                
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                a(a < 0) = NaN;%replace missing with NaN
                a=a(515:590,1330:1405);
                a = reshape(a, lnth,1);
                %a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=day;
                out = single(out);%converts to single to save memory
    
                index1=((day-1)*lnth)+1;
                index2=index1+(lnth-1);
                out2(index1:index2,:)=out;

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 's' num2str(year) ' = out2;' ]);
    
end;


sdata = [s2002; s2003; s2004; s2005; s2006; s2007; s2008; s2009; s2010; s2011; s2012; s2013; s2014; s2015; s2016];%2008 contains missing data so not included
cd('N:\data\chl_phenology\data')
csvwrite('meris_pp_8day_9km_spera.csv',sdata);

scatter(sdata(:,1),sdata(:,2));



























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%czcs 8 day data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('I:\CZCS_chl_8day')

l3m_data = hdfread('I:\CZCZ_chl_daily\mapped\C1978303.L3m_DAY_CHL_chlor_a_9km', '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});
l3m_data = hdfread('I:\CZCS_chl_8day\C19782971978304.L3m_8D_CHL_chlor_a_9km', '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});

a=l3m_data;
a(a < 0) = NaN;%replace missing with NaN
a(isnan(a)) = 0 ;
spy (a); figure(gcf)

%a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});

a(a < 0) = NaN;%replace missing with NaN
                
subplot(2,2,3)
spy (a); figure(gcf)



   
year=1978;
day=305;

for year=1978:1986;
    
out2= zeros(lnth*365,5);
out2(:,:) = NaN;

        for day=1:365;
            
            
            if (day<=9); 
            eval(['day1 = ''', '00',num2str(day),'''';]);
            eval(['day2 = ''', '00',num2str(day+7),'''';]);
            elseif (day>9 && day<100)   
            eval(['day1 = ''', '0',num2str(day),'''';]);
            eval(['day2 = ''', '0',num2str(day+7),'''';]);
            else
            eval(['day1 = ''', '',num2str(day),'''';]);
            eval(['day2 = ''', '',num2str(day+7),'''';]);
            end;
            
            eval(['file = ''C', num2str(year),day1,num2str(year),day2,'.L3m_8D_CHL_chlor_a_9km''';]);
                        

                %try to extract chl but inserts NaN if file missing
                try
                a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});
                % a(isnan(a)) = 0 ;
                
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                a(a < 0) = NaN;%replace missing with NaN
                
                if (sum(sum(isnan(a)))<9200000)
                a = resizem(a, [180 360]);%resize to 5 degree
                a = reshape(a, lnth,1);
                a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=day;
                out = single(out);%converts to single to save memory
    
                index1=((day-1)*lnth)+1;
                index2=index1+(lnth-1);
                out2(index1:index2,:)=out;
                
                else display('sparse matrix');
                end;

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 'c' num2str(year) ' = out2;' ]);
    
end;

cdata = [c1978];

cdata = [c1978; c1979; c1980; c1981; c1982; c1983; c1984; c1985; c1986];
cd('M:\data\chl_phenology\data')
csvwrite('czcs_chla_8day_1deg.csv',cdata);

cdata = c1986;
scatter(cdata(:,1),cdata(:,2),'.')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5 DEGREE



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('I:\MERIS_chl_daily')
cd('C:\Users\sailfish\Downloads\meris')

ncdisp('L3m_20020429__GLOB_4_AV-MER_CHL1_DAY_00.nc');
file='L3m_20020429__GLOB_4_AV-MER_CHL1_DAY_00.nc';

subplot(2,2,3)
spy (a); figure(gcf)
aa=flipud(rot90(a));
image(aa);

            
cd('I:\MERIS_chl_daily\chl1\4km')
year=2003;
month=4;
day=29;
     
lat2 = [87.5:-5:-87.5]';
lon2 = [-177.5:5:177.5];
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);




tic
for year=2002:2012;
    
out2= zeros(2592*365,5);
out2(:,:) = NaN;

for month=1:12;
    
            if (month<=9); 
            eval(['month1 = ''', '0',num2str(month),'''';]);
            else
            eval(['month1 = ''', '',num2str(month),'''';]);
            end;
            
               for day=1:31;
            
               if (day<=9); 
               eval(['day1 = ''', '0',num2str(day),'''';]);
               else
               eval(['day1 = ''', '',num2str(day),'''';]);
               end;          
            
                eval(['file = ''L3m_', num2str(year),month1,day1,'__GLOB_100_AV-MER_CHL1_DAY_00.nc''';]);

                %try to extract chl but inserts NaN if file missing
                try
                  
                %calculates day of the year
                DC=datevec(datestr([year,month,day,0,0,0]));
                DV  = DC(:, 1:3);   % [N x 3] array, no time
                DV2 = DV;
                DV2(:, 2:3) = 0;    % [N x 3], day before 01.Jan
                doy = datenum(DV) - datenum(DV2);

                a=flipud(rot90(ncread(file,'CHL1_mean')));
                a(a < 0) = NaN;%replace missing with NaN
                a = resizem(a, [36 72]);%resize to 5 degree
                a = reshape(a, 2592,1);
                a(a==1) = NaN;%replace 1's with NaN; some years have weird 1's
                a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=doy;
                %out = single(out);%converts to single to save memory
    
                index1=((doy-1)*2592)+1;
                index2=index1+2591;
                out2(index1:index2,:)=out;
                clear doy;
                
                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 'a' num2str(year) ' = out2;' ]);
    
end;

toc


mdata = [a2002; a2003; a2004; a2005; a2006; a2007; a2009; a2010; a2011; a2012];%2008 contains missing data so not included

scatter(mdata(:,1),mdata(:,2),'.')

cd('C:\Users\sailfish\Documents\aalldocuments\literature\postdoc_2013\chl_phenology\data')
csvwrite('meris_chla_1day_5deg.csv',mdata);

                














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODIS
lat2 = [87.5:-5:-87.5]';
lon2 = [-177.5:5:177.5];
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);

   
for year=2003:2014;
    
out2= zeros(2592*365,5);
out2(:,:) = NaN;

        for day=1:365;
            
            if (day<=9); 
            eval(['file = ''A', num2str(year),'00',num2str(day),'.L3m_DAY_CHL_chlor_a_4km''';]);
            elseif (day>9 && day<100)
            eval(['file = ''A', num2str(year),'0',num2str(day),'.L3m_DAY_CHL_chlor_a_4km''';]);
            else
            eval(['file = ''A', num2str(year),num2str(day),'.L3m_DAY_CHL_chlor_a_4km''';]);
            end;

                %try to extract chl but inserts NaN if file missing
                try
                a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                %a(a < 0) = 0;%replace missing with NaN
                % a(isnan(a)) = 0 ;
                
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                a(a < 0) = NaN;%replace missing with NaN
                a = resizem(a, [36 72]);%resize to 5 degree
                a = reshape(a, 2592,1);
                a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=day;
                out = single(out);%converts to single to save memory
    
                index1=((day-1)*2592)+1;
                index2=index1+2591;
                out2(index1:index2,:)=out;

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 'a' num2str(year) ' = out2;' ]);
    
end;


mdata = [a2003; a2004; a2005; a2006; a2007; a2009; a2010; a2011; a2012; a2013; a2014];%2008 contains missing data so not included
cd('C:\Users\sailfish\Documents\aalldocuments\literature\postdoc_2013\chl_phenology\data')
csvwrite('modis_chla_1day_5deg.csv',mdata);

scatter(mdata(:,1),mdata(:,2),'.')

data = [a2003; a2004];
scatter(data(:,1),data(:,2),'.')

data=out2;

out(any(isnan(out),2),:) = [];
scatter(out(:,1),out(:,2),'.')

scatter(out2(:,1),out2(:,2),'.')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('I:\SWFS_chl_daily')
year=1998;
day=1;
day=186;

file='S1998001.L3m_DAY_CHL_chlor_a_9km'
A2014365.L3m_DAY_CHL_chlor_a_4km

subplot(2,2,3)
spy (a); figure(gcf)

lat2 = [87.5:-5:-87.5]';
lon2 = [-177.5:5:177.5];
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);

   

for year=2008;
    
out2= zeros(2592*365,5);
out2(:,:) = NaN;

        for day=95;
            
            if (day<=9); 
            eval(['file = ''S', num2str(year),'00',num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            elseif (day>9 && day<100)
            eval(['file = ''S', num2str(year),'0',num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            else
            eval(['file = ''S', num2str(year),num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            end;

                %try to extract chl but inserts NaN if file missing
                try
                a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});
                % a(isnan(a)) = 0 ;
                
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                a(a < 0) = NaN;%replace missing with NaN
                a = resizem(a, [36 72]);%resize to 5 degree
                a = reshape(a, 2592,1);
                a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=day;
                out = single(out);%converts to single to save memory
    
                index1=((day-1)*2592)+1;
                index2=index1+2591;
                out2(index1:index2,:)=out;

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 's' num2str(year) ' = out2;' ]);
    
end;


sdata = [s1998; s1999; s2000; s2001; s2002; s2000; s2003; s2004; s2005; s2006; s2007; s2009];%2008 contains missing data so not included
cd('C:\Users\sailfish\Documents\aalldocuments\literature\postdoc_2013\chl_phenology\data')
csvwrite('swfs_chla_1day_5deg.csv',sdata);

scatter(sdata(:,1),sdata(:,2),'.')
scatter(s2010(:,1),s2010(:,2),'.')
2008
data=out2;













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('I:\CZCZ_chl_daily\mapped')

l3m_data = hdfread('I:\CZCZ_chl_daily\mapped\C1978303.L3m_DAY_CHL_chlor_a_9km', '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});

a=l3m_data;
a(a < 0) = NaN;%replace missing with NaN
a(isnan(a)) = 0 ;
spy (a); figure(gcf)

%a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});

a(a < 0) = NaN;%replace missing with NaN
                
subplot(2,2,3)
spy (a); figure(gcf)

lat2 = [87.5:-5:-87.5]';
lon2 = [-177.5:5:177.5];
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);

   

for year=1986;
    
out2= zeros(2592*365,5);
out2(:,:) = NaN;

        for day=1:365;
            
            if (day<=9); 
            eval(['file = ''C', num2str(year),'00',num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            elseif (day>9 && day<100)
            eval(['file = ''C', num2str(year),'0',num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            else
            eval(['file = ''C', num2str(year),num2str(day),'.L3m_DAY_CHL_chlor_a_9km''';]);
            end;

                %try to extract chl but inserts NaN if file missing
                try
                a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});
                % a(isnan(a)) = 0 ;
                
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                a(a < 0) = NaN;%replace missing with NaN
                
                if (sum(sum(isnan(a)))<9200000)
                a = resizem(a, [36 72]);%resize to 5 degree
                a = reshape(a, 2592,1);
                a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=day;
                out = single(out);%converts to single to save memory
    
                index1=((day-1)*2592)+1;
                index2=index1+2591;
                out2(index1:index2,:)=out;
                
                else display('sparse matrix');
                end;

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 'c' num2str(year) ' = out2;' ]);
    
end;


cdata = [c1978; c1979; c1980; c1981; c1982; c1983; c1984; c1985; c1986];
cd('C:\Users\sailfish\Documents\aalldocuments\literature\postdoc_2013\chl_phenology\data')
csvwrite('czcs_chla_1day_5deg.csv',cdata);

cdata = c1986;
scatter(cdata(:,1),cdata(:,2),'.')


n = 10;
m = 5;
a = randn(n, m);
a(rand(numel(a), 1) < 1) = nan;
aa=resizem(a,[1,1]);









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%czcs 8 day data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('I:\CZCS_chl_8day')

l3m_data = hdfread('I:\CZCZ_chl_daily\mapped\C1978303.L3m_DAY_CHL_chlor_a_9km', '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});
l3m_data = hdfread('I:\CZCS_chl_8day\C19782971978304.L3m_8D_CHL_chlor_a_9km', '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});

a=l3m_data;
a(a < 0) = NaN;%replace missing with NaN
a(isnan(a)) = 0 ;
spy (a); figure(gcf)

%a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});

a(a < 0) = NaN;%replace missing with NaN
                
subplot(2,2,3)
spy (a); figure(gcf)

lat2 = [87.5:-5:-87.5]';
lon2 = [-177.5:5:177.5];
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);

   
year=1978;
day=305;

for year=1978:1986;
    
out2= zeros(2592*365,5);
out2(:,:) = NaN;

        for day=1:365;
            
            
            if (day<=9); 
            eval(['day1 = ''', '00',num2str(day),'''';]);
            eval(['day2 = ''', '00',num2str(day+7),'''';]);
            elseif (day>9 && day<100)   
            eval(['day1 = ''', '0',num2str(day),'''';]);
            eval(['day2 = ''', '0',num2str(day+7),'''';]);
            else
            eval(['day1 = ''', '',num2str(day),'''';]);
            eval(['day2 = ''', '',num2str(day+7),'''';]);
            end;
            
            eval(['file = ''C', num2str(year),day1,num2str(year),day2,'.L3m_8D_CHL_chlor_a_9km''';]);
                        

                %try to extract chl but inserts NaN if file missing
                try
                a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[2160  4320]});
                % a(isnan(a)) = 0 ;
                
                %a = hdfread(file, '/l3m_data', 'Index', {[1  1],[1  1],[4320  8640]});
                a(a < 0) = NaN;%replace missing with NaN
                
                if (sum(sum(isnan(a)))<9200000)
                a = resizem(a, [36 72]);%resize to 5 degree
                a = reshape(a, 2592,1);
                a=reshape(a,length(lonout),1);
    
                out=[lonout latout a];
                out(:,4)=year;
                out(:,5)=day;
                out = single(out);%converts to single to save memory
    
                index1=((day-1)*2592)+1;
                index2=index1+2591;
                out2(index1:index2,:)=out;
                
                else display('sparse matrix');
                end;

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
          
        end;
        
out2(any(isnan(out2),2),:) = [];
eval([ 'c' num2str(year) ' = out2;' ]);
    
end;

cdata = [c1978];

cdata = [c1978; c1979; c1980; c1981; c1982; c1983; c1984; c1985; c1986];
cd('C:\Users\sailfish\Documents\aalldocuments\literature\postdoc_2013\chl_phenology\data')
csvwrite('czcs_chla_8day_5deg.csv',cdata);

cdata = c1986;
scatter(cdata(:,1),cdata(:,2),'.')


































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extracts czcs 8day chl data to 5 degree resolution
cd('I:\Data\phytoplankton\remote.sensing\ocean_toolbox\czcs_chl\8day\mat_files')
cd('E:\Data\phytoplankton\remote.sensing\ocean_toolbox\czcs_chl\8day\mat_files')

n=1;

lat2 = [87.5:-5:-87.5]';
lon2 = [-177.5:5:177.5];
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);


outp=[];

for n=1:341;
    
        eval(['file = ''czcs_8day (', num2str(n),').mat''';]);
       
load(file);
chl=D.Chlor;
time=D.user_friendly_time;
year = str2num(time(1:4));
month = str2num(time(6:7));
day = str2num(time(9:10));
dayofyear = ((month*30.33)+day)-30.33;

%lon=D.longitude;
%lat=D.latitude;

display(time);
chl = resizem(chl, [36,72]);
chl=reshape(chl,2592,1);

out = [lonout latout chl];
out(:,4) = year;
out(:,5) = month;
out(:,6) = day;
out(:,7) = dayofyear;

out(any(isnan(out),2),:) = [];
out=single(out);
outp = vertcat(outp, out);

end;

eval(['save ', 'czcs_5deg_8day', ' outp -ascii -tabs'])








            
%reads in hermes GSM chl-a data and determines seasonal min/max
cd('C:\Users\sailfish\Documents\data\phytoplankton\remote.sensing\hermes\daily_100km')
year=2000
month=5
day=20
%for year=1998:2000;
%for year=2001:2004;
for year=2005:2008;

        for month=1:12;
                       
            if ((month==1)); daya=1:31;
            elseif ((month==2)); daya=1:28;
            elseif ((month==3 || month==5 || month==7 || month==8 || month==10 || month==12)); daya=1:31;
            else daya=1:30;
            end;
           
                    for day=min(daya):max(daya);
                        
         if (year<2002) && (month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_AV-SWF_CHL1_DAY_00.nc''';]);
        elseif (year<2002) && (month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_AV-SWF_CHL1_DAY_00.nc''';]);
        elseif (year<2002) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_AV-SWF_CHL1_DAY_00.nc''';]); 
        elseif (year<2002) && (month>=10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_AV-SWF_CHL1_DAY_00.nc''';]); 
                
        elseif (year==2002) && (month<=3) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_AV-SWF_CHL1_DAY_00.nc''';]);
        elseif (year==2002) && (month<=3) && (day>=10)            
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_AV-SWF_CHL1_DAY_00.nc''';]);
               
        %read in swf/mer data
        elseif (year==2002) && (month>=4 && month<=5) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_AVW-MERSWF_CHL1_DAY_00.nc''';]); 
        elseif (year==2002) && (month>=4 && month<=5) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_AVW-MERSWF_CHL1_DAY_00.nc''';]); 
                         
        %read in swf/mer/mod data
        elseif (year==2002) && (month>=6 && month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_AVW-MERMODSWF_CHL1_DAY_00.nc''';]); 
        elseif (year==2002) && (month>=6 && month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_AVW-MERMODSWF_CHL1_DAY_00.nc''';]); 
        
         elseif (year==2002) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_AVW-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year==2002) && (month>=10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_AVW-MERMODSWF_CHL1_DAY_00.nc''';]);
        
        elseif (year>=2003) && (month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_AVW-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year>=2003) && (month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_AVW-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year>=2003) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_AVW-MERMODSWF_CHL1_DAY_00.nc''';]);
        else
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_AVW-MERMODSWF_CHL1_DAY_00.nc''';]);        
        end;
        
          %try to extract chl but inserts NaN if file missing
          try
          c=netcdf.open(file,'NC_WRITE');
          chl = getnc(file,'CHL1_mean');
          catch ME
          display('missing');
          chl = zeros(180,360);
          chl(:,:) = NaN;
          end;
          
          
chl = reshape(chl, 64800,1);

eval([ 'a' num2str(day) ' = chl;' ]);

       end;
                    
  if (month==2);  aa=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28];
  elseif (month==4 || month==6 || month==9 || month==11); aa=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28 a29 a30];     
  else aa=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28 a29 a30 a31];
  end;
                                            
  eval([ 'm' num2str(month) ' = aa;' ]);
 
        end;
        
        aaa=[m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12];
       
          
   %now have to eliminate rows where data is missing
   %smooth series and select min and max
   out= zeros(64800,3);
   out(:,:) = NaN;

   for rowr=1:64800
              
       dat=aaa(rowr,:); 
       %dat=aaa(2221,:); 
       i = find(isnan(dat));
       x = length(i);
       
       
       %if (x>190); out(rowr,1:2)=NaN;
       if (x>275); out(rowr,1:2)=NaN;
       
       else    
       y=dat';     
       yday=[1:1:365]';
       dat=[yday y];
       dat(any(isnan(dat),2),:) = [];
       n=length(dat(:,2));
       dat=interp1(dat(:,1),dat(:,2),yday(:,1),'cubic');
       window = 10;
       dat=filter(ones(1,window)/window,1,dat);
       dat=[yday dat];
       dat=dat(6:359,:);
       maxchl=max(dat(:,2));
       minchl=min(dat(:,2));
       dmax = dat(dat(:,2) == maxchl,:);
       dmax = dmax(1,1);
       dmin = dat(dat(:,2) == minchl,:);
       dmin = dmin(1,1);   
       out(rowr,1)=dmin;
       out(rowr,2)=dmax;
       out(rowr,3)=n;
       end;
   
   end;

   latout = [89.5:-1:-89.5]';
   latout = repmat(latout,360,1);
   %latout = latout*(-1);
   lonout = [-179.5:1:179.5];
   lonout = repmat(lonout,180,1);
   lonout = reshape(lonout,64800,1);
        
    out(:,4) = lonout;
    out(:,5) = latout;
    out(:,3)=year;
    
    out(any(isnan(out),2),:) = [];
    
    eval([ 'a' num2str(year) ' = out;' ]);

           
end;



save('d1998', 'a1998');
save('d1999', 'a1999');
save('d2000', 'a2000');
save('d2001', 'a2001');
save('d2002', 'a2002');
save('d2003', 'a2003');
save('d2004', 'a2004');
save('d2005', 'a2005');
save('d2006', 'a2006');
save('d2007', 'a2007');
save('d2008', 'a2008');

load('d1998.mat');
load('d1999.mat');
load('d2000.mat');
load('d2001.mat');
load('d2002.mat');
load('d2003.mat');
load('d2004.mat');
load('d2005.mat');
load('d2006.mat');
load('d2007.mat');
load('d2008.mat');
load('d2000.mat');

swf_phen = [a1998 ; a1999; a2000; a2001; a2002; a2003; a2004; a2005; a2006; a2007; a2008];
save('swf_phen', 'swf_phen','-ascii','-tabs');

    eval(['save ', 'coords_1_16', ' coords -ascii -tabs'])
a = [a1997; a1998; a1999; a2000; a2001; a2002; a2003];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



stem(dat)   


          
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reads in hermes GSM chl-a data and determines seasonal min/max
       cd('C:\Users\sailfish\Documents\data\phytoplankton\remote.sensing\hermes\daily_100km')


for year=1998:1999;
    
        for month=1:12;
                       
            if ((month==1)); daya=1:31;
            elseif ((month==2)); daya=1:28;
            elseif ((month==3 || month==5 || month==7 || month==8 || month==10 || month==12)); daya=1:31;
            else daya=1:30;
            end;
           
                    for day=min(daya):max(daya);
                        
         if (year<=2001) && (month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]);
        elseif (year<=2001) && (month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]);
        elseif (year<=2001) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]); 
        elseif (year<=2001) && (month>=10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]); 
                
        elseif (year==2002) && (month<=3) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]);
        elseif (year==2002) && (month<=3) && (day>=10)            
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]);
               
        %read in swf/mer data
        elseif (year==2002) && (month>=4 && month<=5) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERSWF_CHL1_DAY_00.nc''';]); 
        elseif (year==2002) && (month>=4 && month<=5) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-MERSWF_CHL1_DAY_00.nc''';]); 
                         
        %read in swf/mer/mod data
        elseif (year==2002) && (month>=6 && month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]); 
        elseif (year==2002) && (month>=6 && month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]); 
        
         elseif (year==2002) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year==2002) && (month>=10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
        
        elseif (year>2003) && (month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year>2003) && (month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year>2003) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
        else
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);        
        end;
        
          %try to extract chl but inserts NaN if file missing
          try
          c=netcdf.open(file,'NC_WRITE');
          chl = getnc(file,'CHL1_mean');
          catch ME
          display('missing');
          chl = zeros(180,360);
          chl(:,:) = NaN;
          end;

chl = reshape(chl, 64800,1);

eval([ 'a' num2str(day) ' = chl;' ]);

       end;
                    
  if (month==2);  aa=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28];
  elseif (month==4 || month==6 || month==9 || month==11); aa=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28 a29 a30];     
  else aa=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28 a29 a30 a31];
  end;
                                            
  eval([ 'm' num2str(month) ' = aa;' ]);
 
        end;
        
        aaa=[m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12];
       
          
   %now have to eliminate rows where data is missing
   %smooth series and select min and max
   out= zeros(64800,2);
   %out(:,:) = NaN;

   for rowr=1:64800
              
       dat=aaa(rowr,:); 
       %dat=aaa(2221,:); 
       i = find(isnan(dat));
       x = length(i);
       
       
       if (x>190); out(rowr,1:2)=NaN;
       %if (x>275); out(rowr,1:2)=NaN;
       
       else    
       y=dat';     
       yday=[1:1:365]';
       dat=[yday y];
       dat(any(isnan(dat),2),:) = [];
       dat=interp1(dat(:,1),dat(:,2),yday(:,1),'cubic');
       window = 10;
       dat=filter(ones(1,window)/window,1,dat);
       dat=[yday dat];
       dat=dat(6:359,:);
       maxchl=max(dat(:,2));
       minchl=min(dat(:,2));
       dmax = dat(dat(:,2) == maxchl,:);
       dmax = dmax(1,1);
       dmin = dat(dat(:,2) == minchl,:);
       dmin = dmin(1,1);   
       out(rowr,1)=dmin;
       out(rowr,2)=dmax;
       end;
   
   end;

   latout = [89.5:-1:-89.5]';
   latout = repmat(latout,360,1);
   %latout = latout*(-1);
   lonout = [-179.5:1:179.5];
   lonout = repmat(lonout,180,1);
   lonout = reshape(lonout,64800,1);
        
    out(:,4) = lonout;
    out(:,5) = latout;
    out(:,3)=year;
    
    out(any(isnan(out),2),:) = [];
    
    eval([ 'a' num2str(year) ' = out;' ]);

           
end;

    eval(['save ', 'coords_1_16', ' coords -ascii -tabs'])
a = [a1997; a1998; a1999; a2000; a2001; a2002; a2003];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



stem(dat) 

eval(['save ', 'coords_1_16', ' coords -ascii -tabs'])

%remove rows with missing data


swfs_chla_mo_1_4deg = [n1997; n1998; n1999; n2000; n2001; n2002; n2003; n2004; n2005; n2006; n2007; n2008; n2009];

latout = [89.875:-.25:-89.875]';
latout = repmat(latout,1440,1);
latout = latout*(-1);
lonout = [-179.875:.25:179.875];
lonout = repmat(lonout,720,1);
lonout = reshape(lonout,1036800,1);
coords=[lonout latout];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extracts seawifs 8day chl data to 5 degree resolution
cd('I:\Data\phytoplankton\remote.sensing\ocean_toolbox\seawifs_chl\8_day')
cd('E:\Data\phytoplankton\remote.sensing\ocean_toolbox\seawifs_chl\8_day')

outp=[];

for n=1:702;
    
        eval(['file = ''swfs_8day_opendap_ (', num2str(n),').mat''';]);
       
load(file);
chl=D.Chlor;
time=D.user_friendly_time;
year = str2num(time(1:4));
month = str2num(time(6:7));
day = str2num(time(9:10));
dayofyear = ((month*30.33)+day)-30.33;


display(time);
%sst(sst==-3)=NaN;
chl = resizem(chl, [36,72]);
%sst = flipud(sst);
chl=reshape(chl,2592,1);

latout = [87.5:-5:-87.5]';
latout = repmat(latout,72,1);

lonout = [-177.5:5:177.5];
lonout = repmat(lonout,36,1);
lonout = reshape(lonout,2592,1);

out = [lonout latout chl];
out(:,4) = year;
out(:,5) = month;
out(:,6) = day;
out(:,7) = dayofyear;

out(any(isnan(out),2),:) = [];
outp = vertcat(outp, out);

end;

eval(['save ', 'swfs_5deg_8day', ' outp -ascii -tabs'])

csvwrite('swfs_chla_8day_5deg.csv',outp);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extracts czcs 8day chl data to 5 degree resolution
cd('I:\Data\phytoplankton\remote.sensing\ocean_toolbox\czcs_chl\8day\mat_files')
cd('E:\Data\phytoplankton\remote.sensing\ocean_toolbox\czcs_chl\8day\mat_files')

n=1;

outp=[];

for n=1:341;
    
        eval(['file = ''czcs_8day (', num2str(n),').mat''';]);
       
load(file);
chl=D.Chlor;
time=D.user_friendly_time;
year = str2num(time(1:4));
month = str2num(time(6:7));
day = str2num(time(9:10));
dayofyear = ((month*30.33)+day)-30.33;

%lon=D.longitude;
%lat=D.latitude;

display(time);
%sst(sst==-3)=NaN;
chl = resizem(chl, [36,72]);
%sst = flipud(sst);
chl=reshape(chl,2592,1);

latout = [87.5:-5:-87.5]';
latout = repmat(latout,72,1);

lonout = [-177.5:5:177.5];
lonout = repmat(lonout,36,1);
lonout = reshape(lonout,2592,1);

out = [lonout latout chl];
out(:,4) = year;
out(:,5) = month;
out(:,6) = day;
out(:,7) = dayofyear;

out(any(isnan(out),2),:) = [];
outp = vertcat(outp, out);

end;

eval(['save ', 'czcs_5deg_8day', ' outp -ascii -tabs'])







%%extracts 8day czcs data to 1 degree resoluton
outp=[];

for n=1:341;
    
        eval(['file = ''czcs_8day (', num2str(n),').mat''';]);
       
load(file);
chl=D.Chlor;
time=D.user_friendly_time;
year = str2num(time(1:4));
month = str2num(time(6:7));
day = str2num(time(9:10));
dayofyear = ((month*30.33)+day)-30.33;

chl = resizem(chl, [180,360]);
chl=reshape(chl,64800,1);

%set up geographic coordinates
latout = [89.5:-1:-89.5]';
latout = repmat(latout,360,1);
lonout = [-179.5:1:179.5];
lonout = repmat(lonout,180,1);
lonout = reshape(lonout,64800,1);

%lon=D.longitude;
%lat=D.latitude;

display(time);

out = [lonout latout chl];
out(:,4) = year;
out(:,5) = month;
out(:,6) = day;
out(:,7) = dayofyear;

out(any(isnan(out),2),:) = [];
outp = vertcat(outp, out);

end;

eval(['save ', 'swfs_5deg_8day', ' outp -ascii -tabs'])

csvwrite('czcs_chla_8day_5deg.csv',outp);






%%extracts 8day seawifs data to 1 degree resoluton
cd('C:\Users\sailfish\Data\phytoplankton\remote.sensing\seawifs_chl_nc\8_day');
outp=[];

for n=1:662;
    
        eval(['file = ''swfs_8day_opendap_ (', num2str(n),').mat''';]);
       
load(file);
chl=D.Chlor;
time=D.user_friendly_time;
year = str2num(time(1:4));
month = str2num(time(6:7));
day = str2num(time(9:10));
dayofyear = ((month*30.33)+day)-30.33;

chl = resizem(chl, [180,360]);
chl=reshape(chl,64800,1);

%set up geographic coordinates
latout = [89.5:-1:-89.5]';
latout = repmat(latout,360,1);
lonout = [-179.5:1:179.5];
lonout = repmat(lonout,180,1);
lonout = reshape(lonout,64800,1);

%lon=D.longitude;
%lat=D.latitude;

display(time);

out = [lonout latout chl];
out(:,4) = year;
out(:,5) = month;
out(:,6) = day;
out(:,7) = dayofyear;

out(any(isnan(out),2),:) = [];
outp = vertcat(outp, out);

end;

eval(['save ', 'swfs_5deg_8day', ' outp -ascii -tabs'])

csvwrite('swfs_chla_8day_1deg.csv',outp);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reads in hermes GSM chl-a data and determines seasonal min/max
cd('I:\Data\phytoplankton\remote.sensing\hermes\chl_gsm\daily_100km_case1')

for year=1997:2011;
    
        for month=1:12;
                       
            if ((month==1)); daya=1:31;
            elseif ((month==2)); daya=1:28;
            elseif ((month==3 || month==5 || month==7 || month==8 || month==10 || month==12)); daya=1:31;
            else daya=1:30;
            end;
           
                    for day=min(daya):max(daya);
                        
         if (year<=2001) && (month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]);
        elseif (year<=2001) && (month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]);
        elseif (year<=2001) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]); 
        elseif (year<=2001) && (month>=10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]); 
                
        elseif (year==2002) && (month<=3) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]);
        elseif (year==2002) && (month<=3) && (day>=10)            
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-SWF_CHL1_DAY_00.nc''';]);
               
        %read in swf/mer data
        elseif (year==2002) && (month>=4 && month<=5) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERSWF_CHL1_DAY_00.nc''';]); 
        elseif (year==2002) && (month>=4 && month<=5) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-MERSWF_CHL1_DAY_00.nc''';]); 
                         
        %read in swf/mer/mod data
        elseif (year==2002) && (month>=6 && month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]); 
        elseif (year==2002) && (month>=6 && month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]); 
        
         elseif (year==2002) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year==2002) && (month>=10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
        
        elseif (year>2003) && (month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year>2003) && (month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year>2003) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);
         elseif (year>2003) && (month>=10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_GSM-MERMODSWF_CHL1_DAY_00.nc''';]);        
        
        elseif (year>2010) && (month<10) && (day<10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMOD_CHL1_DAY_00.nc''';]);
         elseif (year>2010) && (month<10) && (day>=10)
        eval(['file = ''L3m_', num2str(year),'0',num2str(month),num2str(day),'__GLOB_100_GSM-MERMOD_CHL1_DAY_00.nc''';]);
         elseif (year>2010) && (month>=10) && (day<10)
        eval(['file = ''L3m_', num2str(year),num2str(month),'0',num2str(day),'__GLOB_100_GSM-MERMOD_CHL1_DAY_00.nc''';]);
        else
        eval(['file = ''L3m_', num2str(year),num2str(month),num2str(day),'__GLOB_100_GSM-MERMOD_CHL1_DAY_00.nc''';]);        
        end;
        
          %try to extract chl but inserts NaN if file missing
          try
          c=netcdf.open(file,'NC_WRITE');
          chl = getnc(file,'CHL1_mean');
          netcdf.close(c);
          catch ME
          display('missing');
          chl = zeros(180,360);
          chl(:,:) = NaN;
          end;

chl = resizem(chl, [36,72]);
chl=reshape(chl,2592,1);
latout = [87.5:-5:-87.5]';
latout = repmat(latout,72,1);
lonout = [-177.5:5:177.5];
lonout = repmat(lonout,36,1);
lonout = reshape(lonout,2592,1);
dayofyear = ((month*30.33)+day)-30.33;
       
out = [lonout latout chl];
out(:,4) = year;
out(:,5) = month;
out(:,6) = day;
out(:,7) = dayofyear;
out(any(isnan(out),2),:) = [];
%outp = vertcat(outp, out);
eval([ 'a' num2str(day) ' = out;' ]);

       end;
  
  if (month==2);  aa=[a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; a13; a14; a15; a16; a17; a18; a19; a20; a21; a22; a23; a24; a25; a26; a27; a28];
  elseif (month==4 || month==6 || month==9 || month==11); aa=[a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; a13; a14; a15; a16; a17; a18; a19; a20; a21; a22; a23; a24; a25; a26; a27; a28; a29; a30];     
  else aa=[a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; a13; a14; a15; a16; a17; a18; a19; a20; a21; a22; a23; a24; a25; a26; a27; a28; a29; a30; a31];
  end;
  
  eval([ 'm' num2str(month) ' = aa;' ]);
 
        end;
        
        aaa=[m1; m2; m3; m4; m5; m6; m7; m8; m9; m10; m11; m12];
        eval([ 'a' num2str(year) ' = aaa;' ]);
           
end;

    eval(['save ', 'coords_1_16', ' coords -ascii -tabs'])
a = [a1997; a1998; a1999; a2000; a2001; a2002; a2003; a2004; a2005; a2006; a2007; a2008; a2009; a2010; a2011];
csvwrite('gsm_chla_1day_5deg.csv',a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








