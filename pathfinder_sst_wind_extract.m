cd('C:\Users\sailfish\Downloads\pathfinder_SST\v5.2v2')

%file='pathfinder_sst_daily (6162).nc'
file='pathfinder_sst_daily (6163).nc'
tm=ncreadatt(file,'/','start_time');
                yeart=tm(1:4);
                
ncdisp(file)
sst = getnc(file,'sea_surface_temperature');
wnd = getnc(file,'wind_speed');
tm=ncreadatt(file,'/','start_time');
year=tm(1:4);
month=tm(5:6);
day=tm(7:8);

time = getnc(file,'time');

id=746; 
day=300;
   
lon2 = getnc(file,'lon')';
lat2 = getnc(file,'lat');
latlen=length(lat2);
lonlen=length(lon2);
latout = repmat(lat2,lonlen,1);
lonout = repmat(lon2,latlen,1);
lonout = reshape(lonout,length(latout),1);
lnth=length(lonout);

outp=[NaN NaN NaN NaN NaN NaN NaN];

tic
for id=1:11360;
                
                eval(['file = ''pathfinder_sst_daily (',num2str(id),').nc''';]);          
               
                try
                tm=ncreadatt(file,'/','start_time');
                yeart=tm(1:4);
                montht=tm(5:6);
                dayt=tm(7:8);
                sst = getnc(file,'sea_surface_temperature');
                sst(sst<=-30000) = NaN;
                sst=sst-273.15;
                wnd = getnc(file,'wind_speed');
                wnd(wnd<0) = NaN;

                %sst = resizem(sst, [180 360]);%resize to 5 degree
                sst = reshape(sst, lnth,1);
                sst=reshape(sst,length(lonout),1);
                %wnd = resizem(wnd, [180 360]);%resize to 5 degree
                wnd = reshape(wnd, lnth,1);
                wnd=reshape(wnd,length(lonout),1);
                
                out=[lonout latout sst wnd];
                out(:,5)=str2num(yeart);
                out(:,6)=str2num(montht);
                out(:,7)=str2num(dayt);
                out=out(out(:,2)> -67 & out(:,1)< -65.5 & out(:,2)< 45 & out(:,2)> 42.75,:);    
                outp = vertcat(outp, out);
                clear out sst wnd 
                outp(any(isnan(outp),2),:) = [];

                catch ME
                display('missing');
                %a = zeros(4320, 8640);
                %a(:,:) = NaN;
                end;
      
%outp(any(isnan(outp),2),:) = [];
%eval([ 'a' num2str(year) ' = outp;' ]);
%clear outp out a;
    
end;

toc





mdata = [a2005; a2006];%2008 contains missing data so not included

scatter(mdata(:,1),mdata(:,2),'.')

%cd('N:\data\chl_phenology\data')
csvwrite('pathfinder_sst_daily_2005_2006.csv',mdata);
max(mdata(:,3)

