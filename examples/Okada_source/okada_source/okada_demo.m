% GRIDS
clear all

file_name='depth.asc';

fid = fopen(file_name,'r');
 numLines = 6;
 text_line = cell(numLines,1);
 for ii = 1:numLines
     text_line(ii) = {fgetl(fid)}; 
 end
 fclose(fid);
m=str2double(regexp(text_line{1},'\d*','Match'));
n=str2double(regexp(text_line{2},'\d*','Match'));
lon1=str2double(regexp(text_line{3},'\d+(\.)?(\d+)?','Match'));
lat1=str2double(regexp(text_line{4},'\d+(\.)?(\d+)?','Match'));
dx0=str2double(regexp(text_line{5},'\d+(\.)?(\d+)?','Match'));

Dep0=dlmread('depth.asc','',6,0);
Dep0=-flipud(Dep0); % for slope calculation

lon0=[0:m-1]*dx0+lon1;
lat0=[0:n-1]*dx0+lat1;

[Lon0,Lat0] = meshgrid(lon0,lat0);

% slopes 
R_earth = 6371e3;
[Dep0x,Dep0y] = gradient(Dep0,dx0);
Dep0x = Dep0x ./ (R_earth/360*2*pi*sin(Lat0*pi/180));
Dep0y = Dep0y / (R_earth/360*2*pi);

% SOURCES

source_data=dlmread('aleutians_source.txt','',1,0);

fault_lon0 =source_data(:,1);
fault_lat0 = source_data(:,2);
fault_length = source_data(:,3)*1000;
fault_width = source_data(:,4)*1000;
fault_dip = source_data(:,5);
fault_rake = source_data(:,6);
fault_strike = source_data(:,7);
fault_slip = source_data(:,8); 
fault_depth0 = source_data(:,9)*1000;

uE = 0*Lat0;
uN = 0*Lat0;
uZ = 0*Lat0;
for i=1:1%length(fault_lat0)
  [dE,dN] = ll2tm(Lat0,Lon0,fault_lat0(i),fault_lon0(i));

  [uEi,uNi,uZi] = okada85(dE,dN,...
    fault_depth0(i),...
    fault_strike(i),fault_dip(i),...
    fault_length(i),fault_width(i),fault_rake(i),fault_slip(i),0.0);
  uE = uE + uEi;
  uN = uN + uNi;
  uZ = uZ + uZi;
end
Z0 = uZ + uE.*Dep0x + uN.*Dep0y;
U0 = Z0*0;

figure(1)
clf
colormap jet
pcolor(Lon0,Lat0,Z0),shading flat
hold on
contour(Lon0,Lat0,Dep0,[-100:10:0],'w-')
caxis([-5 5])
print -djpeg100 aleutians.jpg

save('aleutians.txt','Z0','-ascii','-double','-tabs');

