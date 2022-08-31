clear all
%fdir1='/Volumes/Seagate Backup Plus Drive/Babak/CresentCity/2arcmin/';
fdir1='/Volumes/Seagate Backup Plus Drive/Babak/CresentCity/4arcsec/';
%fdir3='/Volumes/Seagate Backup Plus Drive/Babak/CresentCity/16arcsec/';

%dep=load('../external_files/depth_30min.txt');
m=1232;
n=1000;

%[n,m]=size(dep);
dx=4/3600.;;
dy=4/3600.;;
x=[0:m-1]*dx+234.75;
y=[0:n-1]*dy+41.4;



figure(3)
wid=7;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

clf


    
%fnum=sprintf('%.5d',nfile(num));
eta=load([fdir1 'eta_0180']);

eta(eta>2)=NaN;

%subplot(length(nfile),1, num)

pcolor(x,y,eta),shading flat
hold on
caxis([-0.5 0.5])
%title([' Time = ' hr{num} ' hr '])

%caxis([-0.05 0.05])

ylabel(' Lat (deg) ')
xlabel(' Lon (deg) ')
cbar=colorbar;
set(get(cbar,'ylabel'),'String','\eta (m) ')
set(gcf,'Renderer','zbuffer');

%print -djpeg eta_inlet_shoal_irr.jpg