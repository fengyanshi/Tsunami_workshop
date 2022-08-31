clear all
fdir='output/';

dep=load('../okada_source/model_depth.txt');

[n,m]=size(dep);
dx=0.03333;
dy=0.03333;
x=[0:m-1]*dx+132.01667;
y=[0:n-1]*dy-59.98333;


nfile=[1 15];

figure(1)
clf
wid=5;
len=7;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
colormap jet


for num=1:length(nfile)
    
fnum=sprintf('%.5d',nfile(num));
eta=load([fdir 'eta_' fnum]);

eta(dep<0)=NaN;

subplot(length(nfile),1, num)

pcolor(x,y,eta),shading flat
caxis([-1 1])
title([' Time = ' num2str((nfile(num)-1)*2) ' min '])

ylabel(' Lat (deg) ')
xlabel(' Lon (deg) ')
cbar=colorbar;
set(get(cbar,'ylabel'),'String','\eta (m) ')
set(gcf,'Renderer','zbuffer');

end
print -djpeg eta.jpg