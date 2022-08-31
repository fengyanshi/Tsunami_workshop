clear all
dep=load('../external_files/depth_30min.txt');

[n,m]=size(dep);
dx=0.5;
dy=0.5;
x=[0:m-1]*dx+132.01667;
y=[0:n-1]*dy-59.98333;


wid=6;
len=4;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf

pcolor(x,y,-dep),shading interp
demcmap(dep)

ylabel(' Lat (deg) ')
xlabel(' Lon (deg) ')
%cbar=colorbar;
%set(get(cbar,'ylabel'),'String','-depth (m) ')
set(gcf,'Renderer','zbuffer');
