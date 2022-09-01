clear all
fdir='output/';

m=1000;
dx=1.0;
SLP=0.0125;
Xslp = 500.0;

% bathy
x=[0:m-1]*dx;
dep=zeros(m)+5.0;
dep(x>Xslp)=5.0-(x(x>Xslp)-Xslp)*SLP;

% wavemaker and sponge

files=[0,50,100];

wid=6;
len=8;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf

for num=1:length(files)
subplot(length(files),1,num)
fnum=sprintf('%.5d',files(num));
eta=load([fdir 'eta_' fnum]);
plot(x,-dep,'k',x,eta,'b','LineWidth',2)
hold on
axis([0 1000 -5 1.0])
grid
xlabel('x(m)')
ylabel('eta(m)')
title([' Time = ' num2str(files(num)*1) ' sec '])
end
print -djpeg100 snapshot.jpg
