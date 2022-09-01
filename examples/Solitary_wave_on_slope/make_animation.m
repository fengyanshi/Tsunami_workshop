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

files=[1:2:199];

myVideo = VideoWriter('videoOut.mp4','MPEG-4');
myVideo.FrameRate = 10;  
myVideo.Quality = 100;
vidHeight = 400; %this is the value in which it should reproduce
vidWidth = 800; %this is the value in which it should reproduce
open(myVideo);

wid=8;
len=4;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

for num=1:length(files)
fnum=sprintf('%.5d',files(num));
eta=load([fdir 'eta_' fnum]);
clf
plot(x,-dep,'k',x,eta,'b','LineWidth',2)
hold on
axis([0 1000 -5 1.0])
grid
xlabel('x(m)')
ylabel('eta(m)')
title([' Time = ' num2str(files(num)*1) ' sec '])
pause(0.1)
F = print('-RGBImage','-r150');
J = imresize(F,[vidHeight vidWidth]);
mov(num).cdata = J;
writeVideo(myVideo,mov(num).cdata);

end
close(myVideo)

