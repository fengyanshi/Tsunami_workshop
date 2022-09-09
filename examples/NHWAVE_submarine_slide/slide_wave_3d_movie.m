clear all

fdir1='output/';
dep=load('depth.txt');

[ny,nx]=size(dep);
dx=20;
dy=20;
view_scale=4;
depth_threshold = 1.0;

x1=(0:1:nx-1)*dx;
y1=(0:1:ny-1)*dy;

[x,y]=meshgrid(x1,y1);

nfiles=[10:1:50];

set(gcf,'units','inches','paperunits','inches','papersize', [8 4],'position',[1 1 8 4],'paperposition',[0 0 8 4]);

myVideo = VideoWriter('videoOut.mp4','MPEG-4');
myVideo.FrameRate = 10;  
myVideo.Quality = 100;
vidHeight = 576; %this is the value in which it should reproduce
vidWidth = 864; %this is the value in which it should reproduce
open(myVideo);

for num=1:length(nfiles)

fnum=sprintf('%.4d',nfiles(num));

eta=load([fdir1 'eta_' fnum]);
slide=load([fdir1 'slide_' fnum]);
eta(eta+dep-slide<depth_threshold)=NaN;

clf

surface(x,y,-dep, 'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','CDataMapping','direct','DiffuseStrength',0.5)

hold on

slide_bottom=-dep+slide;
slide_bottom(slide<depth_threshold)=NaN;
surface(x,y,slide_bottom, 'FaceColor',[0.54 0.44 0.35],'EdgeColor','none','CDataMapping','direct','DiffuseStrength',0.5)

hsurf=surface(x,y,eta*view_scale,'FaceColor',[0.1 0.6 1],'EdgeColor','none','CDataMapping','direct','DiffuseStrength',0.5);

alpha(hsurf,0.75);
lightangle(-35,70)
view([47 38])
grid

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

pause(0.1)
F = print('-RGBImage','-r600');
J = imresize(F,[vidHeight vidWidth]);
mov(num).cdata = J;


writeVideo(myVideo,mov(num).cdata);

end

close(myVideo)