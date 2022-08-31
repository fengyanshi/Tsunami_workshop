clear all;

%==============
dx(1)=10;     dy(1)=10;
Nghost=4;
fdir='output/';
%==============


% base grid -------------------------
grid=load('subgrid_info.txt');
dep=load([fdir 'Grd01_dep.out']);
num_grid = size(grid,1)+1;
[n(1), m(1)]=size(dep);
x{1} = 0 : dx(1) : (m(1)-1)*dx(1);
y{1} = 0 : dy(1) : (n(1)-1)*dy(1);
[Xgrid{1}, Ygrid{1}]=meshgrid(x{1},y{1});
%-----------------------------------------

% subgrid ---------------------------
for k=2:num_grid
m(k)=grid(k-1,1);
n(k)=grid(k-1,2);
isk(k)=grid(k-1,3);
trackx(k)=grid(k-1,4);
tracky(k)=grid(k-1,5);

dx(k)=dx(k-1)/isk(k);
dy(k)=dy(k-1)/isk(k);

x{k} = 0 : dx(k) : (m(k)-1)*dx(k);
y{k} = 0 : dy(k) : (n(k)-1)*dy(k);

x{k} = x{k} + x{k-1}(1)+(trackx(k)-Nghost-1)*dx(k-1) + Nghost*dx(k);
y{k} = y{k} + y{k-1}(1)+(tracky(k)-Nghost-1)*dy(k-1) + Nghost*dy(k);

[Xgrid{k}, Ygrid{k}]=meshgrid(x{k},y{k});
end
%----------------------------------------

time=load([fdir 'Grd01_track.txt']);
nstart=input('nstart=');
nend=input('nend=');

icount=0;
for num=nstart:1:nend
icount=icount+1;

fname=sprintf('%.5d',num);

myVideo = VideoWriter('videoOut.mp4','MPEG-4');
myVideo.FrameRate = 25;  
myVideo.Quality = 100;
vidHeight = 576; %this is the value in which it should reproduce
vidWidth = 576; %this is the value in which it should reproduce
open(myVideo);

clf
colormap jet;

for k=1:num_grid
fgrid=sprintf('%.2d',k);
eta=load([fdir 'Grd' fgrid '_eta_' fname]);
pcolor(Xgrid{k},Ygrid{k},eta),shading interp

hold on

colorbar

plot( [Xgrid{k}(1,1) Xgrid{k}(1,end) Xgrid{k}(end,end) Xgrid{k}(end,1) Xgrid{k}(1,1)], ...
    [Ygrid{k}(1,1) Ygrid{k}(1,end) Ygrid{k}(end,end) Ygrid{k}(end,1) Ygrid{k}(1,1)], ...
    'k--', 'LineWidth', 2 )

end

title (['Time = ',num2str(time(num), '%5.2f'), ' (s)'],'fontsize',11);
box on;
caxis([-0.5 0.5])

pause(0.5)

F = print('-RGBImage','-r150');
J = imresize(F,[vidHeight vidWidth]);
mov(k).cdata = J;


writeVideo(myVideo,mov(k).cdata);

end

close(myVideo)

