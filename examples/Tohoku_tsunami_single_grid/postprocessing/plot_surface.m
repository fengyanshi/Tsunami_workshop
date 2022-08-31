% Clearing Matlab workspace
clear all

% -----------------------
% ----- User Input ------
% -----------------------

% Directory of output data files
fdir='../work/output/';

% Directory of depth file
depth_file_dir = '../external_files/';

% Time series files to plot
nfile=[2 9];

% Time values for series in hours
time={'1' '8'};

% -----------------------
% -- End of user input --
% -----------------------


% Getting depth file and determining domain dimensions
dep=load( [ depth_file_dir 'depth_30min.txt' ]);
[n,m]=size(dep);

% Setting up partition
dx=0.5;
dy=0.5;
x=[0:m-1]*dx+132.01667;
y=[0:n-1]*dy-59.98333;


% Dimensions of plot window 
wid=5;
len=7;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf
colormap jet

for num=1:length(nfile)

% Padding integer values with zeros
% to be 5 letters long e.g. 1 -> 00001
fnum=sprintf('%.5d',nfile(num));

% Loading data from files
eta=load([fdir 'eta_' fnum]);

% Removing masked regions from plot
eta(dep<0)=NaN;

% Plotting data on different subplot regions
subplot(length(nfile),1, num)

% Plotting wave displacement
pcolor(x,y,eta),shading flat

title([' Time = ' time{num} ' hr '])

% Reducing contour plot range magnitude
% for subsequent subplots
if num==1
caxis([-0.5 0.5])
else
caxis([-0.05 0.05])
end

ylabel(' Lat (deg) ')
xlabel(' Lon (deg) ')

cbar=colorbar;
set(get(cbar,'ylabel'),'String','\eta (m) ')
set(gcf,'Renderer','zbuffer');

end
