% Clearing Matlab workspace
clear all

% -----------------------
% ----- User Input ------
% -----------------------

% Directory of depth file
fdir = '../external_files/'

% Change to plot color bars
plot_color_bars = true;

% -----------------------
% -- End of user input --
% -----------------------


% Getting depth file and determining domain dimensions
dep=load( [fdir 'depth_30min.txt']);
[n,m]=size(dep);

% Setting up partition
dx=0.5;
dy=0.5;
x=[0:m-1]*dx+132.01667;
y=[0:n-1]*dy-59.98333;

% Plot window dimensions
wid=6;
len=4;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf

% Plotting depth
pcolor(x,y,-dep),shading interp

% Using demcmap if Mapping Toolbox is installed
has_demcmap = ~isempty(which('demcmap'));
if (has_demcmap)
demcmap(dep)
end

ylabel(' Lat (deg) ')
xlabel(' Lon (deg) ')

% Ploting color bars if toggled
if ( plot_color_bars )
cbar=colorbar;
set(get(cbar,'ylabel'),'String','-depth (m) ') 
end

set(gcf,'Renderer','zbuffer');
