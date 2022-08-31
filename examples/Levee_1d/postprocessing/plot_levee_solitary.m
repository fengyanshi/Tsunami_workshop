% Clearing Matlab workspace
clear all

% -----------------------
% ----- User Input ------
% -----------------------

% Directory of output data files
fdir='../input_files/output/';

% Directory for plot files
plotDir = [fdir,'plots/'];

% Time series files to plot
nfile=1:301;

% Change to plot color bars
plot_color_bars = true;

% -----------------------
% -- End of user input --
% -----------------------


% Getting domain dimensions from depth file
dep=-load([fdir 'dep.out']);
[m,n]=size(dep);
% Removing artifical y-axis points for 2D simulation
dep = dep(2,:);

% Setting up partition
dx = 0.08;
x = (0:n-1)*dx;

% Location of wavemaker
x_wavemaker=[4,4];
y_wavemaker=[-1,1];

% Plot window dimensions
wid=11;
len=2.5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf

% Making plot directory if it doesn't exsist
if ~isdir(plotDir)
    mkdir( plotDir )
end

for i = 1:301

    num=nfile(i);
    % Padding integer values with zeros
    % to be 5 letters long e.g. 1 -> 00001
    fnum=sprintf('%.5d',nfile(num));

    % Loading data from files
    fileDir = [fdir 'eta_' fnum];
    eta=load(fileDir);
    
    % Removing artifical y-axis points for 2D simulation
    eta = eta(2,:);
    
    fprintf( ['READING IN: eta_' fnum , '\n'])

    % Plotting depth
    fill([x fliplr(x)] , [dep (dep*0 -0.5)] ,[0.8 0.8 0.8])
    
    % Plotting wave displacement and wave maker location
    hold on
        fill([x fliplr(x)], [eta fliplr(dep) ] , 'c' )
        plot(x_wavemaker,y_wavemaker,'--r','LineWidth',2)
    hold off

    h2=text(3,.15,'Wavemaker','Color','r');
    set(h2, 'rotation', 90)

    xlabel('x (m)')
    ylabel('Height (m)' )

    ylim([-0.45,1])

    % Printing plot to png file
    pngFile = [plotDir 'surf_' fnum '.png'];
    print('-dpng', pngFile)

end





