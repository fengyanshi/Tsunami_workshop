### PLOT Tohoku_tsunami WAVE ###

# import necessary modules
import numpy as np               
import matplotlib.pyplot as plt
import matplotlib as mat


# write your OWN PC folder path for fdir and dep
# Remember that we use for Mac & Linux machines '/', while on windows '\'
dep=np.loadtxt('/Users/Gaby/Desktop/Postprocessing-Workshop/simple_cases_output/Tohoku_tsunami/depth_30min.txt')

fdir = '/Users/Gaby/Desktop/Postprocessing-Workshop/simple_cases_output/Tohoku_tsunami/tohoku_tsunami/'

# define bathy location
n,m = np.shape(dep)
dx = 0.5
dy = 0.5

x = np.asarray([float(xa)*dx+132.01667 for xa in range(m)])
y = np.asarray([float(ya)*dy-59.98333 for ya in range(n)])

nfile = [2, 9]    # range of eta files you want to plot
hr = ['1','8']  # time  you want to plot

# figure size option 
wid=5    # width
length=7 # length

# Plot figure
fig = plt.figure(figsize=(wid,length),dpi=200)

for num in range(len(nfile)):
    fnum= '%.5d' % nfile[num]
    eta = np.loadtxt(fdir+'eta_'+fnum)
    eta_masked = np.ma.masked_where(dep<0,eta) # do nt plot where dep<0

    ax = fig.add_subplot(len(nfile),1,num+1)
    fig.subplots_adjust(hspace=.45)
    plt.pcolor(x, y, eta_masked,cmap='jet')

    title = 'Time = '+hr[num]+ ' hr'
    plt.title(title)
    plt.axis('tight')

    plt.ylabel('Lat (deg)')
    plt.xlabel('Lon (deg)')
    cbar=plt.colorbar()
    cbar.set_label(r'$\eta$'+' (m)', rotation=90)

    if num == 0:
        plt.clim(-0.5, 0.5)
    else:
        plt.clim(-0.05, 0.05)

    mat.rcParams.update({'font.size': 10})
# save figure  
fig.savefig('eta_tsunami.png', dpi=fig.dpi)
