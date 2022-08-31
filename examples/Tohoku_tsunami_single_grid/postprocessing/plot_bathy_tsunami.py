### PLOT Tohoku_tsunami bathymetry ###

# import necessary modules
import numpy as np               
import matplotlib.pyplot as plt


# write your OWN PC folder path for fdir and dep.
# Remember that we use for Mac & Linux machines '/', while on windows '\'
dep=np.loadtxt('/Users/Gaby/Desktop/Postprocessing-Workshop/simple_cases_output/Tohoku_tsunami/depth_30min.txt')

# define bathy location
n,m = np.shape(dep)
dx = 0.5
dy = 0.5

x = np.asarray([float(xa)*dx+132.01667 for xa in range(m)])
y = np.asarray([float(ya)*dy-59.98333 for ya in range(n)])

# figure size option 
wid=6    # width
length=4 # length

# Plot figure
fig = plt.figure(figsize=(wid,length),dpi=200)
ax = fig.add_subplot(1,1,1)
fig.subplots_adjust(hspace=1,wspace=.25)

plt.pcolor(x, y, -1*dep,cmap='terrain')
plt.axis('tight')  
plt.ylabel('Lat (deg)')
plt.xlabel('Lon (deg)')

# figure colorbar
cbar=plt.colorbar()
cbar.set_label('Bathymetry (m)', rotation=90)

# save figure
fig.savefig('tsunami_bathy.png', dpi=fig.dpi)
