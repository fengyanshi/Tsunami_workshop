import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

# Setup Directory Paths
simulationDir = '..'+'/'+'work_solitary'
outputDir = simulationDir+'/'+'output'
plotsDir = simulationDir+'/'+'plots'

# check if plots directory exists if not create it
if not os.path.exists(plotsDir):
    os.makedirs(plotsDir)

# Read in the depth bathymetry file and setup dimensions
depth = (np.loadtxt('../work_solitary/depth_levee.txt'))*-1  # Change Depth to (-) under MWL and (+) over MWL.
[Nglob,Mglob]  = depth.shape   # recall that Fotran/FUNWAVE and C/C++/Matlab/Python are reversed! 

# Setting up field dimensions and time loops
numOfETA = 15 # number of ETA text files
Lt = 40.0 # total horizontal lenght (m)
WaveMaker = 4.0 # wavemaker x location (m)
plotInt = 30.0 # plot interval (sec)

def readETAData(num,numOfETA,postprocessDir):
    """Function that reads in the eta_00### ASCII file and returns
    it in the 2D NumPy array of size [Nglob,Mglob]."""

    assert num in range(1,numOfETA+1), "File Index:%d is not in range of station numbers." %(num,)
 
    #fileIndex = str(num)
    fileIndex = num
    fileName = outputDir+'/'+'eta_{0:05d}'.format(fileIndex)
    print("READING IN: ", fileName)
    freeSurface = np.loadtxt(fileName)

    return freeSurface

# Plot and save every ETA file:
for i in range(1, numOfETA+1):
    surf = readETAData(i,numOfETA,outputDir)        # Y axis = eta data
    x = np.linspace(0, Lt, Mglob)                     # X axis = Total Length (m)
    
    fig  = plt.figure(figsize=(18,4), dpi=600)
    ax = fig.add_subplot(1,1,1)
    plt.plot(x, surf[1,:], 'c-', linewidth = 0.2)
    plt.axis([-0.5,Lt,min(depth[1,:])-.05,1])

    plt.xlabel('Length (%4.2f m)' % (Lt), fontsize = 12, fontweight = 'bold')
    plt.ylabel('Height (m)', fontsize = 12, fontweight = 'bold')
    

    # Water Fill:
    plt.fill_between(x, depth[1,:], surf[1,:],
                     where = surf[1,:] > depth[1,:],
                     facecolor = 'cyan', interpolate =True)
    
    # Bottom Fill:
    plt.fill_between(x, min(depth[1,:])-.05, depth[1,:], 
                               where= depth[1,:] > (depth[1,:]-.05),facecolor = '0.35',
                               hatch = 'X')

    #Annotations:
    an = plt.annotate('\n~\n|', xy=(WaveMaker,.05), xycoords='data',fontsize = 20,
                                 ha='center', va='center')
    Time = plt.annotate('Time: %4.2f sec'%(i*plotInt-plotInt), xy=(Lt*.85,3),fontsize = 16,
                                 ha='center', va='center')

    fileIndex = i
    fileName = plotsDir+'/'+'surf_{0:05d}.png'.format(fileIndex)
    plt.savefig(fileName, ext="png", bbox_inches='tight')
    plt.close()

