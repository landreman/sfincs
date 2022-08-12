#!/usr/bin/env python

# This python script plots the distribution functions from a single SFINCS run.
# The script checks to see which dimensions have the highest resolution, and will plot contours of f versus those 2 coordinates.
# In the other coordinates, the first index is used.

ispecies = 0
# Which species to plot (0-based)

iteration = 0
# Which iteration to plot

import matplotlib.pyplot as plt
import h5py
import numpy
import inspect, os

#print ("This is "+ inspect.getfile(inspect.currentframe()))
this_filename = "sfincsPlotF"
print ("This is "+ this_filename)

try:
    f = h5py.File('sfincsOutput.h5','r')
except:
    print ("Unable to open sfincsOutput.h5")
    exit(1)

# The expression [()] converts from an h5py dataset to a numpy ndarray:

full_f_available = True
try:
    full_f = f['full_f'][()]
    print ("full_f is available.")
except KeyError:
    full_f_available = False
    print ("full_f is not available.")

delta_f_available = True
try:
    delta_f = f['delta_f'][()]
    print ("delta_f is available.")
except KeyError:
    delta_f_available = False
    print ("delta_f is not available.")

if (not full_f_available) and (not delta_f_available):
    print ("Stopping since neither full_f nor delta_f are available.")
    exit(0)

integerToRepresentTrue = int(f['integerToRepresentTrue'][()])

theta = f['export_f_theta'][()]
Ntheta = int(f['N_export_f_theta'][()])

zeta = f['export_f_zeta'][()]
Nzeta = int(f['N_export_f_zeta'][()])

x = f['export_f_x'][()]
Nx = int(f['N_export_f_x'][()])

export_f_xi_option = int(f['export_f_xi_option'][()])
if export_f_xi_option==0:
    Nxi = int(f['Nxi'][()])
    xi = range(0,Nxi)
    xi_label = 'Legendre mode number'
else:
    Nxi = int(f['N_export_f_xi'][()])
    xi = f['export_f_xi'][()]
    xi_label = 'xi'

numRepeatedDimensions = 0
if Ntheta>1:
    numRepeatedDimensions += 1
if Nzeta>1:
    numRepeatedDimensions += 1
if Nxi>1:
    numRepeatedDimensions += 1
if Nx>1:
    numRepeatedDimensions += 1

if numRepeatedDimensions < 2:
    print ("Error: This script is designed for plots of the distribution function versus 2 phase space coordinates, so you need to export f with >1 grid points in at least 2 of the coordinates.")
    exit(1)

print ("Resolution on which the distribution function is available:")
print ("N_export_f_theta: ",Ntheta)
print ("N_export_f_zeta:  ",Nzeta)
print ("N_export_f_xi:    ",Nxi)
print ("N_export_f_x:     ",Nx)
# Determine which 2 of the 4 coordinates have the highest resolution:
resolutions = sorted([Ntheta, Nzeta, Nxi, Nx])
print ("resolutions: ",resolutions)
plotThisCoordinate = [False,False,False,False]

temp = resolutions[3]
if temp==Ntheta:
    plotThisCoordinate[0]=True
elif temp==Nzeta:
    plotThisCoordinate[1]=True
elif temp==Nxi:
    plotThisCoordinate[2]=True
else:
    plotThisCoordinate[3]=True

temp = resolutions[2]
if temp==Ntheta and (not plotThisCoordinate[0]):
    plotThisCoordinate[0]=True
elif temp==Nzeta and (not plotThisCoordinate[1]):
    plotThisCoordinate[1]=True
elif temp==Nxi and (not plotThisCoordinate[2]):
    plotThisCoordinate[2]=True
else:
    plotThisCoordinate[3]=True

print ("plotThisCoordinate: ",plotThisCoordinate)

# Dimensions of delta_f and full_f are x, xi, zeta, theta, species, iteration

N=0
if plotThisCoordinate == [False,False,True,True]:
    # xi vs x
    N += 1
    abscissa = xi
    abscissa_label = xi_label
    ordinate = x
    ordinate_label = "x"
    if full_f_available:
        full_f_data  =  full_f[:,:,0,0,ispecies,iteration]
    if delta_f_available:
        delta_f_data = delta_f[:,:,0,0,ispecies,iteration]

if plotThisCoordinate == [False,True,False,True]:
    # zeta vs x
    N += 1
    abscissa = zeta
    abscissa_label = 'zeta'
    ordinate = x
    ordinate_label = "x"
    if full_f_available:
        full_f_data  =  full_f[:,0,:,0,ispecies,iteration]
    if delta_f_available:
        delta_f_data = delta_f[:,0,:,0,ispecies,iteration]

if plotThisCoordinate == [True,False,False,True]:
    # theta vs x
    N += 1
    abscissa = x
    abscissa_label = 'x'
    ordinate = theta
    ordinate_label = "theta"
    if full_f_available:
        full_f_data  =  full_f[:,0,0,:,ispecies,iteration].transpose()
    if delta_f_available:
        delta_f_data = delta_f[:,0,0,:,ispecies,iteration].transpose()

if plotThisCoordinate == [False,True,True,False]:
    # zeta vs xi
    N += 1
    abscissa = zeta
    abscissa_label = 'zeta'
    ordinate = xi
    ordinate_label = xi_label
    if full_f_available:
        full_f_data  =  full_f[0,:,:,0,ispecies,iteration]
    if delta_f_available:
        delta_f_data = delta_f[0,:,:,0,ispecies,iteration]

if plotThisCoordinate == [True,False,True,False]:
    # theta vs xi
    N += 1
    abscissa = xi
    abscissa_label = xi_label
    ordinate = theta
    ordinate_label = "theta"
    if full_f_available:
        full_f_data  =  full_f[0,:,0,:,ispecies,iteration].transpose()
    if delta_f_available:
        delta_f_data = delta_f[0,:,0,:,ispecies,iteration].transpose()

if plotThisCoordinate == [True,True,False,False]:
    # theta vs zeta
    N += 1
    abscissa = zeta
    abscissa_label = "zeta"
    ordinate = theta
    ordinate_label = "theta"
    if full_f_available:
        full_f_data  =  full_f[0,0,:,:,ispecies,iteration].transpose()
    if delta_f_available:
        delta_f_data = delta_f[0,0,:,:,ispecies,iteration].transpose()

if N != 1:
    print ("Error! Something went wrong. N = ",N)
    exit(1)


fig = plt.figure(1)
fig.patch.set_facecolor('white')

numRows = 1
numCols = 1
if full_f_available and delta_f_available:
        numCols = 3
plotNum = 1
numContours = 20

# Plot full f
if full_f_available:
    plt.subplot(numRows, numCols, plotNum)
    plt.contourf(abscissa,ordinate,full_f_data, numContours)
    plt.colorbar()
    plt.title('full f')
    plt.xlabel(abscissa_label)
    plt.ylabel(ordinate_label)
    plotNum += 1

# Plot delta f
if delta_f_available:
    plt.subplot(numRows, numCols, plotNum)
    plt.contourf(abscissa,ordinate,delta_f_data, numContours)
    plt.colorbar()
    plt.title('delta f')
    plt.xlabel(abscissa_label)
    plt.ylabel(ordinate_label)
    plotNum += 1

# Plot f0
if delta_f_available and full_f_available:
    plt.subplot(numRows, numCols, plotNum)
    plt.contourf(abscissa,ordinate,full_f_data - delta_f_data, numContours)
    plt.colorbar()
    plt.title('f_0 = full f - delta f')
    plt.xlabel(abscissa_label)
    plt.ylabel(ordinate_label)
    plotNum += 1


ax = fig.add_axes([0,0,1,1], frameon=False)
#titleString = "Plot generated by "+ inspect.getfile(inspect.currentframe()) + "\nRun in "+os.getcwd()
titleString = "Plot generated by "+ this_filename + "\nRun in "+os.getcwd()
ax.text(0.5,0.99,titleString,horizontalalignment='center',verticalalignment='top')

###############################################3
plt.show()
