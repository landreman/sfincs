#!/usr/bin/env python

# This python script plots the output of a single SFINCS run.

import matplotlib.pyplot as plt
import h5py
import numpy

try:
    f = h5py.File('sfincsOutput.h5','r')
except:
    print "Unable to open sfincsOutput.h5"
    exit(1)

# The expression [()] converts from an h5py dataset to a numpy ndarray:
theta = f['theta'][()]
zeta = f['zeta'][()]
BHat = f['BHat'][()]
DHat = f['DHat'][()]
BHat_sup_theta = f['BHat_sup_theta'][()]
BHat_sup_zeta = f['BHat_sup_zeta'][()]
BHat_sub_psi = f['BHat_sub_psi'][()]
BHat_sub_theta = f['BHat_sub_theta'][()]
BHat_sub_zeta = f['BHat_sub_zeta'][()]
iota = f['iota'][()]
Nspecies = f['Nspecies'][()]
totalDensity = f['totalDensity'][()]
Mach = f['MachUsingFSAThermalSpeed'][()]
totalPressure = f['totalPressure'][()]
Zs = f['Zs'][()]
jHat = f['jHat'][()]

NspeciesToPlot = Nspecies
if Nspecies>3:
    NspeciesToPlot = 3

###############################################3
# Plot the input magnetic field
###############################################3
plt.figure(1)

numRows = 2
numCols = 4

plt.subplot(numRows, numCols, 1)
plt.contourf(zeta,theta,BHat.transpose(),20)
plt.title('BHat')
plt.xlabel('zeta')
plt.ylabel('theta')
plt.colorbar()
# Plot a field line:
if iota>0:
    plt.plot([0,zeta.max()],[0,zeta.max()*iota])
else:
    plt.plot([0,zeta.max()],[-zeta.max()*iota,0])

plt.subplot(numRows, numCols, 2)
plt.contourf(zeta,theta,BHat.transpose(),20)
#plt.pcolormesh(zeta,theta,BHat.transpose())
plt.title('BHat')
plt.xlabel('zeta')
plt.ylabel('theta')
plt.colorbar()
for theta0 in theta:
    plt.plot(zeta, [theta0]*zeta.size,'.k')

plt.subplot(numRows, numCols, 3)
plt.contourf(zeta,theta,BHat_sup_theta.transpose(),20)
plt.title('BHat_sup_theta')
plt.xlabel('zeta')
plt.ylabel('theta')
plt.colorbar()

plt.subplot(numRows, numCols, 4)
plt.contourf(zeta,theta,BHat_sup_zeta.transpose(),20)
plt.title('BHat_sup_zeta')
plt.xlabel('zeta')
plt.ylabel('theta')
plt.colorbar()

plt.subplot(numRows, numCols, 5)
plt.contourf(zeta,theta,DHat.transpose(),20)
plt.title('DHat (inverse Jacobian)')
plt.xlabel('zeta')
plt.ylabel('theta')
plt.colorbar()

plt.subplot(numRows, numCols, 6)
plt.contourf(zeta,theta,BHat_sub_psi.transpose(),20)
plt.title('BHat_sub_psi')
plt.xlabel('zeta')
plt.ylabel('theta')
plt.colorbar()

plt.subplot(numRows, numCols, 7)
plt.contourf(zeta,theta,BHat_sub_theta.transpose(),20)
plt.title('BHat_sub_theta')
plt.xlabel('zeta')
plt.ylabel('theta')
plt.colorbar()

plt.subplot(numRows, numCols, 8)
plt.contourf(zeta,theta,BHat_sub_zeta.transpose(),20)
plt.title('BHat_sub_zeta')
plt.xlabel('zeta')
plt.ylabel('theta')
plt.colorbar()


###############################################3
# Plot the outputs for iteration 1
###############################################3
plt.figure(2)

numRows = NspeciesToPlot
numCols = 4

#if zeta.size>1:
iteration = 0
for ispecies in range(0,NspeciesToPlot):

    # Plot density
    plt.subplot(numRows, numCols, ispecies*numCols+1)
    plt.contourf(zeta,theta,totalDensity[:,:,ispecies,iteration].transpose(),20)
    plt.title('totalDensity for species {} (Z={})'.format(ispecies+1,Zs[ispecies]))
    plt.xlabel('zeta')
    plt.ylabel('theta')
    plt.colorbar()

    # Plot Mach #
    plt.subplot(numRows, numCols, ispecies*numCols+2)
    plt.contourf(zeta,theta,Mach[:,:,ispecies,iteration].transpose(),20)
    plt.title('Mach # for species {} (Z={})'.format(ispecies+1,Zs[ispecies]))
    plt.xlabel('zeta')
    plt.ylabel('theta')
    plt.colorbar()

    # Plot pressure
    plt.subplot(numRows, numCols, ispecies*numCols+3)
    plt.contourf(zeta,theta,totalPressure[:,:,ispecies,iteration].transpose(),20)
    plt.title('totalPressure for species {} (Z={})'.format(ispecies+1,Zs[ispecies]))
    plt.xlabel('zeta')
    plt.ylabel('theta')
    plt.colorbar()

# Plot parallel current
plt.subplot(numRows, numCols, 4)
plt.contourf(zeta,theta,jHat[:,:,iteration].transpose(),20)
plt.title('jHat (parallel current)')
plt.xlabel('zeta')
plt.ylabel('theta')
plt.colorbar()

plt.subplot(1,1,1)
plt.text(0.5,0.9,'Results for first iteration')

###############################################3
plt.show()
