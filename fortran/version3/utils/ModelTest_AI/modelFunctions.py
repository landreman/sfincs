#! 
# Description:
#************* 
# Python script containing various functions used to calculate the sine and cosine parts of the perturbed 
# impurity density. All functions are based on the following paper
# http://aip.scitation.org/doi/abs/10.1063/1.873593
#
# Created by:  Aylwin Iantchenko (07-12-2017)
##########################################################################################################################################
# Import packages
##########################################################################################################################################
import sys
from numpy import multiply
from termcolor import colored
from math import log,sqrt
from usefulFun import *
##########################################################################################################################################
# Initialize constants
##########################################################################################################################################

# General constants
pi=3.14159265359
eps0 = 8.8541878176203e-12
c = 2.99792458e08
keVinJoule =  1000*1.6021766e-19
e = 1.6021766e-19
me = 0.0002723085


# SFINCS Reference values
Bbar = 1.0 # T
Rbar = 1.0 # m
nbar = 1e19 # m^-3
mbar = 3.345244e-27 # kg
Tbar = 1.0 #in keV
phibar = 1.0 # in kV

# Other parameters
ni1 = 4.0 # reference density
R0 = 3.0
Rovera = 3.0
vbar = sqrt(2.0*keVinJoule*Tbar/mbar)
iIDX = 0 # ion index
zIDX = 1 # impurity index

M0 = 0 # assuming no rotation

##########################################################################################################################################
# Functions to calculate analythical formula for the perturbed density
##########################################################################################################################################
def getPert(path,nun, epsilont, mHats,THats,Zs,IHat,GHat,B0OverBBar,rN,nHats,dnHatdrN,dTHatdrN,adiabaticTHat,adiabaticNHat,Delta_SFINCS,psiLCFS,r,dPhiHatdrN,varPlot,q):
# Description: 
# This is the main function that calls all the other
# to calculate the sine and cosine parts of the perturbed impurity density

# Initialize array
 Vars = []

 # read variables
 # If path provided, read variables from file, otherwise use the values provided in input. 
 if path:
  nun, epsilont, mHats,THats,Zs,IHat,GHat,B0OverBBar,rN,nHats,dnHatdrN,dTHatdrN,adiabaticTHat,adiabaticNHat,Delta_SFINCS,dPhiHatdrN,iota = readH5(path,['nu_n', 'epsilon_t', 'mHats','THats','Zs','IHat','GHat','B0OverBBar','rN','nHats','dnHatdrN','dTHatdrN','adiabaticTHat','adiabaticNHat','Delta','dPhiHatdrN','iota'],True) # get normalised quantities
  Ntheta,Nx = readH5(path,['Ntheta', 'Nx'],True)

  q = 1.0/iota

 # To simplify further calculations
 rovera = epsilont*Rovera
 r = epsilont*R0
 a = r/rovera

 # convert to SI units,
 nubar,ms,Ts,Zs,I,G,B,rN,ns,dndPsi,dTdPsi,Te,ne,r,psiLCFS,dPhidrN =  norm2SI(nun, mHats,THats,Zs,IHat,GHat,B0OverBBar,rN,nHats,dnHatdrN,dTHatdrN,adiabaticTHat,adiabaticNHat,r,a,dPhiHatdrN)

 # Separate between main ion and impurity parameters for simplicity
 Zi,Ti,mi,ni,dndPsii,dTdPsii = Zs[iIDX],Ts[iIDX], ms[iIDX], ns[iIDX], dnHatdrN[iIDX], dTHatdrN[iIDX] # ion
 Zz,Tz,mz,nz,dndPsiz,dTdPsiz = Zs[zIDX],Ts[zIDX], ms[zIDX], ns[zIDX], dnHatdrN[zIDX], dTHatdrN[zIDX] # impurity
 
 # Calculate parameters
 B_poloidal = getBpol(I,r) # Poloidal magnetic field
 Itor = getItor(R0,Rbar,r,rN,B,q) # get I_t
 tauIZ,nuIz = getTauIZ(nubar,Zi,Zz,nHats[zIDX],mHats[iIDX],THats[iIDX]) # collision frequency
 Lperp = getLperp(dndPsii,Ti,dTdPsii,nHats[iIDX],Itor) # perpendicular gradient scale length
 g = calcg(mi,ni,Ti,I,tauIZ,nz,dndPsii,dTdPsii,r,B_poloidal,Lperp) # get g (Eq. 24 in the mentioned paper)
 alpha = getAlpha(Zz,nz,Te,Ti,ni,ne) # get alpha
 gamma = getGamma(Ti,nHats[iIDX],dTdPsii,dndPsii,Itor,Lperp) # get gamma, (Eq. 24 in the mentioned paper)

# get sinus and cosinus parts of the perturbed density
 ns,nc = calcPert(epsilont,g,alpha,gamma,M0) 
 vti = getVthermal(Ti,mi)

# If provided, execute the expression given in "varPlot". This is used to get the x-variable when plotting using the script "testModel.py"
 for val in varPlot:
  exec('Vars.append(' + val + ')') 

 return tuple(Vars)

   
def norm2SI(nu_n, ms,Ts,Zs,I,G,B,rN,ns,dnHatdrN,dTHatdrN,Te,ne,r,a,dPhiHatdrN):
# Description: 
# Convert dimensionless input to SI
 psiLCFS = pow(a,2.0)/2.0 * B*Bbar
 dPhidrN = dPhiHatdrN*phibar*1000 # times 1000 since phibar in kV
 nubar = nu_n*vbar/Rbar
 ms = np.multiply(ms,mbar)
 Ts = np.multiply(Ts,Tbar)
 ns = np.multiply(ns,nbar)
 I = I * Rbar*Bbar
 G = G * Rbar*Bbar
 B = B*Bbar
 Te = Te*Tbar
 ne = ne*nbar
 r = r*Rbar
 dndPsi = dnHatdrN 
 dTdPsi = dTHatdrN 
 RbarOvera = Rovera/R0 
 return nubar, ms,Ts,Zs,I,G,B,rN,ns,dndPsi,dTdPsi,Te,ne,r,psiLCFS,dPhidrN

def getBpol(I,r):
 return I/r
 
def getItor(R0,Rbar,r,rN,B0,q):
 Itor = R0/Rbar*B0*rN/r*q/(r*B0) # toroidal component
 return Itor

def getLperp(dndPsii,Ti,dTdPsii,ni,Itor):
 dPdPsii = dndPsii*Ti + dTdPsii*ni 
 Lperp = 1.0/(-Itor*(dPdPsii/(ni*Ti) - 3.0/2.0*dTdPsii/Ti))
 return Lperp
 
def getTauIZ(nubar,Zi,Zs,nsHat,miHat,TiHat):
 nu_s = nubar*Zi**2.0*Zs**2.0*nsHat/(sqrt(miHat*(TiHat)**3.0))
 tau_s = 1.0/nu_s
 return tau_s,1.0/tau_s

def calcg(mi,ni,Ti,I,tauIZ,nz,dnHatdPsii,dTHatdPsii,r,B_poloidal,Lperp): 
 dlnPidPsi = 1.0/Ti*dTHatdPsii + 1.0/ni*dnHatdPsii 
 dlnTidPsi = 1.0/Ti*dTHatdPsii 
 g =  mi*ni/(Lperp*e*tauIZ*nz*B_poloidal/r)
 return g

def getAlpha(Zz,nz,Te,Ti,ni,ne):
 alpha = nz*Zz*Zz/((ne/Te+ni/Ti)*Ti)
 return alpha
 
def getGamma(Ti,ni,dTidPsi,dnidPsi,Itor,Lperp):
 c0 = 0.33 # See Phys. Plasmas 1999
 gamma = -Lperp*c0*Itor/Ti*dTidPsi
 return gamma 
 
def r2PsiDerivative(data,r):
 # function to transform r derivatives to psi derivatives
 out = []
 for val in data:
  out.append(val*2.0*r)
 return out 
 
def calcPert(eps,g,alpha,gamma,M0): 
 # function to compute the sine and cosine parts of the perturbed density
 ns = 2*eps*g*((1+alpha) + (1+gamma)*M0**2)/((1+alpha)**2 + (1+gamma)**2*g**2)
 nc = 2*eps*((1+alpha)*M0**2 - (1+gamma)*g**2)/((1+alpha)**2 + (1+gamma)**2*g**2)
 return ns,nc      
 
def getVthermal(T,m):
 vt = sqrt(2.0*T*keVinJoule/m)
 return vt
  
##########################################################################################################################################
# Functions to check SFINCS output with the analythical model
##########################################################################################################################################  

def CorrectNorm(dn,nIn):
 # function to correctly normalise the SFINCS perturbed density output
 # in order to be able to compare with the model
 nz = nIn[1]
 dnNormed = np.multiply(dn,1.0/nz)
 return dnNormed 

def getDensityParts(path,theta,dn):
 sign = 1
 # function to calculate the sine and
 # cosine parts of the perturbed impurity density, based on 
 # SFINCS output

 if path:
  theta,dn,nInput = readH5(path,['theta','densityPerturbation','nHats'],False)
  dn=dn[0][:,zIDX,-1]
  dn = CorrectNorm(dn,nInput)
  sign = -1 # if using SFINCS, ns has to be multiplied with -1
# Get cosine and sine functions
  sine,cosine = getTrigFuns(theta)

# Integrate dn*sin to get pi*ns and integrate dn*cos to get pi*nc
 ns = sign*OnedInt(multiply(dn,sine),theta)*1.0/(pi) # sine part
 nc = OnedInt(multiply(dn,cosine),theta)*1.0/(pi) # cosine part
 
 return ns, nc
