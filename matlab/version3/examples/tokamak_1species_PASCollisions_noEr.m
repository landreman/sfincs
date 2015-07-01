clear
addpath('../src')

sfincs_defaults()

% Change any parameters from the defaults here:

global inputRadialCoordinate inputRadialCoordinateForGradients rN_wish
inputRadialCoordinate = 3;
inputRadialCoordinateForGradients = 1;
rN_wish = 0.3;

global B0OverBBar GHat IHat iota epsilon_t epsilon_h helicity_l helicity_n
global psiAHat aHat 

B0OverBBar = 1;
GHat = 1;
IHat = 0;
iota = 1.31;
epsilon_t = 0.1;
epsilon_h = 0;
helicity_l = 1;
helicity_n = 1;
psiAHat = 0.045;
aHat = 0.1;

global Zs mHats nHats THats dnHatdpsiNs dTHatdpsiNs
Zs = 1;
mHats = 1;
nHats = 1;
THats = 0.5;
dnHatdpsiNs = -1;
dTHatdpsiNs = -0.5;

global Delta alpha nu_n collisionOperator
Delta = 4.5694e-3;
alpha = 1;
nu_n = 8.4774e-3;
collisionOperator = 1;

global Ntheta Nzeta Nxi Nx
Ntheta = 21;
Nzeta = 1;
Nxi = 31;
Nx = 8;

sfincs_main()

directory = ['../../../fortran/version3/examples/',mfilename];
sfincs_compareToFortran(fullfile(directory,'sfincsOutput.h5'))
sfincs_compareMatricesAndVectorsToFortran(directory)
