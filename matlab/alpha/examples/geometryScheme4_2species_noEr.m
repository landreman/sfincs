clear
addpath('../src')

sfincs_defaults()

% Change any parameters from the defaults here:

global geometryScheme inputRadialCoordinateForGradients
geometryScheme = 4;
inputRadialCoordinateForGradients = 2;

global Zs mHats nHats THats dnHatdrHats dTHatdrHats
Zs = [1,-1];
mHats = [1,0.000545509];
nHats = [1,1];
THats = [1,1];
dnHatdrHats = [-0.5,-0.5];
dTHatdrHats = [-2,-2];

global Delta gamma nu_n
Delta = 4.5694e-3;
gamma = 1;
nu_n = 8.31565e-3;

global Nalpha Nzeta Nxi Nx
Nalpha = 13;
Nzeta = 23;
Nxi = 48;
Nx = 5;

sfincs_main()


directory = ['../../../fortran/alpha/examples/',mfilename];
sfincs_compareToFortran(fullfile(directory,'sfincsOutput.h5'))
sfincs_compareMatricesAndVectorsToFortran(directory)
