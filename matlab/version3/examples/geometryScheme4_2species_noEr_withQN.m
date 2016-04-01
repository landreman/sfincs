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

global Delta alpha nu_n includePhi1 includePhi1InKineticEquation
Delta = 4.5694e-3;
alpha = 1;
nu_n = 8.31565e-3;

includePhi1 = true;
includePhi1InKineticEquation = false;

global Ntheta Nzeta Nxi Nx
Ntheta = 13;
Nzeta = 19;
Nxi = 48;
Nx = 5;

sfincs_main()


directory = ['../../../fortran/version3/examples/',mfilename];
sfincs_compareToFortran(fullfile(directory,'sfincsOutput.h5'))
sfincs_compareMatricesAndVectorsToFortran(directory)
