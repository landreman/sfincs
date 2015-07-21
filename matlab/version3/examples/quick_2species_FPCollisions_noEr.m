clear
addpath('../src')

sfincs_defaults()

% Change any parameters from the defaults here:

global geometryScheme inputRadialCoordinateForGradients
geometryScheme = 4;
inputRadialCoordinateForGradients = 1;

global Zs mHats nHats THats dnHatdpsiNs dTHatdpsiNs
Zs = [1,6];
mHats = [1,6];
nHats = [0.6, 0.009];
THats = [0.5, 0.8];
dnHatdpsiNs = [-0.3 -0.001];
dTHatdpsiNs = [-0.3, -0.2];

global Delta alpha nu_n
Delta = 4.5694e-3;
alpha = 1;
nu_n = 8.4774e-3;

global Ntheta Nzeta Nxi Nx
Ntheta = 5;
Nzeta = 7;
Nxi = 8;
Nx = 5;

%global xGridScheme
%xGridScheme = 5;
%xGrid_k = 1;

sfincs_main()

directory = ['../../../fortran/version3/examples/',mfilename];
sfincs_compareToFortran(fullfile(directory,'sfincsOutput.h5'))
sfincs_compareMatricesAndVectorsToFortran(directory)
