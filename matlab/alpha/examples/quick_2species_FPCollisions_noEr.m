clear
addpath('../src')

sfincs_defaults()

% Change any parameters from the defaults here:

global geometryScheme 
geometryScheme = 4;

global Zs mHats nHats THats dnHatdrHats dTHatdrHats
Zs = [1,6];
mHats = [1,6];
nHats = [0.6, 0.009];
THats = [0.5, 0.8];
dnHatdrHats = [-0.587199, -0.00195733];
dTHatdrHats = [-0.587199, -0.391466];

global Delta gamma nu_n
Delta = 4.5694e-3;
gamma = 1;
nu_n = 8.4774e-3;

global Nalpha Nzeta Nxi Nx
Nalpha = 5;
Nzeta = 7;
Nxi = 8;
Nx = 5;

%global xGridScheme
%xGridScheme = 5;
%xGrid_k = 1;

sfincs_main()

directory = ['../../../fortran/alpha/examples/',mfilename];
sfincs_compareToFortran(fullfile(directory,'sfincsOutput.h5'))
sfincs_compareMatricesAndVectorsToFortran(directory)
