clear
addpath('../src')

sfincs_defaults()

% Change any parameters from the defaults here:

global geometryScheme inputRadialCoordinate rN_wish inputRadialCoordinateForGradients equilibriumFile
geometryScheme = 5;
inputRadialCoordinate = 3;
inputRadialCoordinateForGradients = 2;
rN_wish = 0.88;
equilibriumFile = '../../../equilibria/wout_w7x_standardConfig.nc';

global Zs mHats nHats THats dnHatdrHats dTHatdrHats
Zs = [1,-1,20];
mHats = [1,0.000545509,40];
nHats = [0.62,0.66,0.002];
THats = [1.1,1.3,1.6];
dnHatdrHats = [-15.0,-15.5,-0.025];
dTHatdrHats = [-12.0,-14.0,-16.0];

global Delta alpha nu_n dPhiHatdrHat includePhi1 includePhi1InKineticEquation
Delta = 4.5694e-3;
alpha = 1;
nu_n = 8.31565e-3;

dPhiHatdrHat = 8.5897;

%includePhi1 = true;
%includePhi1InKineticEquation=true;

global Ntheta Nzeta Nxi Nx solverTolerance maxNumNonlinearIterations
Ntheta = 9;
Nzeta = 17;
Nxi = 18;
Nx = 4;
%solverTolerance = 1e-7;
maxNumNonlinearIterations = 20;

global preconditioner_x preconditioner_xi
preconditioner_x = 0;
preconditioner_xi = 0;

sfincs_main()


directory = ['../../../fortran/version3/examples/',mfilename];
sfincs_compareToFortran(fullfile(directory,'sfincsOutput.h5'))
sfincs_compareMatricesAndVectorsToFortran(directory)
