clear
addpath('../src')

sfincs_defaults()

% Change any parameters from the defaults here:
global RHSMode plotB
RHSMode = 3;
plotB = true;

global geometryScheme equilibriumFile rN_wish
%geometryScheme = 5;
%equilibriumFile = 'C:\Users\landreman\Box Sync\MATLAB\20150601-01 Sfincs version 3\equilibria\wout_w7x_standardConfig.nc';
geometryScheme = 11;
%equilibriumFile = 'C:\Users\landreman\Box Sync\MATLAB\20150601-01 Sfincs version 3\equilibria\w7x_standardConfig.bc';
equilibriumFile = '../../../equilibria/w7x_standardConfig.bc';
rN_wish = 0.5;

global nuPrime EStar collisionOperator
nuPrime = 1.0;
EStar = 0.2;
collisionOperator = 1;

global Ntheta Nzeta Nxi
Ntheta = 17;
Nzeta = 31;
Nxi = 24;

sfincs_main()

sfincs_compareToFortran('C:\Users\landreman\Documents\Ubuntu\sfincsOutput11.h5')
