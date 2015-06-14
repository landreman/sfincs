function sfincs_defaults()

% See the user manual for documentation of these parameters.

% *************************************************************************
% general namelist
% *************************************************************************

global RHSMode solveSystem plotB

RHSMode = 1;
solveSystem = true;
plotB = false;

% *************************************************************************
% geometryParameters namelist
% *************************************************************************

global geometryScheme inputRadialCoordinate inputRadialCoordinateForGradients
global psiHat_wish psiN_wish rHat_wish rN_wish
global B0OverBBar GHat IHat iota epsilon_t epsilon_h epsilon_antisymm
global helicity_l helicity_n helicity_antisymm_l helicity_antisymm_n
global psiAHat aHat
global equilibriumFile VMECRadialOption min_Bmn_to_load
global force0RadialCurrentInEquilibrium

geometryScheme = 1;
inputRadialCoordinate = 3;
inputRadialCoordinateForGradients = 2;
psiHat_wish = -1;
psiN_wish = 0.25;
rHat_wish = -1;
rN_wish = 0.5;

B0OverBBar = 1;
GHat = 3.7481;
IHat = 0;
iota = 0.4542;
epsilon_t = -0.07053;
epsilon_h = 0.05067;
epsilon_antisymm = 0;
helicity_l = 2;
helicity_n = 10;
helicity_antisymm_l = 1;
helicity_antisymm_n = 0;
psiAHat = 0.15596;
aHat = 0.5585;

equilibriumFile = '';
VMECRadialOption = 1;
min_Bmn_to_load = 0;
force0RadialCurrentInEquilibrium = true;

% *************************************************************************
% speciesParameters namelist
% *************************************************************************

global Zs mHats nHats THats
global dnHatdpsiHats dTHatdpsiHats dnHatdpsiNs dTHatdpsiNs
global dnHatdrHats dTHatdrHats dnHatdrNs dTHatdrNs

Zs = 1;
mHats = 1;
nHats = 1;
THats = 1;

dnHatdpsiHats = 0;
dTHatdpsiHats = 0;
dnHatdpsiNs = 0;
dTHatdpsiNs = 0;
dnHatdrHats = 0;
dTHatdrHats = 0;
dnHatdrNs = 0;
dTHatdrNs = 0;

% *************************************************************************
% physicsParameters namelist
% *************************************************************************

global Delta alpha nu_n nuPrime EStar EParallelHat
global dPhiHatdpsiHat dPhiHatdpsiN dPhiHatdrHat dPhiHatdrN
global collisionOperator constraintScheme
global includeXDotTerm includeElectricFieldTermInXiDot useDKESExBDrift include_fDivVE_term
global includePhi1 nonlinear includeTemperatureEquilibrationTerm magneticDriftScheme

Delta = 4.5694e-3;
alpha = 1;
nu_n = 8.33e-3;
nuPrime = 1;
EStar = 0;
EParallelHat = 0;

dPhiHatdpsiHat = 0;
dPhiHatdpsiN = 0;
dPhiHatdrHat = 0;
dPhiHatdrN = 0;

collisionOperator = 0;
constraintScheme = -1;

includeXDotTerm = true;
includeElectricFieldTermInXiDot = true;
useDKESExBDrift = false;
include_fDivVE_term = false;
includePhi1 = false;
nonlinear = false;
includeTemperatureEquilibrationTerm = false;
magneticDriftScheme = 0;

% *************************************************************************
% resolutionParameters namelist
% *************************************************************************

global Ntheta Nzeta Nxi Nx solverTolerance NL 
global NxPotentialsPerVth xMax forceOddNthetaAndNzeta

Ntheta = 15;
Nzeta = 15;
Nxi = 16;
Nx = 5;
solverTolerance = 1e-6;
NL = 4;
NxPotentialsPerVth = 40;
xMax = 5;

forceOddNthetaAndNzeta = true;

% *************************************************************************
% otherNumericalParameters namelist
% *************************************************************************

global thetaDerivativeScheme zetaDerivativeScheme xGridScheme xGrid_k
global xPotentialsGridScheme useIterativeLinearSolver
global GMRES_restart GMRES_maxIterations
global nonlinearTolerance

thetaDerivativeScheme = 2;
zetaDerivativeScheme = 2;
xGridScheme = 5;
xGrid_k = 0;
xPotentialsGridScheme = 2;
useIterativeLinearSolver = true;

GMRES_restart = 200;
GMRES_maxIterations = 200;
nonlinearTolerance = 1e-6;
% *************************************************************************
% preconditionerOptions namelist
% *************************************************************************

global preconditioner_species preconditioner_x preconditioner_x_min_L
global preconditioner_theta preconditioner_theta_min_L preconditioner_zeta preconditioner_zeta_min_L
global preconditioner_xi reusePreconditioner

preconditioner_species = 1;
preconditioner_x = 1;
preconditioner_x_min_L = 0;
preconditioner_theta = 0;
preconditioner_theta_min_L = 0;
preconditioner_zeta = 0;
preconditioner_zeta_min_L = 0;
preconditioner_xi = 0;
reusePreconditioner = true;

% *************************************************************************
% constants
% *************************************************************************

global BLOCK_F BLOCK_QN BLOCK_PHI1_CONSTRAINT BLOCK_DENSITY_CONSTRAINT BLOCK_PRESSURE_CONSTRAINT BLOCK_F_CONSTRAINT

BLOCK_F = 9990;
BLOCK_QN = 9991;
BLOCK_PHI1_CONSTRAINT = 9992;
BLOCK_DENSITY_CONSTRAINT = 9993;
BLOCK_PRESSURE_CONSTRAINT = 9994;
BLOCK_F_CONSTRAINT = 9995;

end