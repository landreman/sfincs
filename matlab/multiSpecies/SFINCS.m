function SFINCS()

% SFINCS:
% The Stellarator Fokker-Planck Iterative Neoclassical Conservative Solver.
% Multiple-species version.
% Written in 2013 by Matt Landreman
% Massachusetts Institute of Technology
% Plasma Science & Fusion Center

% Dimensional quantities in this program are normalized to "reference" values:
% \bar{B} = reference magnetic field, typically 1 Tesla.
% \bar{R} = reference length, typically 1 meter.
% \bar{n} = reference density, typically 10^19 m^{-3}, 10^20 m^{-3}, or something similar.
% \bar{m} = reference mass, typically either the mass of hydrogen or deuterium.
% \bar{T} = reference temperature in energy units, typically 1 eV or 1 keV.
% \bar{v} = \sqrt{2 * \bar{T} / \bar{m}} = reference speed
% \bar{Phi} = reference electrostatic potential, typically 1 V or 1 kV.

% You can choose any reference parameters you like, not just the values
% suggested here. The code "knows" about the reference values only through
% the 3 combinations Delta, alpha, and nu_n, input below.

% Radial gradients of density, temperature, and electrostatic potential are
% specified as derivatives with respect to psi_N, where psi_N is the 
% toroidal flux normalized to the value at the last closed flux surface. 
% (psi_N=0 is the magnetic axis, and psi_N=1 is the last closed flux 
% surface.)

% --------------------------------------------------
% Program control parameters:
% --------------------------------------------------

programMode = 1;
% 1 = single run.
% 2 = Do a convergence scan and save the results.
% 3 = Load a previous convergence scan and plot the results. (Doesn't do any new solves.)
% 4 = Do a nuPrime scan and save the results.
% 5 = Load a previous nuPrime scan and plot the results. (Doesn't do any new solves.)

% The setting below matters for programMode=3 or programMode=5 only:
dataFileToPlot = 'm20130318_02_SFINCS_2013-03-18_14-47_convergenceScan_convergenceScan.mat';

RHSMode = 1;
% 1 = Use a single right-hand side.
% 2 = Use multiple right-hand sides to compute the transport matrix.
% At present, RHSMode==2 is only allowed when a single species is used.

% The variable below is set to true only for rare testing
% circumstances. Typically it should be false.
%testQuasisymmetryIsomorphism = true;
testQuasisymmetryIsomorphism = false;

%saveAllVariablesInAFileUponCompletion = true;
saveAllVariablesInAFileUponCompletion = false;

% The string below is appended to the filename of the data file
% (but before the .mat extension.)
filenameNote = '';

% --------------------------------------------------
% Geometry parameters:
% --------------------------------------------------

geometryScheme = 12;
% 1 = Two-helicity model
% 2 = Three-helicity approximation of the LHD standard configuration
% 3 = Four-helicity approximation of the LHD inward-shifted configuration
% 4 = Three-helicity approximation of the W7-X Standard configuration
% 10= Read the boozer coordinate data from the file specified as "fort996boozer_file" below
% 11= Read the boozer coordinate data from the file specified as "JGboozer_file" below (stellarator symmetric file)
% 12= Read the boozer coordinate data from the file specified as "JGboozer_file_NonStelSym" below (non-stellarator symmetric file)

% Additional parameters used only when geometryScheme=1:
% B = BBar * B0OverBBar * [1 + epsilon_t * cos(theta) + epsilon_h * cos(helicity_l * theta - helicity_n * zeta)]
B0OverBBar = 0.7;
epsilon_t = 0.13;
epsilon_h = 0.1;
helicity_l = 2;
helicity_n = 5;
% iota is the rotational transform = 1 / (safety factor q)
iota = 1.31;
% G is c/2 * the poloidal current outside the flux
% surface. Equivalently, G is the coefficient of grad zeta in the
% covariant representation of vector B. GHat is G normalized by \bar{B}\bar{R}.
GHat = 1.0;
% I is c/2 * the toroidal current inside the flux
% surface. Equivalently, I is the coefficient of grad theta in the
% covariant representation of vector B. IHat is I normalized by \bar{B}\bar{R}.
IHat = 0.8;
% dGdpHat = G'/(\mu_0 p') \bar{B}/\bar{R} is only used for NTV calculation, see notes_NTV.pdf.
dGdpHat=NaN;

% geometryScheme=10 parameters:
fort996boozer_file='TJII-midradius_example_s_0493_fort.996';
% Note that PsiA is not stored in the fort.996 file, so we use the
% PsiAHat setting below

% geometryScheme=11 and 12 parameters:
JGboozer_file='w7x-sc1.bc'; % stellarator symmetric example, geometryScheme=11
JGboozer_file_NonStelSym='out_neo-2_2_axisym'; % non-stellarator symmetric example, geometryScheme=12, requires Nzeta=1
%JGboozer_file_NonStelSym='out_neo-2_n2_sym_c_m64_n16';
normradius_wish=0.5;   %The calculation will be performed for the radius
                       %closest to this one in the JGboozer_file(_NonStelSym)
min_Bmn_to_load=1e-6;  %Filter out any Bmn components smaller than this

% --------------------------------------------------
% Species parameters:
% --------------------------------------------------

% Zs          = charges of each species, in units of the proton charge e
% mHats       = masses of each species, normalized to the reference mass \bar{m}
% nHats       = densities of each species, normalized to the reference density \bar{n}
% THats       = temperatures of each species, normalized to the reference temperature \bar{T}
% dnHatdpsiNs = radial gradient of the density of each species with respect to the normalized toroidal flux psi_N, normalized to the reference density \bar{n}
% dTHatdpsiNs = radial gradient of the temperature of each species with respect to the normalized toroidal flux psi_N, normalized to the reference temperature \bar{T}

%{
% Here is an example for 1 species:
Zs = 1;
mHats = 1;
nHats = 1.0;
dnHatdpsiNs = -0.5;
THats = 0.1;
dTHatdpsiNs = -0.7;
%}

%{
% Here is an example for 2 species:
Zs = [1, 6];
mHats = [1, 6];
nHats = [0.6, 0.009];
dnHatdpsiNs = [-0.3, -0.001];
THats = [0.8, 0.8];
dTHatdpsiNs = [-0.2, -0.2];
%}

% Here is an example for 2 species:
Zs = [1, -1];
mHats = [1, 1/1836];
nHats = [0.6, 0.6];
dnHatdpsiNs = [-0.3, -0.3];
THats = [0.8, 0.8];
dTHatdpsiNs = [-0.2, -0.2];
%


% --------------------------------------------------
% Other physics parameters:
% --------------------------------------------------

% Roughly speaking, Delta is rho_* at the reference parameters.
% More precisely, 
% Delta = c * \bar{m} * \bar{v} / (e * \bar{B} * \bar{R}) in Gaussian units,
% Delta =     \bar{m} * \bar{v} / (e * \bar{B} * \bar{R}) in SI units,
% where
% c = speed of light
% e = proton charge
%Delta = 0.001;
Delta = 4.5694e-3; %reference values: \bar{T}=1 keV, \bar{n}=10^20 m^-3,
                   %\bar{Phi}=1 kV, \bar{B}=1 T, \bar{R}=1 m, \bar{m}=proton mass

% alpha = e * \bar{Phi} / \bar{T} (in both Gaussian and SI units)
% where again e = proton charge.
alpha = 1.0;  % \bar{T}=1 keV, \bar{Phi}=1 kV

% psiAHat = psi_a / (\bar{B} * \bar{R}^2) (in both Gaussian and SI units)
% where 2*pi*psi_a is the toroidal flux at the last closed flux surface
% (the surface where psi_N = 1.)
% The value of psiAHat here is over-written for geometryScheme = 2, 3, 4, 11 and 12.
psiAHat = 0.03;

% Inductive electric field (often 0).
% EParallelHat = <E dot B> * \bar{R} / (\bar{Phi} * \bar{B})  (in both Gaussian and SI units)
% where
% E = electric field vector
% B = magnetic field vector
% < ... > denotes a flux surface average.
EParallelHat = 0;

% Radial electric field.
dPhiHatdpsiN = 0.0;

% nu_n is the collisionality at the reference parameters.
% More precisely, nu_n = \bar{nu} * \bar{R} / \bar{v} (in both Gaussian and SI units)
% where \bar{nu} is the dimensional collision frequency at the reference parameters:
%
%                  4 * sqrt{2*pi} * \bar{n} * e^4 * ln(Lambda)
% \bar{nu} = -----------------------------------------------------------   (SI units)
%             3 * (4 * pi * epsilon_0)^2 * sqrt(\bar{m}} * \bar{T}^(3/2)
%
% or, equivalently,
%
%                  4 * sqrt{2*pi} * \bar{n} * e^4 * ln(Lambda)
% \bar{nu} = -----------------------------------------------------------   (Gaussian units)
%                       3 * sqrt(\bar{m}} * \bar{T}^(3/2)
%
%
% Notice that collisionality is defined differently in the single-species code!
%nu_n = 0.00831565d+0;
nu_n = 8.4774e-3;  %reference values: \bar{T}=1 keV, \bar{n}=10^20 m^-3, \bar{Phi}=1 kV, 
                   %\bar{B}=1 T, \bar{R}=1 m, \bar{m}=proton mass, ln(Lambda)=17.3

% If testQuasisymmetryIsomorphism is true, the value of nu_n is changed so
% the physical collisionality stays constant as the helicity is changed.

collisionOperator = 0;
% 0 = Full linearized Fokker-Planck operator
% 1 = Pitch angle scattering, with no momentum conservation

% Unless you know what you are doing, keep constraintScheme = -1.
constraintScheme = -1;
% -1 = Automatic: if collisionOperator==0 then use constraintScheme=1, otherwise use constraintScheme=2.
%  0 = No constraints
%  1 = 2 constraints per species: <n_1> = 0 and <p_1> = 0.
%  2 = Nx constraints per species: <f_1>=0 at each x.

% To use one of the 4 most common trajectory models, the remaining parameters
% in this section should be set as follows:
%
% Full trajectories:
%   includeXDotTerm = true
%   includeElectricFieldTermInXiDot = true
%   useDKESExBDrift = false
%   include_fDivVE_term = false
%
% Partial trajectories: (non-conservative, as defined in the paper.)
%   includeXDotTerm = false
%   includeElectricFieldTermInXiDot = false
%   useDKESExBDrift = false
%   include_fDivVE_term = false
%
% Conservative partial trajectories: (Not discussed in the paper.)
%   includeXDotTerm = false
%   includeElectricFieldTermInXiDot = false
%   useDKESExBDrift = false
%   include_fDivVE_term = true
%
% DKES trajectories:
%   includeXDotTerm = false
%   includeElectricFieldTermInXiDot = false
%   useDKESExBDrift = true
%   include_fDivVE_term = false

includeXDotTerm = true;
%includeXDotTerm = false;

includeElectricFieldTermInXiDot = true;
%includeElectricFieldTermInXiDot = false;

% If useDKESExBDrift=true, the ExB drift term in the df/dtheta and df/dzeta terms is taken
% to be E x B / <B^2> instead of E x B / B^2.
%useDKESExBDrift = true;
useDKESExBDrift = false;

%include_fDivVE_term = true;
include_fDivVE_term = false;
% If true, a term f_{s1} div (v_E) is included in the kinetic equation.
% This term may make sense to include with the partial trajectory model
% as it restores Liouville's theorem (particle conservation) and eliminates
% the need for either a particle or heat source.

% --------------------------------------------------
% Numerical resolution parameters:
% --------------------------------------------------

% For each of the quantities below, the 'Converged' value is used except
% when that quantity is being varied in a convergence scan, in which case
% each value in the array that follows (e.g. Nthetas, NLs, etc.) is used.

% Number of grid points in the poloidal direction.
% Memory and time requirements DO depend strongly on this parameter.
NthetaConverged = 5;
%Nthetas = NthetaConverged;
Nthetas = floor(linspace(5,22,9));

% Number of grid points in the toroidal direction
% (per identical segment of the stellarator.)
% Memory and time requirements DO depend strongly on this parameter.
NzetaConverged = 11;
%Nzetas = NzetaConverged;
Nzetas = floor(linspace(5,20,9));

% Number of Legendre polynomials used to represent the distribution
% function.
% Memory and time requirements DO depend strongly on this parameter.
% The value of this parameter required for convergence depends strongly on
% the collisionality. At high collisionality, this parameter can be as low
% as ~ 5. At low collisionality, this parameter may need to be many 10s or
% even > 100 for convergence.
NxiConverged = 8;
%Nxis = NxiConverged;
Nxis = floor(linspace(4,20,9));

% Number of Legendre polynomials used to represent the Rosenbluth
% potentials: (Typically 2 or 4 is plenty.)
% Memory and time requirements do NOT depend strongly on this parameter.
NLConverged = 4;
%NLs = NLConverged;
NLs = 4:1:6;

% Number of grid points in energy used to represent the distribution
% function.
% Memory and time requirements DO depend strongly on this parameter.
% This parameter almost always needs to be at least 5.
% Sometimes 5-7 is sufficient, but sometimes it needs to be as high as
% 10-15.  Typically this parameter needs to be higher as the number of
% species is increased.
NxConverged = 5;
Nxs=4:15;

% Number of grid points in energy used to represent the Rosenbluth
% potentials.
% Memory and time requirements do NOT depend strongly on this parameter.
NxPotentialsPerVthConverged = 40;
%NxPotentialsPerVths = NxPotentialsPerVthConverged ;
NxPotentialsPerVths = [40, 81];
%NxPotentialsPerVths = floor(linspace(20,80,7));

% Tolerance used to define convergence of the Krylov solver.
% This parameter does not affect memory requirements but it does affect the
% time required for solution.
log10tolConverged = 7;
%log10tols = log10tolConverged;
log10tols = 5:0.5:9;


% --------------------------------------------------
% Other numerical parameters:
% --------------------------------------------------

% For most production runs, you do want to use the Krylov solver rather
% than a direct solver, in which case this parameter should be "true".
% However, for low-resolution problems, or if the Krylov solvers are
% failing to converge, you might want to set this parameter to "false".
tryIterativeSolver = true;
%tryIterativeSolver = false;

orderOfSolversToTry = [1, 3, 4];
% 1 = GMRES
% 2 = BiCGStab
% 3 = BiCGStab(l)
% 4 = TFQMR
% 5 = CGS

% Below are some setting for the Krylov solvers:
maxIterations = 200;
restart = maxIterations; % Used only for GMRES.

%tryDirectSolverIfIterativeSolversFail = true;
tryDirectSolverIfIterativeSolversFail = false;

thetaGridMode = 2;
% 0 = uniform periodic spectral
% 1 = 2nd order uniform finite-difference
% 2 = 4th order uniform finite-difference
% 3 = 6th order uniform finite-difference
% This parameter should almost always be 2.

forceThetaParity = 1;
% 0 = either even or odd Ntheta is fine.
% 1 = force Ntheta to be odd.
% 2 = force Ntheta to be even.
% This parameter should almost always be 1.

% The variable below is used only for testing, and it should
% otherwise be 1.
factorToIncludeInFNormalization = 1;

% --------------------------------------------------
% Settings for the preconditioner:
% --------------------------------------------------

preconditioner_species = 1;
% 0 = keep full species coupling
% 1 = drop all cross-species coupling
% Recommended value: 1

preconditioner_x = 1;
% 0 = keep full x coupling.
% 1 = keep only diagonal in x.
% 2 = keep upper-triangular part of x.
% 3 = Keep tridiagonal terms in x.
% 4 = Keep diagonal and superdiagonal in x.
% This parameter should almost always be 1.

preconditioner_x_min_L = 2;
% The simplified x coupling is used in the preconditioner only when the
% Legendre index L is >= this value; otherwise the full x coupling is used
% in the preconditioner.  Set to 0 to precondition at every L.
% Usually, good values for this parameter are 0, 1, or 2.

preconditioner_xi = 0;
% 0 = Keep full xi coupling
% 1 = keep only tridiagonal terms in xi.
% Either 0 or 1 may be appropriate for this parameter.

preconditioner_xi_max_L = Inf;
% All L coupling is dropped for L >= this value.
% Recommended value: Inf

preconditioner_theta_min_L = Inf;
% The full d/dtheta matrix is used for L < this value.
% Set this to 0 if you don't want to use the full d/dtheta matrix in the
% preconditioner for any L.
% Recommended values: 0, 1, 2, or Inf

preconditioner_theta_max_L = Inf;
% All theta coupling is dropped for L >= this value.
% Recommended value: Inf

%preconditioner_theta_remove_cyclic = true;
preconditioner_theta_remove_cyclic = false;
% If true, the (1,end) and (end,1) elements of the d/dtheta matrix are
% removed in the preconditioner.
% Recommended value: false

preconditioner_zeta_min_L = Inf;
% The full d/dzeta matrix is used for L < this value.
% Set this to 0 if you don't want to use the full d/dzeta matrix in the
% preconditioner for any L.
% Recommended values: 0, 1, 2, or Inf

preconditioner_zeta_max_L = Inf;
% All theta coupling is dropped for L >= this value.
% Recommended value: Inf

%preconditioner_zeta_remove_cyclic = true;
preconditioner_zeta_remove_cyclic = false;
% If true, the (1,end) and (end,1) elements of the d/dzeta matrix are
% removed in the preconditioner.
% Recommended value: false

% --------------------------------------------------
% Plotting options:
% --------------------------------------------------

% The following offset is added to all figure numbers. It can sometimes
% be convenient to change this number if you want to save figures rather
% than over-write them when re-running the code.
figureOffset=20;

% If smoothFigures==true, contours plots of quantities on the (theta,zeta) grid will
% use interpolation to reduce jagged contours.
smoothFigures = true;
%smoothFigures = false;

plotSpeedGrid = true;
%plotSpeedGrid = false;

% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% End of the input parameters.
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------

% nuPrime is not used in the multi-species code.
nuPrime = 1.0;


%{
fprintf('For an analogous run in the single-species version of SFINCS, you should set the following parameters in that code:\n')
fprintf('Delta = %g\n',Delta*sqrt(mHats)/Zs)
fprintf('omega = %g\n',alpha*Delta*sqrt(mHats)/2)
fprintf('nu_n = %g\n',nu_n*nHats*Zs^4/(THats^(3/2)))
%}

% Validate some of the input parameters:

Nspecies = numel(Zs);
if numel(mHats) ~= Nspecies
    error('Number of mHats does not equal the number of Zs')
end
if numel(nHats) ~= Nspecies
    error('Number of nHats does not equal the number of Zs')
end
if numel(dnHatdpsiNs) ~= Nspecies
    error('Number of dnHatdpsiNs does not equal the number of Zs')
end
if numel(THats) ~= Nspecies
    error('Number of THats does not equal the number of Zs')
end
if numel(dTHatdpsiNs) ~= Nspecies
    error('Number of dTHatdpsiNs does not equal the number of Zs')
end

if Nspecies > 1 && RHSMode == 2
    error('RHSMode==2 is only allowed when a single species is used')
end

switch preconditioner_x
    case {0,1,2,3,4}
    otherwise
        error('Invalid setting for preconditioner_x')
end
switch preconditioner_xi
    case {0,1}
    otherwise
        error('Invalid setting for preconditioner_xi')
end

if constraintScheme < 0
    if collisionOperator == 0
        constraintScheme = 1;
    else
        constraintScheme = 2;
    end
end

if testQuasisymmetryIsomorphism
    if geometryScheme ~= 1
        error('To test the quasisymmetry isomorphism, you should set geometryScheme=1.')
    end
    
    if epsilon_t ~= 0
        error('To test the quasisymmetry isomorphism, you should set epsilon_t=0.')
    end
    
    if RHSMode ~= 1
        error('testQuasisymmetryIsomorphism and RHSMode=2 are not presently compatible.')
    end
    
    nuStarS = nu_n;
    nu_n = nuStarS*abs(helicity_n/iota - helicity_l);
end

FSADensityPerturbation = 0;
FSABFlow = 0;
FSAPressurePerturbation = 0;
particleFlux = 0;
momentumFlux = 0;
heatFlux = 0;
jHat = 0;
FSABjHat = 0;
transportMatrix = zeros(3);

speedGridFigureHandle = 0;
KrylovFigureHandle = 0;

iteration=0;

if plotSpeedGrid
    figure(figureOffset+7)
    clf
end

fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
fprintf('SFINCS: The Stellarator Fokker-Planck Iterative Neoclassical Conservative Solver.\n')

switch programMode
    case 1
        
        fprintf('Beginning a single run.\n')
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        Ntheta=NthetaConverged;
        Nzeta=NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        solveDKE();

    case 2
        fprintf('Beginning convergence scans.\n')
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        startTime=clock;
        
        switch RHSMode
            case 1
                quantitiesToRecord = {'FSABFlow','particleFlux','momentumFlux','heatFlux','NTV'};
            case 2
                quantitiesToRecord = {'L11','L12=L21','L13=L31','L12=L21','L22','L23=L32','L13=L31','L23=L32','L33'};
            otherwise
                error('Invalid RHSMode')
        end
            
        linespecs = {'.-b','.-r','.-g','.-m','.:c','.-r','.:k','.:b','.-m'};
        
        parametersToVary = {'N\theta','N\zeta','NL','N\xi','Nx','NxPotentialsPerVth','-log_{10}tol'};
        abscissae = {Nthetas, Nzetas, NLs, Nxis, Nxs, NxPotentialsPerVths, log10tols};
        convergeds = {NthetaConverged, NzetaConverged, NLConverged, NxiConverged, NxConverged, NxPotentialsPerVthConverged, log10tolConverged};
        
        numQuantities = numel(quantitiesToRecord);
        numParameters = numel(parametersToVary);
        quantities = cell(numParameters,1);
        quantities{1} = zeros(numel(Nthetas), numQuantities, Nspecies);
        quantities{2} = zeros(numel(Nzetas), numQuantities, Nspecies);
        quantities{3} = zeros(numel(NLs), numQuantities, Nspecies);
        quantities{4} = zeros(numel(Nxis), numQuantities, Nspecies);
        quantities{5} = zeros(numel(Nxs), numQuantities, Nspecies);
        quantities{6} = zeros(numel(NxPotentialsPerVths), numQuantities, Nspecies);
        quantities{7} = zeros(numel(log10tols), numQuantities, Nspecies);
        parameterScanNum = 1;
        
        % Vary Ntheta, keeping other numerical parameters fixed.
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(Nthetas)
            Ntheta=Nthetas(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1,:)=FSABFlow;
                    quantities{parameterScanNum}(iii,2,:)=particleFlux;
                    quantities{parameterScanNum}(iii,3,:)=momentumFlux;
                    quantities{parameterScanNum}(iii,4,:)=heatFlux;
                    quantities{parameterScanNum}(iii,5,:)=NTV;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary Nzeta, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(Nzetas)
            Nzeta=Nzetas(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1,:)=FSABFlow;
                    quantities{parameterScanNum}(iii,2,:)=particleFlux;
                    quantities{parameterScanNum}(iii,3,:)=momentumFlux;
                    quantities{parameterScanNum}(iii,4,:)=heatFlux;
                    quantities{parameterScanNum}(iii,5,:)=NTV;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary NL, keeping other numerical parameters fixed.
        Nzeta = NzetaConverged;
        Ntheta=NthetaConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(NLs)
            NL=NLs(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1,:)=FSABFlow;
                    quantities{parameterScanNum}(iii,2,:)=particleFlux;
                    quantities{parameterScanNum}(iii,3,:)=momentumFlux;
                    quantities{parameterScanNum}(iii,4,:)=heatFlux;
                    quantities{parameterScanNum}(iii,5,:)=NTV;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary Nxi, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(Nxis)
            Nxi=Nxis(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1,:)=FSABFlow;
                    quantities{parameterScanNum}(iii,2,:)=particleFlux;
                    quantities{parameterScanNum}(iii,3,:)=momentumFlux;
                    quantities{parameterScanNum}(iii,4,:)=heatFlux;
                    quantities{parameterScanNum}(iii,5,:)=NTV;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        
        % Vary Nx, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(Nxs)
            Nx = Nxs(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1,:)=FSABFlow;
                    quantities{parameterScanNum}(iii,2,:)=particleFlux;
                    quantities{parameterScanNum}(iii,3,:)=momentumFlux;
                    quantities{parameterScanNum}(iii,4,:)=heatFlux;
                    quantities{parameterScanNum}(iii,5,:)=NTV;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary NxPotentialsPerVth, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        tol = 10^(-log10tolConverged);
        for iii = 1:numel(NxPotentialsPerVths)
            NxPotentialsPerVth = NxPotentialsPerVths(iii);
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1,:)=FSABFlow;
                    quantities{parameterScanNum}(iii,2,:)=particleFlux;
                    quantities{parameterScanNum}(iii,3,:)=momentumFlux;
                    quantities{parameterScanNum}(iii,4,:)=heatFlux;
                    quantities{parameterScanNum}(iii,5,:)=NTV;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        % Vary tol, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        Nzeta = NzetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        for iii = 1:numel(log10tols)
            tol = 10^(-log10tols(iii));
            solveDKE()
            switch RHSMode
                case 1
                    quantities{parameterScanNum}(iii,1,:)=FSABFlow;
                    quantities{parameterScanNum}(iii,2,:)=particleFlux;
                    quantities{parameterScanNum}(iii,3,:)=momentumFlux;
                    quantities{parameterScanNum}(iii,4,:)=heatFlux;
                    quantities{parameterScanNum}(iii,5,:)=NTV;
                case 2
                    quantities{parameterScanNum}(iii,:)=reshape(transportMatrix,[9,1]);
            end
        end
        parameterScanNum = parameterScanNum+1;
        
        temp=dbstack;
        nameOfThisProgram=sprintf('%s',temp.file);
        filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_convergenceScan_',filenameNote];
        outputFilename=[filenameBase,'.mat'];
        if saveAllVariablesInAFileUponCompletion
            save(outputFilename)
        end

        plotConvergenceScan()

    case 3
        % Plot a previous convergence scan:
        
        load(dataFileToPlot)
        plotConvergenceScan()

    case 4
        % Do a nuPrime scan:
        
        if RHSMode ~= 2
            error('programMode=4 requires RHSMode=2')
        end
        
        fprintf('Beginning nuPrime scans.\n')
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        %numNus = 13;
        %nuPrimes = logspace(-1,2,numNus);
        
        numNus = 17;
        nuPrimes = logspace(-2,2,numNus);

        
        referenceNuPrimes = [0.01, 0.1, 0.3,  1, 10, 100];
        referenceNthetas =  [  15,  15,  15, 15, 15,  15];
        referenceNzetas =   [  15,  13,  13, 13, 13,  13];
        referenceNxis =     [  48,  37,  34, 13, 13,  13];
        referenceNxs =      [   5,   5,   6,  8,  8,   8];
        
        %{
        referenceNuPrimes = [0.01, 0.1];
        referenceNthetas =  [  15,  15];
        referenceNzetas =   [  15,  13];
        referenceNxis =     [  48,  37];
        referenceNxs =      [   5,   5];
        %}
        %{
        referenceNuPrimes = [0.1, 0.3,  1, 10, 100];
        referenceNthetas =  [ 19,  19, 15, 15,  15];
        referenceNzetas =   [ 13,  13, 13, 13,  13];
        referenceNxis =     [ 37,  34, 13, 13,  13];
        referenceNxs =      [  7, 7.5,  8,  8,   8];
        %}
        %{
        referenceNuPrimes = [0.1,  1, 10, 100];
        referenceNthetas =  [ 13,  9,  7,   7];
        referenceNzetas =   [  9,  9,  7,   7];
        referenceNxis =     [ 11,  8,  6,   6];
        referenceNxs =      [  4,  4,  5,   5];
        %}
        
        numReferenceNus = numel(referenceNuPrimes);
        if numel(referenceNthetas) ~= numReferenceNus || numel(referenceNzetas) ~= numReferenceNus || numel(referenceNxis) ~= numReferenceNus || numel(referenceNxs) ~= numReferenceNus
            error('Number of reference points for resolution is not consistent')
        end
        
        NthetaMultipliers = [1, 2, 1, 1, 1];
        NzetaMultipliers =  [1, 1, 2, 1, 1];
        NxiMultipliers    = [1, 1, 1, 2, 1];
        NxMultipliers     = [1, 1, 1, 1, 1.7];
    
        NConvergence = numel(NthetaMultipliers);
        if numel(NzetaMultipliers) ~= NConvergence || numel(NxiMultipliers) ~= NConvergence || numel(NxMultipliers) ~= NConvergence
            error('Sizes of multiplier arrays are not consistent')
        end
        
        numRuns = 3*NConvergence*numNus;
        scanResults = zeros(3, NConvergence, numNus, 3, 3);
        
        runNum=0;
        
        Nthetas = interp1(log10(referenceNuPrimes), referenceNthetas, log10(nuPrimes),'cubic');
        Nzetas = interp1(log10(referenceNuPrimes), referenceNzetas, log10(nuPrimes),'cubic');
        Nxis = interp1(log10(referenceNuPrimes), referenceNxis, log10(nuPrimes),'cubic');
        Nxs = interp1(log10(referenceNuPrimes), referenceNxs, log10(nuPrimes),'cubic');
        
        NL=NLConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);

        for iConvergence = 1:NConvergence
            for iCollisionOperator = 1:3
                for iNu = 1:numNus
                    runNum = runNum + 1;
                    fprintf('##################################################################\n')
                    fprintf('Run %d of %d.\n', runNum, numRuns)
                    fprintf('##################################################################\n')
                    
                    nuPrime = nuPrimes(iNu);
                    if nuPrime > 1
                        preconditioner_x_min_L = 1;
                    else
                        preconditioner_x_min_L = 0;
                    end
                    
                    collisionOperator = iCollisionOperator-1;
                    Ntheta=floor(Nthetas(iNu)*NthetaMultipliers(iConvergence));
                    Nzeta=floor(Nzetas(iNu)*NzetaMultipliers(iConvergence));
                    Nxi=floor(Nxis(iNu)*NxiMultipliers(iConvergence));
                    Nx=floor(Nxs(iNu)*NxMultipliers(iConvergence));
                    
                    solveDKE()
                    
                    scanResults(iCollisionOperator, iConvergence, iNu, :, :) = transportMatrix;
                end
            end
        end
        
        temp=dbstack;
        nameOfThisProgram=sprintf('%s',temp.file);
        filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_nuPrimeScan_',filenameNote];
        outputFilename=[filenameBase,'.mat'];
        if saveAllVariablesInAFileUponCompletion
            save(outputFilename)
        end

        plotNuPrimeScan()
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Done with nuPrime scans.\n')
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    case 5
        % Plot results of a previous nuPrime scan:
        
        load(dataFileToPlot)
        plotNuPrimeScan()
        
    otherwise
        error('Invalid setting for programMode.')
end

    function plotConvergenceScan()
        % ------------------------------------------------------
        % Plot results of the convergence scan:
        % ------------------------------------------------------
        
        quantityToRowMap = [1, 2, 3, 2, 4, 5, 3, 5, 6];
        switch RHSMode
            case 1
                numRows = numQuantities;
            case 2
                numRows = 6;
        end
        numCols = numParameters;
        
        for ispecies=1:Nspecies
            maxs=ones(numQuantities,1)*(-Inf);
            mins=ones(numQuantities,1)*(Inf);
            for iParameter = 1:numParameters
                maxs = max([maxs, quantities{iParameter}(:,:,ispecies)'],[],2);
                mins = min([mins, quantities{iParameter}(:,:,ispecies)'],[],2);
            end
            
            figure(1+figureOffset+(ispecies-1)*2)
            clf
            for iQuantity = 1:numQuantities
                if maxs(iQuantity) <= mins(iQuantity)
                    maxs(iQuantity) = mins(iQuantity)+1;
                end
                switch RHSMode
                    case 1
                        iRow = iQuantity;
                    case 2
                        iRow = quantityToRowMap(iQuantity);
                end
                for iParameter = 1:numParameters
                    subplot(numRows, numCols, iParameter  + (iRow - 1)*numParameters)
                    plot(1./abscissae{iParameter}, quantities{iParameter}(:,iQuantity,ispecies)', linespecs{iQuantity})
                    hold on
                    plot(1./[convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
                    ylim([mins(iQuantity), maxs(iQuantity)])
                    xlabel(['1/',parametersToVary{iParameter}])
                    ylabel(quantitiesToRecord{iQuantity})
                end
            end
            switch RHSMode
                case 1
                    stringForTop = sprintf('SFINCS convergence scan: species %d, nu_n=%g. Base case: Ntheta=%d, Nzeta=%d, NL=%d, Nxi=%d, Nx=%d, NxPotentialsPerVth=%g, -log10tol=%g.', ...
                        ispecies, nu_n, NthetaConverged, NzetaConverged, NLConverged, NxiConverged, NxConverged, NxPotentialsPerVthConverged, log10tolConverged);
                case 2
                    stringForTop = sprintf('SFINCS convergence scan: species %d, nuPrime=%g. Base case: Ntheta=%d, Nzeta=%d, NL=%d, Nxi=%d, Nx=%d, NxPotentialsPerVth=%g, -log10tol=%g.', ...
                        ispecies, nuPrime, NthetaConverged, NzetaConverged, NLConverged, NxiConverged, NxConverged, NxPotentialsPerVthConverged, log10tolConverged);
            end
            annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
                'Interpreter','none','VerticalAlignment','bottom',...
                'FontSize',12,'LineStyle','none','String',stringForTop);
            
            figure(2+figureOffset+(ispecies-1)*2)
            clf
            for iQuantity = 1:numQuantities
                if maxs(iQuantity) <= mins(iQuantity)
                    maxs(iQuantity) = mins(iQuantity)+1;
                end
                switch RHSMode
                    case 1
                        iRow = iQuantity;
                    case 2
                        iRow = quantityToRowMap(iQuantity);
                end
                for iParameter = 1:numParameters
                    subplot(numRows, numCols, iParameter  + (iRow - 1)*numParameters)
                    plot(abscissae{iParameter}, quantities{iParameter}(:,iQuantity,ispecies)', linespecs{iQuantity})
                    hold on
                    plot([convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
                    ylim([mins(iQuantity), maxs(iQuantity)])
                    xlabel(parametersToVary{iParameter})
                    ylabel(quantitiesToRecord{iQuantity})
                end
            end
            
            annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
                'Interpreter','none','VerticalAlignment','bottom',...
                'FontSize',12,'LineStyle','none','String',stringForTop);
        end
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time for convergence scans: %g seconds.\n',etime(clock,startTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
    end

    function plotNuPrimeScan()
        legendText = {'base case','2x N\theta','2x N\zeta','2x N\xi','2x Nx'};
        colors = [  1,0,0;
                    0.7,0.5,0;
                    0,0.7,0;
                    0,0,1;
                    1,0,1];
        
        for iCollision = 1:3
            figure(iCollision+figureOffset)
            clf
            numRows = 2;
            numCols = 3;
            plotNum = 1;
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 1, 1))), '.-','Color',colors(iConvergence,:))
                hold on
            end
            xlabel('nuPrime')
            ylabel('-L11')
            legend(legendText)
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 1, 2))), '.-','Color',colors(iConvergence,:))
                hold on
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 2, 1))), '.:','Color',colors(iConvergence,:))
            end
            xlabel('nuPrime')
            ylabel('L12=L21')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                semilogx(nuPrimes, squeeze(scanResults(iCollision,iConvergence, :, 1, 3)), '.-','Color',colors(iConvergence,:))
                hold on
                semilogx(nuPrimes, squeeze(scanResults(iCollision,iConvergence, :, 3, 1)), '.:','Color',colors(iConvergence,:))
            end
            xlabel('nuPrime')
            ylabel('L13=L31')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 2, 2))), '.-','Color',colors(iConvergence,:))
                hold on
            end
            xlabel('nuPrime')
            ylabel('-L22')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                semilogx(nuPrimes, squeeze(scanResults(iCollision,iConvergence, :, 2, 3)), '.-','Color',colors(iConvergence,:))
                hold on
                semilogx(nuPrimes, squeeze(scanResults(iCollision,iConvergence, :, 3, 2)), '.:','Color',colors(iConvergence,:))
            end
            xlabel('nuPrime')
            ylabel('L23=L32')
            
            subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
            for iConvergence = 1:NConvergence
                loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 3, 3))), '.-','Color',colors(iConvergence,:))
                hold on
            end
            xlabel('nuPrime')
            ylabel('L33')
        end
        
        legendText = {'Fokker-Planck','pitch-angle scattering','momentum-conserving model'};
        colors = [  1,0,0;
                    0,0.7,0;
                    0,0,1];

        figure(4+figureOffset)
        clf
        numRows = 2;
        numCols = 3;
        plotNum = 1;
        iConvergence = 1;
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 1, 1))), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('-L11')
        legend(legendText)
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            data = 0.5*(scanResults(iCollision,iConvergence, :, 1, 2) + scanResults(iCollision,iConvergence, :, 2, 1));
            loglog(nuPrimes, abs(squeeze(data)), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('L12=L21')
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            data = 0.5*(scanResults(iCollision,iConvergence, :, 1, 3) + scanResults(iCollision,iConvergence, :, 3, 1));
            %loglog(nuPrimes, abs(squeeze(data)), '.-','Color',colors(iCollision,:))
            semilogx(nuPrimes, squeeze(data), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('L13=L31')
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 2, 2))), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('-L22')
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            data = 0.5*(scanResults(iCollision,iConvergence, :, 3, 2) + scanResults(iCollision,iConvergence, :, 2, 3));
            %loglog(nuPrimes, abs(squeeze(data)), '.-','Color',colors(iCollision,:))
            semilogx(nuPrimes, squeeze(data), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('L23=L32')
        
        subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
        for iCollision = 1:3
            loglog(nuPrimes, abs(squeeze(scanResults(iCollision,iConvergence, :, 3, 3))), '.-','Color',colors(iCollision,:))
            hold on
        end
        xlabel('nuPrime')
        ylabel('L33')
    end

% --------------------------------------------------------
% --------------------------------------------------------
% Done with the routines for convergence scans.
% Next comes the core function of the code.
% --------------------------------------------------------
% --------------------------------------------------------

    function solveDKE()
        
        startTimeForThisRun=tic;
        
        sqrtpi=sqrt(pi);
        iteration = iteration+1;
        
        
        % --------------------------------------------------------
        % --------------------------------------------------------
        % First, set up the grids, differentiation matrices, and
        % integration weights for each coordinate.
        % --------------------------------------------------------
        % --------------------------------------------------------
        
        
        switch forceThetaParity
            case 0
                % Do nothing
            case 1
                % For Ntheta to be odd
                if mod(Ntheta,2)==0
                    Ntheta=Ntheta+1;
                end
            case 2
                % For Ntheta to be even
                if mod(Ntheta,2)==1
                    Ntheta=Ntheta+1;
                end
            otherwise
                error('Invalid forceThetaParity')
        end
        
        if mod(Nzeta,2)==0
            Nzeta=Nzeta+1;
        end
        
        if iteration>1
            fprintf('********************************************************************\n')
        end
        fprintf('Ntheta = %d,  Nzeta = %d,  NL = %d,  Nxi = %d,  Nx = %d, NxPtentialsPerVth = %g, tol = %g\n',Ntheta, Nzeta,NL,Nxi,Nx,NxPotentialsPerVth,tol)
        
        tic
        
        % Generate abscissae, quadrature weights, and derivative matrix for theta grid.
        if Ntheta == 1
            theta = 0;
            thetaWeights = 2*pi;
            ddtheta = 0;
            ddtheta_preconditioner = 0;
        else
            switch thetaGridMode
                case 0
                    % Spectral uniform
                    scheme = 20;
                case 1
                    % Uniform periodic 2nd order FD
                    scheme = 0;
                case 2
                    % Uniform periodic 4th order FD
                    scheme = 10;
                case 3
                    % Uniform periodic 6th order FD
                    scheme = 70;
                otherwise
                    error('Error! Invalid thetaGridMode')
            end
            [theta, thetaWeights, ddtheta, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
            
            scheme = 0;
            [~, ~, ddtheta_preconditioner, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
            if preconditioner_theta_remove_cyclic
                ddtheta_preconditioner(1,end) = 0;
                ddtheta_preconditioner(end,1) = 0;
            end
           
        end
        
        
        % Generate abscissae, quadrature weights, and derivative matrix for zeta grid.
        setNPeriods()
        zetaMax = 2*pi/NPeriods;
        
        if Nzeta==1
            zeta=0;
            zetaWeights=2*pi;
            ddzeta=0;
            ddzeta_preconditioner=0;
        else
            switch thetaGridMode
                case 0
                    % Spectral uniform
                    scheme = 20;
                case 1
                    % Uniform periodic 2nd order FD
                    scheme = 0;
                case 2
                    % Uniform periodic 4th order FD
                    scheme = 10;
                case 3
                    % Uniform periodic 6th order FD
                    scheme = 70;
                otherwise
                    error('Error! Invalid thetaGridMode')
            end
            [zeta, zetaWeights, ddzeta, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Nzeta, 0, zetaMax, scheme);
            zetaWeights = zetaWeights * NPeriods;
            
            scheme = 0;
            [~, ~, ddzeta_preconditioner, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Nzeta, 0, zetaMax, scheme);
            if preconditioner_zeta_remove_cyclic
                ddzeta_preconditioner(1,end) = 0;
                ddzeta_preconditioner(end,1) = 0;                
            end
        end
        
        % Evaluate the magnetic field and its derivatives on the
        % (theta,zeta) grid:
        computeBHat()
        
        % Compute a few quantities related to the magnetic field
        VPrimeHat = thetaWeights' * (1./BHat.^2) * zetaWeights;
        FSABHat2 = 4*pi*pi/VPrimeHat;
        
        % Generate abscissae, quadrature weights, and derivative matrices for
        % the energy (x) grid used to represent the distribution function.
        k=0;
        scale=1;
        pointAtZero = false;
        [x, ddx, d2dx2, xWeights] = m20130312_02_SpectralNodesWeightsAndDifferentiationMatricesForV(Nx, k, scale, pointAtZero);
        
        % Make the energy grid and differentiation matrices for the
        % Rosenbluth potentials:
        function y=weight(x)
            x2temp=x.*x;
            y=exp(-x2temp);
        end
        xMaxPotentials=max([5, max(x)]);
        xMin=0;
        NxPotentials = ceil(xMaxPotentials * NxPotentialsPerVth);
        % Uniform grid with 5-point stencil for derivatives:
        scheme = 12;
        [xPotentials, ~, ddxPotentials, d2dx2Potentials] = m20121125_04_DifferentiationMatricesForUniformGrid(NxPotentials, xMin, xMaxPotentials, scheme);
        
        % Make the matrices for interpolating between the two energy grids:
        regridPolynomialToUniform = m20120703_03_polynomialInterpolationMatrix(x,xPotentials,weight(x),weight(xPotentials));
        %regridUniformToPolynomial = m20121127_02_makeHighOrderInterpolationMatrix(xPotentials,x,0,'f');
        
        
        if plotSpeedGrid
            if iteration == 1
                speedGridFigureHandle = figure(figureOffset+7);
            else
                set(0, 'CurrentFigure', speedGridFigureHandle);
            end
            plot(xPotentials,zeros(size(xPotentials))+iteration,'.r')
            hold on
            plot(x, zeros(size(x))+iteration,'o')
            title('Speed grid for distribution function (blue) and Rosenbluth potentials(red)')
            xlabel('x')
            ylabel('iteration')
        end
        
        % Set the size of the main linear system:
        matrixSize = Nx * Nxi * Ntheta * Nzeta;
        switch constraintScheme
            case 0
                % Nothing to do here.
            case 1
                matrixSize = matrixSize + 2;
            case 2
                matrixSize = matrixSize + Nx;
            otherwise
                error('Invalid setting for constraintScheme')
        end
        matrixSize = matrixSize * Nspecies;
        
        % To build the matrix as efficiently as possible, a reasonably
        % accurate estimate of the number of nonzeros (nnz) is needed beforehand:
        estimated_nnz = 1 * (Nx*Nx*Nspecies*Nspecies*Nxi*Ntheta*Nzeta ...
            + Nspecies*(nnz(ddtheta)*Nx*(3*Nxi)*Nzeta + nnz(ddzeta)*Nx*(3*Nxi)*Ntheta ...
            + Nx*(5*Nxi)*Ntheta*Nzeta + 3*Nx*Nx*Nxi*Ntheta*Nzeta ...
            + 2*2*Nx*1*Ntheta*Nzeta));
        
        estimated_nnz_original = estimated_nnz;
        fprintf('matrixSize: %d.\n',matrixSize)
        
        if iteration==1
            figure(figureOffset+4);
            clf
            numRows=1+Nspecies;
            numCols=3;
            plotNum=1;
            numContours=15;
            
            if smoothFigures
                %[thetaFine, zetaFine] = meshgrid(linspace(0,max(theta)), linspace(0,max(zeta)));
                [zetaFine, thetaFine] = meshgrid(linspace(0,max(zeta)), linspace(0,max(theta)));
                thetaToUse = thetaFine;
                zetaToUse = zetaFine;
                BHatToUse = interp2(zeta,theta,BHat,zetaFine,thetaFine,'cubic');
                dBHatdthetaToUse = interp2(zeta,theta,dBHatdtheta,zetaFine,thetaFine,'cubic');
                dBHatdzetaToUse = interp2(zeta,theta,dBHatdzeta,zetaFine,thetaFine,'cubic');
            else
                thetaToUse = theta;
                zetaToUse = zeta;
                BHatToUse = BHat;
                dBHatdzetaToUse = dBHatdzeta;
                dBHatdthetaToUse = dBHatdtheta;
            end
            
            subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
            contourf(zetaToUse,thetaToUse,BHatToUse,numContours,'EdgeColor','none')
            colorbar
            xlabel('\zeta')
            ylabel('\theta')
            title('BHat')
            
            subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
            contourf(zetaToUse,thetaToUse,dBHatdthetaToUse,numContours,'EdgeColor','none')
            colorbar
            xlabel('\zeta')
            ylabel('\theta')
            title('dBHatdtheta (line indicates B direction)')
            hold on
            if iota>0
                plot([0,max(zeta)], [0, iota*max(zeta)],'k')
            else
                plot([0,max(zeta)], [max(theta), max(theta)+iota*max(zeta)],'k')
            end
            
            subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
            contourf(zetaToUse,thetaToUse,dBHatdzetaToUse,numContours,'EdgeColor','none')
            colorbar
            xlabel('\zeta')
            ylabel('\theta')
            title('dBHatdzeta (dots indicate numerical grid)')
            
            %[theta2D, zeta2D] = meshgrid(theta, zeta);
            [zeta2D, theta2D] = meshgrid(zeta, theta);
            hold on
            plot(zeta2D, theta2D, '.k')
            
            drawnow
        end
        
        if RHSMode == 2
            error('RHSMode = 2 is not correctly implemented yet!')
            nu_n = nuPrime * sqrt(THat) * B0OverBBar / (GHat+iota*IHat);
        end
        
        % Begin timer for matrix construction:
        tic
        
        % ------------------------------------------------------
        % ------------------------------------------------------
        % Build the right-hand side of the main linear system.
        % ------------------------------------------------------
        % ------------------------------------------------------
        
        switch RHSMode
            case 1
                RHSSize = 1;
            case 2
                RHSSize = 3;
            otherwise
                error('Invalid RHSMode')
        end
        
        x2=x.*x;
        expx2=exp(-x2);
        rhs=zeros(matrixSize,RHSSize);
        
        for ispecies = 1:Nspecies
            spatialPartOfRHS_gradients = Delta * nHats(ispecies) * THats(ispecies) ...
                .* ((mHats(ispecies)/THats(ispecies))^(3/2)) ...
                ./ (factorToIncludeInFNormalization*2*pi*sqrtpi*psiAHat * Zs(ispecies) * (BHat.^3)) ...
                .* (GHat*dBHatdtheta - IHat*dBHatdzeta) ;
            
            spatialPartOfRHS_EPar = alpha*Zs(ispecies)*nHats(ispecies)*mHats(ispecies) ...
                * (GHat+iota*IHat)./(factorToIncludeInFNormalization*pi*sqrtpi*THats(ispecies)^2*FSABHat2*BHat);
            
            for col=1:RHSSize
                switch RHSMode
                    case 1
                        dnHatdpsiToUse = dnHatdpsiNs(ispecies);
                        dTHatdpsiToUse = dTHatdpsiNs(ispecies);
                        dPhiHatdpsiToUse = dPhiHatdpsiN;
                        EParallelHatToUse = EParallelHat;
                    case 2
                        dPhiHatdpsiToUse = 0;
                        switch col
                            case 1
                                dnHatdpsiToUse = 1;
                                dTHatdpsiToUse = 0;
                                EParallelHatToUse = 0;
                            case 2
                                % The next 2 lines ensure (1/n)*dn/dpsi + (3/2)*dT/dpsi = 0 while dT/dpsi is nonzero.
                                dnHatdpsiToUse = (3/2)*nHat/THat;
                                dTHatdpsiToUse = 1;
                                EParallelHatToUse = 0;
                            case 3
                                dnHatdpsiToUse = 0;
                                dTHatdpsiToUse = 0;
                                EParallelHatToUse = 1;
                        end
                    otherwise
                        error('Invalid RHSMode')
                end
                
                xPartOfRHS_gradients = x2.*expx2.*(dnHatdpsiToUse/nHats(ispecies) + alpha*Zs(ispecies)*dPhiHatdpsiToUse/THats(ispecies) + (x2-3/2)*dTHatdpsiToUse/THats(ispecies));
                xPartOfRHS_EPar = EParallelHatToUse*x.*expx2;
                
                for ix=1:Nx
                    for itheta=1:Ntheta
                        L=0;
                        indices = getIndices(ispecies, ix, L+1, itheta, 1:Nzeta, 0);
                        rhs(indices, col) = (4/3) * spatialPartOfRHS_gradients(itheta,:)' * xPartOfRHS_gradients(ix);
                        
                        L=1;
                        indices = getIndices(ispecies, ix, L+1, itheta, 1:Nzeta, 0);
                        rhs(indices, col) = spatialPartOfRHS_EPar(itheta,:)' * xPartOfRHS_EPar(ix);
                        
                        L=2;
                        indices = getIndices(ispecies, ix, L+1, itheta, 1:Nzeta, 0);
                        rhs(indices, col) = (2/3) * spatialPartOfRHS_gradients(itheta,:)' * xPartOfRHS_gradients(ix);
                    end
                end
            end
        end
        
        assignin('base','rhsm',rhs)
        
        % ------------------------------------------------------
        % ------------------------------------------------------
        % Build the matrix for the main linear system.
        % ------------------------------------------------------
        % ------------------------------------------------------
        
        sparseCreatorIndex=1;
        sparseCreator_i=0;
        sparseCreator_j=0;
        sparseCreator_s=0;
        resetSparseCreator()
        
        if tryIterativeSolver
            matricesToMake=0:1;
        else
            matricesToMake=1;
        end
        
        for whichMatrixToMake = matricesToMake
            % 0 = preconditioner
            % 1 = main matrix
            
            if whichMatrixToMake==1
                ddthetaToUse = ddtheta;
                ddzetaToUse = ddzeta;
                maxLForThetaDot = Nxi-1;
                maxLForZetaDot = Nxi-1;
                maxLForXiDot = Nxi-1;
            else
                ddthetaToUse = ddtheta_preconditioner;
                ddzetaToUse = ddzeta_preconditioner;
                maxLForThetaDot = min([preconditioner_theta_max_L, Nxi-1]);
                maxLForZetaDot = min([preconditioner_zeta_max_L, Nxi-1]);
                maxLForXiDot = min([preconditioner_xi_max_L, Nxi-1]);
            end
            
            matrixStartTime = tic;
            
            % -----------------------------------------
            % Add collisionless terms:
            % -----------------------------------------
            
            for ispecies = 1:Nspecies
                sqrtTHat = sqrt(THats(ispecies));
                sqrtMass = sqrt(mHats(ispecies));
                
                % -----------------------------------------
                % Add d/dtheta terms:
                % -----------------------------------------
                
                for izeta=1:Nzeta
                    if useDKESExBDrift
                        thetaPartOfExBTerm_lowL = alpha*Delta*GHat*dPhiHatdpsiN/(2*psiAHat*FSABHat2) * ddtheta;
                        thetaPartOfExBTerm_highL = alpha*Delta*GHat*dPhiHatdpsiN/(2*psiAHat*FSABHat2) * ddthetaToUse;
                    else
                        thetaPartOfExBTerm_lowL = alpha*Delta*GHat*dPhiHatdpsiN/(2*psiAHat) * diag(1./BHat(:,izeta).^2)*ddtheta;
                        thetaPartOfExBTerm_highL = alpha*Delta*GHat*dPhiHatdpsiN/(2*psiAHat) * diag(1./BHat(:,izeta).^2)*ddthetaToUse;
                    end
                    thetaPartOfStreamingTerm_lowL = iota*sqrtTHat/sqrtMass*diag(1./BHat(:,izeta))*ddtheta;
                    thetaPartOfStreamingTerm_highL = iota*sqrtTHat/sqrtMass*diag(1./BHat(:,izeta))*ddthetaToUse;
                    for L=0:maxLForThetaDot
                        if L < preconditioner_theta_min_L
                            thetaPartOfStreamingTerm = thetaPartOfStreamingTerm_lowL;
                            thetaPartOfExBTerm = thetaPartOfExBTerm_lowL;
                        else
                            thetaPartOfStreamingTerm = thetaPartOfStreamingTerm_highL;
                            thetaPartOfExBTerm = thetaPartOfExBTerm_highL;
                        end
                        
                        for ix=1:Nx
                            rowIndices = getIndices(ispecies, ix, L+1, 1:Ntheta, izeta, 0);
                            
                            % Diagonal term
                            addSparseBlock(rowIndices, rowIndices, thetaPartOfExBTerm)
                            
                            % Super-diagonal term
                            if (L<maxLForThetaDot)
                                colIndices = getIndices(ispecies, ix, L+1+1, 1:Ntheta, izeta, 0);
                                addSparseBlock(rowIndices, colIndices, x(ix)*(L+1)/(2*L+3)*thetaPartOfStreamingTerm)
                            end
                            
                            % Sub-diagonal term
                            if (L>0)
                                colIndices = getIndices(ispecies, ix, L+1-1, 1:Ntheta, izeta, 0);
                                addSparseBlock(rowIndices, colIndices, x(ix)*L/(2*L-1)*thetaPartOfStreamingTerm)
                            end
                            
                        end
                    end
                end
                
                % -----------------------------------------
                % Add d/dzeta terms:
                % -----------------------------------------
                
                for itheta=1:Ntheta
                    if useDKESExBDrift
                        zetaPartOfExBTerm_lowL = -alpha*Delta*IHat*dPhiHatdpsiN/(2*psiAHat*FSABHat2) *ddzeta;
                        zetaPartOfExBTerm_highL = -alpha*Delta*IHat*dPhiHatdpsiN/(2*psiAHat*FSABHat2) *ddzetaToUse;
                    else
                        zetaPartOfExBTerm_lowL = -alpha*Delta*IHat*dPhiHatdpsiN/(2*psiAHat) * diag(1./BHat(itheta,:).^2)*ddzeta;
                        zetaPartOfExBTerm_highL = -alpha*Delta*IHat*dPhiHatdpsiN/(2*psiAHat) * diag(1./BHat(itheta,:).^2)*ddzetaToUse;
                    end
                    zetaPartOfStreamingTerm_lowL = sqrtTHat/sqrtMass*diag(1./BHat(itheta,:))*ddzeta;
                    zetaPartOfStreamingTerm_highL = sqrtTHat/sqrtMass*diag(1./BHat(itheta,:))*ddzetaToUse;
                    for L=0:maxLForZetaDot
                        if L < preconditioner_zeta_min_L
                            zetaPartOfExBTerm = zetaPartOfExBTerm_lowL;
                            zetaPartOfStreamingTerm = zetaPartOfStreamingTerm_lowL;
                        else
                            zetaPartOfExBTerm = zetaPartOfExBTerm_highL;
                            zetaPartOfStreamingTerm = zetaPartOfStreamingTerm_highL;
                        end
                        
                        for ix=1:Nx
                            rowIndices = getIndices(ispecies, ix, L+1, itheta, 1:Nzeta, 0);
                            
                            % Diagonal term
                            addSparseBlock(rowIndices, rowIndices, zetaPartOfExBTerm)
                            
                            % Super-diagonal term
                            if (L<maxLForZetaDot)
                                colIndices = getIndices(ispecies, ix, L+1+1, itheta, 1:Nzeta, 0);
                                addSparseBlock(rowIndices, colIndices, x(ix)*(L+1)/(2*L+3)*zetaPartOfStreamingTerm)
                            end
                            
                            % Sub-diagonal term
                            if (L>0)
                                colIndices = getIndices(ispecies, ix, L-1+1, itheta, 1:Nzeta, 0);
                                addSparseBlock(rowIndices, colIndices, x(ix)*L/(2*L-1)*zetaPartOfStreamingTerm)
                            end
                            
                        end
                    end
                end
                
                
                % -----------------------------------------
                % Add d/dxi terms:
                % -----------------------------------------
                
                for itheta=1:Ntheta
                    spatialPartOfOldMirrorTerm = -sqrtTHat*(iota*dBHatdtheta(itheta,:)+dBHatdzeta(itheta,:))./(2*sqrtMass*BHat(itheta,:).^2);
                    spatialPartOfNewMirrorTerm = alpha*Delta*dPhiHatdpsiN*(GHat*dBHatdtheta(itheta,:) - IHat*dBHatdzeta(itheta,:))./(4*psiAHat*BHat(itheta,:).^3);
                    for ix=1:Nx
                        for L=0:maxLForXiDot
                            rowIndices = getIndices(ispecies, ix, L+1, itheta, 1:Nzeta, 0);
                            
                            % Super-diagonal term
                            if (L<maxLForXiDot)
                                colIndices = getIndices(ispecies, ix, L+1+1, itheta, 1:Nzeta, 0);
                                addToSparse(rowIndices, colIndices, x(ix)*(L+1)*(L+2)/(2*L+3)*spatialPartOfOldMirrorTerm)
                            end
                            
                            % Sub-diagonal term
                            if (L>0)
                                colIndices = getIndices(ispecies, ix, L-1+1, itheta, 1:Nzeta, 0);
                                addToSparse(rowIndices, colIndices, x(ix)*(-L)*(L-1)/(2*L-1)*spatialPartOfOldMirrorTerm)
                            end
                            
                            if includeElectricFieldTermInXiDot
                                % Diagonal term
                                addToSparse(rowIndices, rowIndices, L*(L+1)/((2*L-1)*(2*L+3))*spatialPartOfNewMirrorTerm)
                                
                                if (whichMatrixToMake==1 || preconditioner_xi==0)
                                    % Super-super-diagonal term:
                                    if (L < maxLForXiDot-1)
                                        colIndices = getIndices(ispecies, ix, L+2+1, itheta, 1:Nzeta, 0);
                                        addToSparse(rowIndices, colIndices, (L+1)*(L+2)*(L+3)/((2*L+5)*(2*L+3))*spatialPartOfNewMirrorTerm)
                                    end
                                    
                                    % Sub-sub-diagonal term:
                                    if (L > 1)
                                        colIndices = getIndices(ispecies, ix, L-2+1, itheta, 1:Nzeta, 0);
                                        addToSparse(rowIndices, colIndices, -L*(L-1)*(L-2)/((2*L-3)*(2*L-1))*spatialPartOfNewMirrorTerm)
                                    end
                                end
                            end
                        end
                    end
                end
                
                
                % -----------------------------------------
                % Add the collisionless d/dx term:
                % -----------------------------------------
                
                if includeXDotTerm
                    xPartOfXDot = diag(x)*ddx;
                    if (whichMatrixToMake==1)
                        xPartOfXDotForLargeL = xPartOfXDot;
                    else
                        % We're making the preconditioner, so simplify matrix
                        % if needed:
                        switch preconditioner_x
                            case 0
                                xPartOfXDotForLargeL = xPartOfXDot;
                            case 1
                                xPartOfXDotForLargeL = diag(diag(xPartOfXDot));
                            case 2
                                xPartOfXDotForLargeL = triu(xPartOfXDot);
                            case 3
                                mask = eye(Nx) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);
                                xPartOfXDotForLargeL = xPartOfXDot .* mask;
                            case 4
                                mask = eye(Nx) + diag(ones(Nx-1,1),1);
                                xPartOfXDotForLargeL = xPartOfXDot .* mask;
                            otherwise
                                error('Invalid setting for preconditioner_x')
                        end
                    end
                    for L=0:(Nxi-1)
                        if L >= preconditioner_x_min_L
                            xPartOfXDotToUse = xPartOfXDotForLargeL;
                        else
                            xPartOfXDotToUse = xPartOfXDot;
                        end
                        for itheta=1:Ntheta
                            for izeta=1:Nzeta
                                spatialPart = alpha*Delta*dPhiHatdpsiN*(GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))/(4*psiAHat*BHat(itheta,izeta)^3);
                                
                                rowIndices = getIndices(ispecies, 1:Nx, L+1, itheta, izeta, 0);
                                
                                % Diagonal term
                                addSparseBlock(rowIndices, rowIndices, 2*(3*L*L+3*L-2)/((2*L+3)*(2*L-1))*spatialPart*xPartOfXDotToUse)
                                
                                if (whichMatrixToMake==1 || preconditioner_xi==0)
                                    % Super-super-diagonal in L
                                    if (L<Nxi-2)
                                        colIndices = getIndices(ispecies, 1:Nx, L+2+1, itheta, izeta, 0);
                                        addSparseBlock(rowIndices, colIndices, (L+1)*(L+2)/((2*L+5)*(2*L+3))*spatialPart*xPartOfXDotToUse)
                                    end
                                    
                                    % Sub-sub-diagonal in L
                                    if (L>1)
                                        colIndices = getIndices(ispecies, 1:Nx, L-2+1, itheta, izeta, 0);
                                        addSparseBlock(rowIndices, colIndices, L*(L-1)/((2*L-3)*(2*L-1))*spatialPart*xPartOfXDotToUse)
                                    end
                                    
                                end
                            end
                        end
                    end
                end
                
                
                % -----------------------------------------
                % Add the optional f div dot v_E term.
                % This term may make sense with the partial trajectories since it
                % restores particle and energy conservation.
                % -----------------------------------------
                
                if include_fDivVE_term
                    for itheta=1:Ntheta
                        for izeta=1:Nzeta
                            elementsToAdd = -dPhiHatdpsiN*Delta*alpha/(psiAHat*(BHat(itheta,izeta)^3))...
                                *(GHat*dBHatdtheta(itheta,izeta) - IHat*dBHatdzeta(itheta,izeta))...
                                *ones(1,Nx);
                            for L=0:(Nxi-1)
                                indices = getIndices(ispecies, 1:Nx, L+1, itheta, izeta, 0);
                                addToSparse(indices, indices, elementsToAdd)
                            end
                        end
                    end
                end
                
                
                
            end
            
            % -----------------------------------------
            % Add collision operator:
            % -----------------------------------------
            
            %{
            % Beginning of original 1-species Fokker-Planck operator.
            erfs=erf(x);
            x2 = x.*x;
            x3 = x2.*x;
            expx2 = exp(-x.*x);
            % Psi is the Chandrasekhar function:
            Psi = (erfs - 2/sqrtpi*x .* expx2) ./ (2*x.*x);
            % Energy-dependent deflection frequency:
            nuD = 3*sqrtpi/4*(erfs - Psi) ./ x3;
            
            switch collisionOperator
                case 0
                    % Full linearized Fokker-Planck operator
                    
                    xWith0s = [0; xPotentials(2:(end-1)); 0];
                    M21 = 4*pi*diag(xWith0s.^2) * regridPolynomialToUniform;
                    M32 = -2*diag(xWith0s.^2);
                    LaplacianTimesX2WithoutL = diag(xPotentials.^2)*d2dx2Potentials + 2*diag(xPotentials)*ddxPotentials;
                    
                    PsiPrime = (-erfs + 2/sqrtpi*x.*(1+x.*x) .* expx2) ./ x3;
                    xPartOfCECD = 3*sqrtpi/4*(diag(Psi./x)*d2dx2   +  diag((PsiPrime.*x  + Psi + 2*Psi.*x2)./x2)*ddx + diag(2*PsiPrime + 4*Psi./x)) + 3*diag(expx2);
                    M12IncludingX0 = nu_n * 3/(2*pi)*diag(expx2)*regridUniformToPolynomial;
                    M13IncludingX0 = -nu_n * 3/(2*pi) * diag(x2.*expx2) * regridUniformToPolynomial* d2dx2Potentials;
                    
                    for L=0:(Nxi-1)
                        M11 = -nu_n * (-0.5*diag(nuD)*L*(L+1) + xPartOfCECD);
                        
                        if L <= (NL-1)
                            % Add Rosenbluth potential stuff
                            
                            M13 = M13IncludingX0;
                            M12 = M12IncludingX0;
                            
                            M22 = LaplacianTimesX2WithoutL-L*(L+1)*eye(NxPotentials);
                            % Add Dirichlet or Neumann boundary condition for
                            % potentials at x=0:
                            if L==0
                                M22(1,:)=ddxPotentials(1,:);
                            else
                                M22(1,:) = 0;
                                M22(1,1) = 1;
                                M12(:,1) = 0;
                                M13(:,1) = 0;
                            end
                            M33 = M22;
                            
                            % Add Robin boundary condition for potentials at x=xMaxPotentials:
                            M22(NxPotentials,:) = xMaxPotentials*ddxPotentials(NxPotentials,:);
                            M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1;
                            
                            M33(NxPotentials,:) = xMaxPotentials*xMaxPotentials*d2dx2Potentials(NxPotentials,:) + (2*L+1)*xMaxPotentials*ddxPotentials(NxPotentials,:);
                            M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1);
                            if L~=0
                                M22(NxPotentials,1)=0;
                                M33(NxPotentials,1)=0;
                            end
                            
                            KWithoutThetaPart = M11 -  (M12 - M13 * (M33 \ M32)) * (M22 \ M21);
                            
                        else
                            KWithoutThetaPart = M11;
                        end
                        
                        if (whichMatrixToMake==0 && L >= preconditioner_x_min_L)
                            % We're making the preconditioner, so simplify
                            % the matrix if needed.
                            switch preconditioner_x
                                case 0
                                    % Do nothing.
                                case 1
                                    KWithoutThetaPart = diag(diag(KWithoutThetaPart));
                                case 2
                                    KWithoutThetaPart = triu(KWithoutThetaPart);
                                case 3
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);
                                    KWithoutThetaPart = KWithoutThetaPart .* mask;
                                case 4
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1);
                                    KWithoutThetaPart = KWithoutThetaPart .* mask;
                                otherwise
                                    error('Invalid setting for preconditioner_x')
                            end
                        end
                        
                        for itheta=1:Ntheta
                            for izeta=1:Nzeta
                                indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                                addSparseBlock(indices, indices, (GHat+iota*IHat)/(BHat(itheta,izeta)^2)*KWithoutThetaPart)
                            end
                        end
                    end
                    
                case {1,2}
                    % Pitch angle scattering operator
                    
                    for itheta=1:Ntheta
                        for izeta=1:Nzeta
                            spatialPart = -nu_n*(GHat+iota*IHat)/(BHat(itheta,izeta)^2);
                            for L=1:(Nxi-1)
                                indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                                addToSparse(indices, indices, -0.5*L*(L+1)*spatialPart*nuD)
                            end
                        end
                    end
                    
                    if collisionOperator==2
                        % Add model field term
                        
                        L=1;
                        fieldTerm = (nuD.*x.*expx2) * ((xWeights.*x.*x2.*nuD)')/0.354162849836926;
                        
                        if (whichMatrixToMake==0 && L >= preconditioner_x_min_L)
                            switch preconditioner_x
                                case 0
                                    % Nothing to do here.
                                case 1
                                    fieldTerm = diag(diag(fieldTerm));
                                case 2
                                    fieldTerm = triu(fieldTerm);
                                case 3
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);
                                    fieldTerm = fieldTerm .* mask;
                                case 4
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1);
                                    fieldTerm = fieldTerm .* mask;
                                otherwise
                                    error('Invalid setting for preconditioner_x')
                            end
                        end
                        
                        for itheta=1:Ntheta
                            for izeta=1:Nzeta
                                indices = ((1:Nx)-1)*Nxi*Ntheta*Nzeta + L*Ntheta*Nzeta + (itheta-1)*Nzeta + izeta;
                                addSparseBlock(indices, indices, -nu_n*(GHat+iota*IHat)/(BHat(itheta,izeta)^2)*fieldTerm)
                            end
                        end
                        
                    end
                otherwise
                    error('Invalid setting for collisionOperator')
            end
            % End of original 1-species collision operator
            %}
            
            % Start of new multi-species collision operators
            switch (collisionOperator)
              case 0
                % Linearized Fokker-Planck operator
                
                xWith0s = [0, xPotentials(2:(end-1))', 0];
                M21 = 4*pi*diag(xWith0s.^2) * regridPolynomialToUniform;
                xWith0s = [0, xPotentials(2:(end-1))', 0];
                M32 = -2*diag(xWith0s.^2);
                LaplacianTimesX2WithoutL = diag(xPotentials.^2)*d2dx2Potentials + 2*diag(xPotentials)*ddxPotentials;
                
                x2 = x.*x;
                expx2 = exp(-x.*x);
                
                CE = zeros(Nx, Nx, Nspecies);
                nuD = zeros(Nx, Nspecies);
                regridSpecies = zeros(Nx, Nx, Nspecies, Nspecies);
                M12IncludingX0 = zeros(Nx, NxPotentials, Nspecies, Nspecies, NL);
                M13IncludingX0 = zeros(Nx, NxPotentials, Nspecies, Nspecies, NL);
                for speciesA = 1:Nspecies
                    for speciesB = 1:Nspecies
                        speciesFactorTest = 3*sqrtpi/4*nHats(speciesB) * Zs(speciesA)*Zs(speciesA)*Zs(speciesB)*Zs(speciesB)/(THats(speciesA)^(3/2)*sqrt(mHats(speciesA)));
                        xb = x * sqrt(THats(speciesA)*mHats(speciesB)/(THats(speciesB)*mHats(speciesA)));
                        erfs = erf(xb);
                        xb2  = xb.*xb;
                        expxb2 = exp(-xb2);
                        Psi = (erfs - 2/sqrtpi*xb .* expxb2) ./ (2*xb2);
                        nuD(:,speciesA) = nuD(:,speciesA) + (speciesFactorTest * (erfs - Psi) ./ (x.^3));
                        coefficientOfd2dx2 = Psi./x;
                        coefficientOfddx = -2*THats(speciesA)*mHats(speciesB)/(THats(speciesB)*mHats(speciesA))*Psi*(1-mHats(speciesA)/mHats(speciesB)) ...
                            + (erfs - Psi)./(x.*x);
                        diagonalPartOfCE = 4/sqrtpi*THats(speciesA)/THats(speciesB)*sqrt(THats(speciesA)*mHats(speciesB)/(THats(speciesB)*mHats(speciesA))) .* expxb2;
                        CE(:,:,speciesA) = CE(:,:,speciesA) + speciesFactorTest*(diag(coefficientOfd2dx2)*d2dx2 + diag(coefficientOfddx)*ddx + diag(diagonalPartOfCE));
                        
                        if speciesA==speciesB
                            regridSpecies(:,:,speciesA,speciesB) = eye(Nx);
                        else
                            regridSpecies(:,:,speciesA,speciesB) = m20120703_03_polynomialInterpolationMatrix(x,xb,weight(x),weight(xb));
                        end
                        
                        speciesFactorField = nHats(speciesA) * Zs(speciesA)*Zs(speciesA)*Zs(speciesB)*Zs(speciesB)...
                            * mHats(speciesA) * THats(speciesB)/(THats(speciesA)^(5/2) * mHats(speciesB) * sqrt(mHats(speciesA)));
                        for L=0:(NL-1)
                            regridUniformToPolynomial = m20120925_09_makeHighOrderUniformRegriddingMatrix(xPotentials,xb,L,'H');
                            M12IncludingX0(:,:,speciesA, speciesB, L+1) = -3/(2*pi)*speciesFactorField*diag(expx2)* regridUniformToPolynomial...
                                * (diag(xPotentials*(1-mHats(speciesA)/mHats(speciesB)))*ddxPotentials + eye(NxPotentials)) ;
                            regridUniformToPolynomial = m20120925_09_makeHighOrderUniformRegriddingMatrix(xPotentials,xb,L,'G');
                            M13IncludingX0(:,:,speciesA, speciesB, L+1) = 3/(2*pi) * speciesFactorField * diag(x2.*expx2) * regridUniformToPolynomial* d2dx2Potentials;
                        end
                    end
                end
                
                for L=0:(Nxi-1)
                    if L <= (NL-1)
                        % Add Rosenbluth potential stuff
                        
                        M22 = LaplacianTimesX2WithoutL-L*(L+1)*eye(NxPotentials);
                        % Add Dirichlet or Neumann boundary condition for
                        % potentials at x=0:
                        if L==0
                            M22(1,:)=ddxPotentials(1,:);
                        else
                            M22(1,:) = 0;
                            M22(1,1) = 1;
                        end
                        M33 = M22;
                        
                        % Add Robin boundary condition for potentials at x=xMaxPotentials:
                        M22(NxPotentials,:) = xMaxPotentials*ddxPotentials(NxPotentials,:);
                        M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1;
                        
                        % Boundary conditions:
                        M33(NxPotentials,:) = xMaxPotentials*xMaxPotentials*d2dx2Potentials(NxPotentials,:) + (2*L+1)*xMaxPotentials*ddxPotentials(NxPotentials,:);
                        M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1);
                        
                        if L~=0
                            M22(NxPotentials,1)=0;
                            M33(NxPotentials,1)=0;
                        end
                        
                        M22BackslashM21 = M22 \ M21;
                        M33BackslashM32 = M33 \ M32;
                        
                    end
                    
                    for speciesA = 1:Nspecies
                        if whichMatrixToMake == 1
                            % We're making the main matrix.
                            speciesBToUse = 1:Nspecies;
                        else
                            % We're making the preconditioner.
                            switch preconditioner_species
                              case 0
                                % Full inter-species coupling
                                speciesBToUse = 1:Nspecies;
                              case 1
                                % No inter-species coupling
                                speciesBToUse = speciesA;
                              otherwise
                                error('Invalid preconditioner_species')
                            end
                        end
                        for speciesB = speciesBToUse
                            % Add CD
                            CD = 3*nHats(speciesA)*Zs(speciesA)*Zs(speciesA)*Zs(speciesB)*Zs(speciesB)...
                                 * mHats(speciesA)/(mHats(speciesB)*THats(speciesA)*sqrt(THats(speciesA)*mHats(speciesA))) ...
                                 * diag(expx2) * regridSpecies(:,:,speciesA, speciesB);
                            
                            if speciesA == speciesB
                                M11 = -0.5*diag(nuD(:,speciesA))*L*(L+1) + CE(:,:,speciesA) + CD;
                            else
                                M11 = CD;
                            end
                            
                            if L <= (NL-1)
                                % Add terms of the collision operator involving
                                % the Rosenbluth potentials.
                                
                                M13 = M13IncludingX0(:,:,speciesA, speciesB, L+1);
                                M12 = M12IncludingX0(:,:,speciesA, speciesB, L+1);
                                
                                % Add Dirichlet or Neumann boundary condition for
                                % potentials at x=0:
                                if L~=0
                                    M12(:,1) = 0;
                                    M13(:,1) = 0;
                                end
                                
                                CHat = M11 -  (M12 - M13 * M33BackslashM32) * M22BackslashM21;
                            else
                                CHat = M11;
                            end
                            
                            % The lines below are invoked to make the local preconditioner.
                            if whichMatrixToMake == 0 && L >= preconditioner_x_min_L
                                switch preconditioner_x
                                  case 0
                                    % Nothing to do here.
                                  case 1
                                    CHat = diag(diag(CHat));
                                  case 2
                                    CHat = triu(CHat);
                                  case 3
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);
                                    CHat = CHat .* mask;
                                  case 4
                                    mask = eye(Nx) + diag(ones(Nx-1,1),1);
                                    CHat = CHat .* mask;
                                  otherwise
                                    error('Invalid preconditionerMethod_x')
                                end
                                
                            end
                            
                            % At this point, CHat holds the collision operator
                            % divided by \bar{nu}
                            
                            for itheta = 1:Ntheta
                                for izeta = 1:Nzeta
                                    rowIndices = getIndices(speciesA, 1:Nx, L+1, itheta, izeta, 0);
                                    colIndices = getIndices(speciesB, 1:Nx, L+1, itheta, izeta, 0);
                                    addSparseBlock(rowIndices, colIndices, -nu_n*(GHat+iota*IHat)/(BHat(itheta,izeta)^2)*CHat)
                                end
                            end
                            
                            
                        end
                    end
                    
                end
                % End of new multi-species Fokker-Planck collision
                % operator.
                
              case (1)
                % Pure pitch angle scattering collision operator
            
                % First, assemble the deflection frequency nuD for
                % species A, which involves a sum over species B:
                nuD = zeros(Nx, Nspecies);
                for speciesA = 1:Nspecies
                    for speciesB = 1:Nspecies
                        speciesFactorTest = 3*sqrtpi/4*nHats(speciesB) * Zs(speciesA)*Zs(speciesA)*Zs(speciesB)*Zs(speciesB)/(THats(speciesA)^(3/2)*sqrt(mHats(speciesA)));
                        xb = x * sqrt(THats(speciesA)*mHats(speciesB)/(THats(speciesB)*mHats(speciesA)));
                        erfs = erf(xb);
                        xb2  = xb.*xb;
                        expxb2 = exp(-xb2);
                        Psi = (erfs - 2/sqrtpi*xb .* expxb2) ./ (2*xb2);
                        nuD(:,speciesA) = nuD(:,speciesA) + (speciesFactorTest * (erfs - Psi) ./ (x.^3));
                    end
                end
                
                % Now that nuD has been assembled, 
                for L=0:(Nxi-1)                    
                    for iSpecies = 1:Nspecies
                        CHat = -0.5*diag(nuD(:,iSpecies))*L*(L+1);
                                                        
                        % At this point, CHat holds the collision operator
                        % divided by \bar{nu}
                            
                        for itheta = 1:Ntheta
                            for izeta = 1:Nzeta
                                indices = getIndices(iSpecies, 1:Nx, L+1, itheta, izeta, 0);
                                addToSparse(indices, indices, -nu_n*(GHat+iota*IHat)/(BHat(itheta,izeta)^2)*CHat)
                            end
                        end
                    end
                end
                
                % End of new multi-species pitch-angle scattering collision
                % operator.
              
              otherwise
                error('collisionOperator must be 0 or 1.')
            end
            
            % --------------------------------------------------
            % Add constraints.
            % --------------------------------------------------
            
            switch constraintScheme
                case 0
                    % Do nothing.
                    
                case 1
                    
                    L=0;
                    for ispecies = 1:Nspecies
                        for itheta=1:Ntheta
                            for izeta=1:Nzeta
                                colIndices = getIndices(ispecies, 1:Nx, L+1, itheta, izeta, 0);
                                
                                rowIndex = getIndices(ispecies, 1, 1, 1, 1, 1);
                                addSparseBlock(rowIndex, colIndices, (x2.*xWeights)' / (BHat(itheta,izeta)^2))
                                
                                rowIndex = getIndices(ispecies, 1, 1, 1, 1, 2);
                                addSparseBlock(rowIndex, colIndices, (x2.*x2.*xWeights)' / (BHat(itheta,izeta)^2))
                            end
                        end
                    end
                    
                case 2
                    L=0;
                    spatialPart = 1./(BHat.*BHat);
                    %spatialPart = reshape(spatialPart',[Ntheta*Nzeta,1])';
                    for ispecies = 1:Nspecies
                        for ix=1:Nx
                            rowIndex = getIndices(ispecies, ix, 1, 1, 1, 3);
                            for itheta=1:Ntheta
                                colIndices = getIndices(ispecies, ix, L+1, itheta, 1:Nzeta, 0);
                                addSparseBlock(rowIndex, colIndices, spatialPart(itheta,:))
                                %rowIndex = Nx*Nxi*Ntheta*Nzeta + ix;
                                %colIndices = (ix-1)*Nxi*Ntheta*Nzeta + (1:(Ntheta*Nzeta));
                                %addSparseBlock(rowIndex, colIndices, spatialPart)
                            end
                        end
                    end
                    
                otherwise
                    error('Invalid constraintScheme')
            end
            
            % --------------------------------------------------
            % Add sources.
            % --------------------------------------------------
            
            % Even though the source is independent of theta and
            % zeta, the kinetic equation is normalized by
            % multiplying through with 1/(B dot grad zeta) ~ 1/B^2
            spatialPart = 1./BHat.^2;

            switch constraintScheme
                case 0
                    % Do nothing
                    
                case 1
                    xPartOfSource1 = (x2-5/2).*expx2;
                    xPartOfSource2 = (x2-3/2).*expx2;
                    
                    L=0;
                    for ispecies = 1:Nspecies
                        for ix=1:Nx
                            for itheta=1:Ntheta
                                rowIndices = getIndices(ispecies, ix, L+1, itheta, 1:Nzeta, 0);
                                
                                colIndex = getIndices(ispecies, 1, 1, 1, 1, 1);
                                addSparseBlock(rowIndices, colIndex, xPartOfSource1(ix)*spatialPart(itheta,:)')
                                
                                colIndex = getIndices(ispecies, 1, 1, 1, 1, 2);
                                addSparseBlock(rowIndices, colIndex, xPartOfSource2(ix)*spatialPart(itheta,:)')
                            end
                        end
                    end
                    
                case 2
                    L=0;
                    for ispecies = 1:Nspecies
                        for ix=1:Nx
                            colIndex = getIndices(ispecies, ix, 1, 1, 1, 3);
                            for itheta=1:Ntheta
                                rowIndices = getIndices(ispecies, ix, L+1, itheta, 1:Nzeta, 0);
                                addSparseBlock(rowIndices, colIndex, spatialPart(itheta,:)')
                            end
                        end
                    end
                otherwise
                    error('Invalid constraintScheme')
            end
            
            % --------------------------------------------------
            % End of adding entries to the matrix.
            % --------------------------------------------------
            
            switch whichMatrixToMake
                case 0
                    fprintf('Time to contruct preconditioner: %g seconds.\n',toc(matrixStartTime))
                    tic
                    preconditionerMatrix = createSparse();
                    fprintf('Time to sparsify preconditioner: %g seconds.\n',toc)
                case 1
                    fprintf('Time to contruct main matrix: %g seconds.\n',toc(matrixStartTime))
                    tic
                    matrix = createSparse();
                    fprintf('Time to sparsify main matrix: %g seconds.\n',toc)
                otherwise
                    error('Program should not get here')
            end
        end
        
        
        % ------------------------------------------------------
        % Finalize matrix
        % ------------------------------------------------------
        
        if ~ tryIterativeSolver
            preconditionerMatrix = matrix;
        end
        fprintf('Actual nnz of matrix: %d,  of preconditioner: %d,  predicted: %d\n',nnz(matrix),nnz(preconditionerMatrix),estimated_nnz_original)
        fprintf('Fraction of nonzero entries in matrix: %g,  in preconditioner: %g\n', nnz(matrix)/numel(matrix), nnz(preconditionerMatrix)/numel(preconditionerMatrix))
        fprintf('nnz(preconditioner)/nnz(matrix): %g\n', nnz(preconditionerMatrix)/nnz(matrix))
        
        assignin('base','mm',matrix)
        if tryIterativeSolver
            assignin('base','pm',preconditionerMatrix)
        end
        
        
        if tryIterativeSolver
            fprintf('LU-factorizing preconditioner...')
            tic
            [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(preconditionerMatrix);
            fprintf('done.  Took %g seconds.\n',toc)
        end
        
        
        function solnVector=preconditioner(rhsVector)
            solnVector = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * rhsVector)));
        end
        
        % ------------------------------------------------------
        % Solve the main linear system
        % ------------------------------------------------------
        
        [soln, didItConverge, residual] ...
            = m20130728_01_solveLinearSystem(matrix, rhs, @preconditioner, ...
            tryIterativeSolver, orderOfSolversToTry, tol, maxIterations, restart, ...
            figureOffset+9, tryDirectSolverIfIterativeSolversFail);        
        
        NTVkernel; %This is a dummy line which is only here to let the variable NTVkernel
                   %from the subroutine computeBHat() be known also in computeOutputs().       
        computeOutputs()
        
        % ------------------------------------------------------
        % Calculate radial heat and particle fluxes
        % ------------------------------------------------------
        function computeOutputs()
            
            for col = 1:RHSSize
                
                if RHSSize > 1
                    fprintf('--- Analyzing solution vector %d of %d. ---\n', col, RHSSize)
                end
                
                switch constraintScheme
                    case 0
                        % Do nothing
                    case 1
                        sources = zeros(Nspecies,2);
                        for ispecies = 1:Nspecies
                            source1 = soln(getIndices(ispecies,1,1,1,1,1), col);
                            source2 = soln(getIndices(ispecies,1,1,1,1,2), col);
                            sources(ispecies,:) = [source1,source2];
                            fprintf('Sources for species %d: %g,  %g\n',ispecies,source1,source2)
                        end
                    case 2
                        sources = zeros(Nspecies,Nx);
                        for ispecies = 1:Nspecies
                            sources(ispecies,:) = soln(getIndices(ispecies,1:Nx,1,1,1,3), col);
                            fprintf('Species %d: min source = %g,   max source = %g\n',min(sources(ispecies,:)),max(sources(ispecies,:)))
                        end
                    otherwise
                        error('Invalid constraintScheme')
                end
                
                
                densityPerturbation = zeros(Ntheta,Nzeta,Nspecies);
                flow = zeros(Ntheta,Nzeta,Nspecies);
                pressurePerturbation = zeros(Ntheta,Nzeta,Nspecies);
                
                particleFluxBeforeSurfaceIntegral = zeros(Ntheta,Nzeta,Nspecies);
                momentumFluxBeforeSurfaceIntegral = zeros(Ntheta,Nzeta,Nspecies);
                heatFluxBeforeSurfaceIntegral = zeros(Ntheta,Nzeta,Nspecies);
                NTVBeforeSurfaceIntegral = zeros(Ntheta,Nzeta,Nspecies);
                
                FSADensityPerturbation = zeros(Nspecies,1);
                FSABFlow = zeros(Nspecies,1);
                FSAPressurePerturbation = zeros(Nspecies,1);
                
                particleFlux = zeros(Nspecies,1);
                momentumFlux = zeros(Nspecies,1);
                heatFlux = zeros(Nspecies,1);
                NTV = zeros(Nspecies,1);
                
                jHat = zeros(Ntheta, Nzeta);
                FSABjHat = 0;
                
                densityPerturbationIntegralWeight = x.^2;
                flowIntegralWeight = x.^3;
                pressurePerturbationIntegralWeight = x.^4;
                
                particleFluxIntegralWeight = x.^4;
                momentumFluxIntegralWeight = x.^5;
                heatFluxIntegralWeight = x.^6;
                NTVIntegralWeight = x.^4;
                
                for ispecies = 1:Nspecies
                    for itheta=1:Ntheta
                        for izeta = 1:Nzeta
                            L=0;
                            indices = getIndices(ispecies, 1:Nx, L+1, itheta, izeta, 0);
                            fSlice = soln(indices, col);
                            densityPerturbation(itheta,izeta,ispecies) = xWeights' * (densityPerturbationIntegralWeight .* fSlice);
                            pressurePerturbation(itheta,izeta,ispecies) = xWeights' * (pressurePerturbationIntegralWeight .* fSlice);
                            particleFluxBeforeSurfaceIntegral(itheta,izeta,ispecies) = (8/3)*xWeights' * (particleFluxIntegralWeight .* fSlice);
                            heatFluxBeforeSurfaceIntegral(itheta,izeta,ispecies) = (8/3)*xWeights' * (heatFluxIntegralWeight .* fSlice);
                            
                            L=1;
                            indices = getIndices(ispecies, 1:Nx, L+1, itheta, izeta, 0);
                            fSlice = soln(indices, col);
                            flow(itheta,izeta,ispecies) = xWeights' * (flowIntegralWeight .* fSlice);
                            momentumFluxBeforeSurfaceIntegral(itheta,izeta,ispecies) = (16/15)*xWeights' * (momentumFluxIntegralWeight .* fSlice);
                            
                            L=2;
                            indices = getIndices(ispecies, 1:Nx, L+1, itheta, izeta, 0);
                            fSlice = soln(indices, col);
                            particleFluxBeforeSurfaceIntegral(itheta,izeta,ispecies) = particleFluxBeforeSurfaceIntegral(itheta,izeta,ispecies) ...
                                + (4/15)*xWeights' * (particleFluxIntegralWeight .* fSlice);
                            heatFluxBeforeSurfaceIntegral(itheta,izeta,ispecies) = heatFluxBeforeSurfaceIntegral(itheta,izeta,ispecies) ...
                                + (4/15)*xWeights' * (heatFluxIntegralWeight .* fSlice);
                            NTVBeforeSurfaceIntegral(itheta,izeta,ispecies) = xWeights' * (NTVIntegralWeight .* fSlice);
                            
                            L=3;
                            indices = getIndices(ispecies, 1:Nx, L+1, itheta, izeta, 0);
                            fSlice = soln(indices, col);
                            momentumFluxBeforeSurfaceIntegral(itheta,izeta,ispecies) = momentumFluxBeforeSurfaceIntegral(itheta,izeta,ispecies) ...
                                + (4/35)*xWeights' * (momentumFluxIntegralWeight .* fSlice);
                            
                        end
                    end
                    
                    densityPerturbation(:,:,ispecies) = factorToIncludeInFNormalization*4*pi*((THats(ispecies)/mHats(ispecies)).^(3/2))*densityPerturbation(:,:,ispecies);
                    flow(:,:,ispecies) = factorToIncludeInFNormalization*(4/3)*pi *((THats(ispecies)/mHats(ispecies)).^2)* flow(:,:,ispecies);
                    pressurePerturbation(:,:,ispecies) = factorToIncludeInFNormalization*(8/3)*pi*mHats(ispecies)*((THats(ispecies)/mHats(ispecies)).^(5/2))*pressurePerturbation(:,:,ispecies);
                    
                    particleFluxBeforeSurfaceIntegral(:,:,ispecies) = -factorToIncludeInFNormalization*pi*Delta*THats(ispecies)*((THats(ispecies)/mHats(ispecies)).^(3/2)) ...
                        ./ (Zs(ispecies)*VPrimeHat*(GHat+iota*IHat)) ...
                        *(GHat*dBHatdtheta-IHat*dBHatdzeta)./(BHat.^3) ...
                        .* particleFluxBeforeSurfaceIntegral(:,:,ispecies);
                    
                    momentumFluxBeforeSurfaceIntegral(:,:,ispecies) = -factorToIncludeInFNormalization*pi*Delta*THats(ispecies)*((THats(ispecies)/mHats(ispecies)).^(2)) ...
                        ./ (Zs(ispecies)*VPrimeHat*(GHat+iota*IHat)) ...
                        *mHats(ispecies)*(GHat*dBHatdtheta-IHat*dBHatdzeta)./(BHat.^3) ...
                        .* momentumFluxBeforeSurfaceIntegral(:,:,ispecies);
                    
                    heatFluxBeforeSurfaceIntegral(:,:,ispecies) = -factorToIncludeInFNormalization*pi*Delta*THats(ispecies)*mHats(ispecies)*((THats(ispecies)/mHats(ispecies)).^(5/2)) ...
                        ./ (2*Zs(ispecies)*VPrimeHat*(GHat+iota*IHat)) ...
                        *(GHat*dBHatdtheta-IHat*dBHatdzeta)./(BHat.^3) ...
                        .* heatFluxBeforeSurfaceIntegral(:,:,ispecies);
                    
                    NTVBeforeSurfaceIntegral(:,:,ispecies) = factorToIncludeInFNormalization * ...
                        2*pi * THats(ispecies) * ((THats(ispecies) / mHats(ispecies)).^(3/2)) ./ (VPrimeHat*(GHat+iota*IHat)) ...
                        .* NTVkernel .* NTVBeforeSurfaceIntegral(:,:,ispecies);

                    FSADensityPerturbation(ispecies) = (1/VPrimeHat) * thetaWeights' * (densityPerturbation(:,:,ispecies)./(BHat.^2)) * zetaWeights;
                    FSABFlow(ispecies) = (1/VPrimeHat) * thetaWeights' * (flow(:,:,ispecies)./BHat) * zetaWeights;
                    FSAPressurePerturbation(ispecies) = (1/VPrimeHat) * thetaWeights' * (pressurePerturbation(:,:,ispecies)./(BHat.^2)) * zetaWeights;
                    
                    particleFlux(ispecies) = thetaWeights' * particleFluxBeforeSurfaceIntegral(:,:,ispecies) * zetaWeights;
                    momentumFlux(ispecies) = thetaWeights' * momentumFluxBeforeSurfaceIntegral(:,:,ispecies) * zetaWeights;
                    heatFlux(ispecies) = thetaWeights' * heatFluxBeforeSurfaceIntegral(:,:,ispecies) * zetaWeights;
                    
                    NTV(ispecies) = thetaWeights' * NTVBeforeSurfaceIntegral(:,:,ispecies) * zetaWeights;
                    
                    jHat = jHat + Zs(ispecies)*flow(:,:,ispecies);
                    FSABjHat = FSABjHat + Zs(ispecies)*FSABFlow(ispecies);
                    
                    fprintf('Results for species %d:\n',ispecies)
                    fprintf('FSADensityPerturbation:  %g\n',FSADensityPerturbation(ispecies))
                    fprintf('FSABFlow:                 %g\n',FSABFlow(ispecies))
                    fprintf('FSAPressurePerturbation: %g\n',FSAPressurePerturbation(ispecies))
                    fprintf('NTV:                     %g\n',NTV(ispecies))
                    fprintf('particleFlux:            %g\n',particleFlux(ispecies))
                    fprintf('momentumFlux:            %g\n',momentumFlux(ispecies))
                    fprintf('heatFlux:                %g\n',heatFlux(ispecies))
                    
                    %{
                    fprintf('FSABFlow in the 1-species code should be     %g\n',FSABFlow*Zs(1)*psiAHat/(Delta*nHats(1)))
                    fprintf('particleFlux in the 1-species code should be %g\n',particleFlux*psiAHat*Zs(1)^2*(GHat+iota*IHat)*VPrimeHat/(Delta^2*nHats(1)*sqrt(mHats(1))))
                    fprintf('heatFlux in the 1-species code should be     %g\n',heatFlux*psiAHat*Zs(1)^2*(GHat+iota*IHat)*VPrimeHat/(Delta^2*nHats(1)*sqrt(mHats(1))))
                    %}
                end
                fprintf('Total bootstrap current (FSABjHat) carried by species in this simulation: %g\n',FSABjHat)
                
                %{
                figure(99)
                clf
                contourf(zeta,theta,jHat,numContours,'EdgeColor','none')
                colorbar
                title('jHat')
                %}
                
                if RHSMode == 2
                    VPrimeHatWithG = VPrimeHat*(GHat+iota*IHat);
                    switch col
                        case 1
                            transportMatrix(1,1) = 4*(GHat+iota*IHat)*particleFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat^(3/2))*GHat);
                            transportMatrix(2,1) = 8*(GHat+iota*IHat)*heatFlux*nHat*B0OverBBar/(GHat*VPrimeHatWithG*(THat^(5/2))*GHat);
                            transportMatrix(3,1) = 2*nHat*FSABFlow/(GHat*THat);
                        case 2
                            transportMatrix(1,2) = 4*(GHat+iota*IHat)*particleFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*GHat);
                            transportMatrix(2,2) = 8*(GHat+iota*IHat)*heatFlux*B0OverBBar/(GHat*VPrimeHatWithG*sqrtTHat*THat*GHat);
                            transportMatrix(3,2) = 2*FSABFlow/(GHat);
                        case 3
                            transportMatrix(1,3) = particleFlux*Delta*Delta*FSABHat2/(VPrimeHatWithG*GHat*psiAHat*omega);
                            transportMatrix(2,3) = 2*Delta*Delta*heatFlux*FSABHat2/(GHat*VPrimeHatWithG*psiAHat*THat*omega);
                            transportMatrix(3,3) = FSABFlow*Delta*Delta*sqrtTHat*FSABHat2/((GHat+iota*IHat)*2*psiAHat*omega*B0OverBBar);
                    end
                end
            end
            
            if RHSMode == 2
                format longg
                transportMatrix
            end
            
            if testQuasisymmetryIsomorphism
                modifiedHeatFluxThatShouldBeConstant = heatFlux*abs(helicity_n/iota-helicity_l)/((helicity_n*IHat+helicity_l*GHat)^2);
                modifiedFSAFlowThatShouldBeConstant = FSABFlow*(helicity_n/iota-helicity_l)/(helicity_n*IHat+helicity_l*GHat);
                fprintf('   > Testing quasisymmetry isomorphism.\n')
                fprintf('   > Below are the two modified quantities that should be independent of helicity:\n')
                fprintf('   > Modified heat flux: %g\n',modifiedHeatFluxThatShouldBeConstant)
                fprintf('   > Modified FSA flow: %g\n',modifiedFSAFlowThatShouldBeConstant)
            end
            
            if programMode == 1
                figure(4+figureOffset)
                
                for ispecies = 1:Nspecies
                    
                    if smoothFigures
                        densityToUse = interp2(zeta,theta,densityPerturbation(:,:,ispecies),zetaFine,thetaFine,'cubic');
                        flowToUse = interp2(zeta,theta,flow(:,:,ispecies),zetaFine,thetaFine,'cubic');
                        pressureToUse = interp2(zeta,theta,pressurePerturbation(:,:,ispecies),zetaFine,thetaFine,'cubic');
                    else
                        densityToUse = densityPerturbation(:,:,ispecies);
                        flowToUse = flow(:,:,ispecies);
                        pressureToUse = pressurePerturbation(:,:,ispecies);
                    end
                    
                    subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
                    contourf(zetaToUse,thetaToUse,densityToUse,numContours,'EdgeColor','none')
                    colorbar
                    xlabel('\zeta')
                    ylabel('\theta')
                    title(['densityPerturbation, species ',num2str(ispecies)])
                    
                    subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
                    contourf(zetaToUse,thetaToUse,flowToUse,numContours,'EdgeColor','none')
                    colorbar
                    xlabel('\zeta')
                    ylabel('\theta')
                    title(['flow, species ',num2str(ispecies)])
                    
                    subplot(numRows,numCols,plotNum); plotNum=plotNum+1;
                    contourf(zetaToUse,thetaToUse,pressureToUse,numContours,'EdgeColor','none')
                    colorbar
                    xlabel('\zeta')
                    ylabel('\theta')
                    title(['pressurePerturbation, species ',num2str(ispecies)])
                end
            end
        end
        
        
        
        
        % --------------------------------------------------------
        % Below are some utilities for building sparse matrices.
        % --------------------------------------------------------
        
        function resetSparseCreator()
            sparseCreatorIndex=1;
            sparseCreator_i=zeros(estimated_nnz,1);
            sparseCreator_j=zeros(estimated_nnz,1);
            sparseCreator_s=zeros(estimated_nnz,1);
        end
        
        function addToSparse(i,j,s)
            n=numel(i);
            if n ~= numel(j)
                error('Error A');
            end
            if n ~= numel(s)
                error('Error B');
            end
            if any(i<1)
                error('Error Q: i<1');
            end
            if any(j<1)
                error('Error Q: j<1');
            end
            sparseCreator_i(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = i;
            sparseCreator_j(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = j;
            sparseCreator_s(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = s;
            sparseCreatorIndex = sparseCreatorIndex+n;
            if sparseCreatorIndex > estimated_nnz
                fprintf('Error! estimated_nnz is too small.\n')
            end
        end
        
        function addSparseBlock(rowIndices, colIndices, block)
            s=size(block);
            if (s(1) ~= numel(rowIndices)) || (s(2) ~= numel(colIndices))
                s
                size(rowIndices)
                size(colIndices)
                error('Error in addSparseBlock!')
            end
            [rows, cols, values] = find(block);
            addToSparse(rowIndices(rows),colIndices(cols),values)
        end
        
        function sparseMatrix = createSparse()
            fprintf('estimated nnz: %d   Actual value required: %d\n',estimated_nnz_original, sparseCreatorIndex)
            sparseMatrix = sparse(sparseCreator_i(1:(sparseCreatorIndex-1)), sparseCreator_j(1:(sparseCreatorIndex-1)), sparseCreator_s(1:(sparseCreatorIndex-1)), matrixSize, matrixSize);
            resetSparseCreator()
        end
        
        % *********************************************************
        % *********************************************************
        %
        % Below are utilities related to indexing the matrices:
        %
        % *********************************************************
        % *********************************************************
        
        
        % *********************************************************
        % For constraintScheme==0,
        % *********************************************************
        
        % Order of the rows of the matrix and of the RHS:
        % --------------------------------
        % for iSpecies = 1:Nspecies
        %   for ix = 1:Nx
        %     for L = 0:(NL-1)
        %       for itheta = 1:Ntheta
        %         for izeta = 1:Nzeta
        %           Enforce the drift-kinetic equation
        
        % Order of the vector of unknowns & of columns in the matrix:
        % --------------------------------
        % for iSpecies = 1:Nspecies
        %   for ix = 1:Nx
        %     for L = 0:(NL-1)
        %       for itheta = 1:Ntheta
        %         for izeta = 1:Nzeta
        %           Enforce the drift-kinetic equation
        
        % *********************************************************
        % For constraintScheme==1,
        % *********************************************************
        
        % Order of the rows of the matrix and of the RHS:
        % --------------------------------
        % for iSpecies = 1:Nspecies
        %   for ix = 1:Nx
        %     for L = 0:(NL-1)
        %       for itheta = 1:Ntheta
        %         for izeta = 1:Nzeta
        %           Enforce the drift-kinetic equation
        % for iSpecies = 1:Nspecies
        %   Force <n_1> = 0
        %   Force <p_1> = 0
        
        
        % Order of the vector of unknowns & of columns in the matrix:
        % --------------------------------
        % for iSpecies = 1:Nspecies
        %   for ix = 1:Nx
        %     for L = 0:(NL-1)
        %       for itheta = 1:Ntheta
        %         for izeta = 1:Nzeta
        %           Enforce the drift-kinetic equation
        % for iSpecies = 1:Nspecies
        %   particle source
        %   energy source
        
        % *********************************************************
        % For constraintScheme==2,
        % *********************************************************
        
        % Order of the rows of the matrix and of the RHS:
        % --------------------------------
        % for iSpecies = 1:Nspecies
        %   for ix = 1:Nx
        %     for L = 0:(NL-1)
        %       for itheta = 1:Ntheta
        %         for izeta = 1:Nzeta
        %           Enforce the drift-kinetic equation
        % for iSpecies = 1:Nspecies
        %   for ix = 1:Nx
        %     Force <f_1> = 0 at that x
        
        
        % Order of the vector of unknowns & of columns in the matrix:
        % --------------------------------
        % for iSpecies = 1:Nspecies
        %   for ix = 1:Nx
        %     for L = 0:(NL-1)
        %       for itheta = 1:Ntheta
        %         for izeta = 1:Nzeta
        %           Enforce the drift-kinetic equation
        % for iSpecies = 1:Nspecies
        %   for ix = 1:Nx
        %     source at that x
        
        function ii = getIndices(i_species, i_x, i_xi, i_theta, i_zeta, f_or_sources)
            % This is a very commonly used function which returns the
            % indices in the master matrix for given indices in the 4
            % coordinate grids (and species.) The basic idea for this
            % function is that of the input indices (i_species, i_x, i_xi,
            % i_theta, and i_zeta), each should either be a single value or
            % an array, but no more than one of these input indices can
            % be an array.
            
            % Allowed values for f_or_sources:
            % 0: distribution function or DKE
            % 1: particle source or <n_1> = 0
            % 2: heat source or <p_1> = 0
            % 3: <f_1> = 0 at each x
            
            moreThan1 = [numel(i_species)>1, numel(i_x)>1, numel(i_xi)>1, numel(i_theta)>1, numel(i_zeta)>1];
            if nnz(moreThan1)<0 || nnz(moreThan1)>1
                error('Either 0 or 1 of the index inputs to getIndices can have >1 elements.')
            end
            if any(i_species <=0)
                error('i_species must be positive')
            end
            if any(i_x <=0)
                error('i_x must be positive')
            end
            if any(i_xi <=0)
                error('i_xi must be positive')
            end
            if any(i_theta <=0)
                error('i_theta must be positive')
            end
            if any(i_zeta <=0)
                error('i_theta must be positive')
            end
            if any(i_species > Nspecies)
                error('i_species must be <= Nspecies')
            end
            if any(i_x > Nx)
                error('i_x must be <= Nx')
            end
            if any(i_xi > Nxi)
                error('i_xi must be <= Nxi')
            end
            if any(i_theta > Ntheta)
                error('i_theta must be <= Ntheta')
            end
            if any(i_zeta > Nzeta)
                error('i_zeta must be <= Nzeta')
            end
            switch f_or_sources
                case 0
                    % Distribution function or DKE:
                    ii = (i_species-1)*Nx*Nxi*Ntheta*Nzeta ...
                        + (i_x-1)*Nxi*Ntheta*Nzeta ...
                        + (i_xi-1)*Ntheta*Nzeta ...
                        + (i_theta-1)*Nzeta ...
                        + i_zeta;
                case 1
                    if constraintScheme ~= 1
                        error('f_or_sources should only be 1 if constraintScheme==1')
                    end
                    % particle source or <n_1> = 0:
                    ii = Nspecies*Nx*Nxi*Ntheta*Nzeta ...
                        + (i_species-1)*2 ...
                        + 1;
                case 2
                    if constraintScheme ~= 1
                        error('f_or_sources should only be 2 if constraintScheme==1')
                    end
                    % heat source or <p_1> = 0:
                    ii = Nspecies*Nx*Nxi*Ntheta*Nzeta ...
                        + (i_species-1)*2 ...
                        + 2;
                case 3
                    if constraintScheme ~= 2
                        error('f_or_sources should only be 3 if constraintScheme==2')
                    end
                    % <f_1> = 0 at each x:
                    ii = Nspecies*Nx*Nxi*Ntheta*Nzeta ...
                        + (i_species-1)*Nx ...
                        + i_x;
                otherwise
                    error('Invalid f_or_sources')
            end
            if any(ii > matrixSize)
                error('Something has gone wrong with indexing: index too big.')
            end
            if any(ii < 1)
                error('Something has gone wrong with indexing: index too small')
            end
        end
        
        % ------------------------------------------------------
        % ------------------------------------------------------
        % Below are routines to set the magnetic geometry.
        % ------------------------------------------------------
        % ------------------------------------------------------
        
        function setNPeriods()
            switch geometryScheme
                case 1
                    NPeriods = max([1, helicity_n]);
                case {2,3}
                    NPeriods = 10;
                case 4
                    NPeriods = 5;
                case 10
                    fid = fopen(fort996boozer_file);
                    if fid<0
                        error('Unable to open file %s\n',fort996boozer_file)
                    end
                    try
                      NPeriods = fscanf(fid,'%d',1);
                      fclose(fid);
                    catch me
                      error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec fort.996 output file.\n',...
                            me.message, fort996boozer_file)
                    end
                case 11
                    fid = fopen(JGboozer_file);
                    if fid<0
                        error('Unable to open file %s\n',JGboozer_file)
                    end
                    try
                        tmp_str=fgetl(fid);       %Skip comment line
                        while strcmp(tmp_str(1:2),'CC');
                            tmp_str=fgetl(fid);     %Skip comment line
                        end
                        header=fscanf(fid,'%d %d %d %d %f %f %f\n',7);
                        NPeriods = header(4);
                        fclose(fid);
                    catch me
                      error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec .bc output file.\n',...
                            me.message, JGboozer_file)
                    end
                case 12
                    fid = fopen(JGboozer_file_NonStelSym);
                    if fid<0
                        error('Unable to open file %s\n',JGboozer_file_NonStelSym)
                    end
                    try
                        tmp_str=fgetl(fid);       %Skip comment line
                        while strcmp(tmp_str(1:2),'CC');
                            tmp_str=fgetl(fid);     %Skip comment line
                        end
                        header=fscanf(fid,'%d %d %d %d %f %f %f\n',7);
                        NPeriods = header(4);
                        fclose(fid);
                    catch me
                      error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec .bc output file.\n',...
                            me.message, JGboozer_file_NonStelSym)
                    end
                otherwise
                    error('Invalid setting for geometryScheme')
            end
        end
        
        function computeBHat()
            % This subroutine does most of the work loading the magnetic geometry.
            
            [zeta2D, theta2D] = meshgrid(zeta,theta);
            
            switch geometryScheme
               case 1
                  % 2-helicity model:
                  BHarmonics_l = [1, helicity_l];
                  if helicity_n==0
                    BHarmonics_n = [0, 0];
                  else
                    BHarmonics_n = [0, 1];
                  end
                  BHarmonics_amplitudes = [epsilon_t, epsilon_h];
                  BHarmonics_parity = [1, 1];
                  
                  %{
                  % Note: the next line should probably be removed
                  % eventually:
                  omega = alpha*Delta*sqrt(mHats(1))/2;
                  dPhiHatdpsiN = EStar * iota * sqrt(THats) * psiAHat * B0OverBBar / (omega * GHat);
                  %}
               case 2
                  % LHD standard configuration.
                  % Values taken from Table 1 of
                  % Beidler et al, Nuclear Fusion 51, 076001 (2011).
                  iota = 0.4542;
                  BHarmonics_l = [1, 2, 1];
                  BHarmonics_n = [0, 1, 1];
                  BHarmonics_amplitudes = [-0.07053, 0.05067, -0.01476];
                  BHarmonics_parity = [1, 1, 1];
                  
                  B0OverBBar = 1; % (Tesla)
                  R0 = 3.7481; % (meters)
                  a = 0.5585; % (meters)
                  GHat = B0OverBBar * R0;
                  %IHat = GHat*3; % Change this to 0 eventually.
                  IHat = 0;
                  psiAHat = B0OverBBar*a^2/2;
                  dGdpHat=NaN;
                  
               case 3
                  % LHD inward-shifted configuration.
                  % Values taken from Table 1 of
                  % Beidler et al, Nuclear Fusion 51, 076001 (2011).
                  iota = 0.4692;
                  BHarmonics_l = [1, 2, 1, 0];
                  BHarmonics_n = [0, 1, 1, 1];
                  BHarmonics_amplitudes = [-0.05927, 0.05267, -0.04956, 0.01045];
                  BHarmonics_parity = [1, 1, 1, 1];
                  
                  B0OverBBar = 1; % (Tesla)
                  R0 = 3.6024; % (meters)
                  a = 0.5400; % (meters)
                  GHat = B0OverBBar * R0;
                  IHat = 0;
                  psiAHat = B0OverBBar*a^2/2;
                  dGdpHat=NaN;
                  
               case 4
                  % W7-X Standard configuration
                  % Values taken from Table 1 of
                  % Beidler et al, Nuclear Fusion 51, 076001 (2011).
                  iota=0.8700;
                  BHarmonics_l = [0, 1, 1];
                  BHarmonics_n = [1, 1, 0];
                  BHarmonics_amplitudes = [0.04645, -0.04351, -0.01902];
                  BHarmonics_parity = [1, 1, 1];
                  
                  B0OverBBar = 3.089; % (Tesla)
                  %R0 = 5.5267; % (meters)
                  %a = 0.5109; % (meters)
                  %psiAHat = pi*B0OverBBar*a^2;
                  GHat = -17.885;%B0OverBBar * R0;
                  IHat = 0;
                  psiAHat = -0.384935;
                  dGdpHat=NaN;
                  
               case 10
                  fid = fopen(fort996boozer_file);
                  % File description:
                  % 1st line: 2 integers:     nfp,ns
                  % 2nd line: 4 real numbers: aspect,rmax,rmin,betaxis
                  % 3rd line: 3 integers:     mboz, nboz, mnboz
                  % 4th line: 7 real numbers: iota,pres,beta,phip,phi,bvco,buco
                  %
                  % Then, you have 'mnboz' lines.
                  % If 'mn' is a dummy integer variable that goes from 1 to mnboz,
                  % for each value of mn you read
                  %
                  % m(mn),n(mn),bmn(mn),rmnc(mn),zmns(mn)pmns(m,n),gmn(mn)
                  try
                    header=fscanf(fid,'%d %d\n %f %f %f %f\n %d %d %d %f %f %f %f %f %f %f',16);
                    mnboz=header(9);
                    modes =fscanf(fid,'%d %d %g %g %g %g %g',[7,mnboz]);
                    fclose(fid);
                    
                    % scalar values
                    %Nper = header(1); %number of field periods
                    iota = header(10);
                    Ihat = header(16);  % Covariant theta comp. of B, known as I in sfincs (in meter * Tesla)
                    Ghat = header(15);  % Covariant phi comp. of B, known as G in sfincs (in meter * Tesla)
                    % Note that the flux at the separatrix is not stored in the
                    % file, so we set PsiAHat in the Physics parameters
                    % section in the beginning of the program
                    
                    % mode amplitudes
                    if modes(1,1)==0 && modes(2,1)==0 
                      B0OverBBar=modes(3,1); %The B00 component in Tesla
                    else
                      error('The first fort996boozer_file entry is not the B00 component')
                    end
                    BHarmonics_l = modes(1,2:end);
                    BHarmonics_n = modes(2,2:end);
                    BHarmonics_amplitudes = modes(3,2:end)/B0OverBBar; % Store the values normalised to the B00 component. 
                    BHarmonics_parity = ones(1,length(BHarmonics_amplitudes));
                    dGdpHat=NaN; %Not implemented yet                    
                  catch me
                    error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec fort.996 output file.\n',...
                        me.message, fort996boozer_file)
                  end
               case 11
              
                  fid = fopen(JGboozer_file);
                  if fid<0
                      error('Unable to open file %s\n',JGboozer_file)
                  end
                  
                  try
                      tmp_str=fgetl(fid);
                      while strcmp(tmp_str(1:2),'CC');
                          tmp_str=fgetl(fid); %Skip comment line
                      end
                      header=fscanf(fid,'%d %d %d %d %f %f %f\n',7);
                      fgetl(fid);  %Skip variable name line
                      
                      NPeriods = header(4);
                      psiAHat  = header(5)/2/pi; %Convert the flux from Tm^2 to Tm^2/rad

                      max_no_of_modes=500;
                      modesm_new=NaN*zeros(1,max_no_of_modes);
                      modesn_new=NaN*zeros(1,max_no_of_modes);
                      modesb_new=NaN*zeros(1,max_no_of_modes);
                      normradius_new=-inf;
                      no_of_modes_new=NaN;
                      iota_new=NaN;
                      G_new=NaN;
                      I_new=NaN;
                      pPrimeHat_new=NaN;
                      end_of_file=0;

                      while (normradius_new<normradius_wish) && not(end_of_file)
                          normradius_old=normradius_new;
                          no_of_modes_old=no_of_modes_new;
                          modesm_old=modesm_new;
                          modesn_old=modesn_new;
                          modesb_old=modesb_new;
                          iota_old=iota_new;
                          G_old=G_new;
                          I_old=I_new;
                          pPrimeHat_old=pPrimeHat_new;
                          
                          fgetl(fid);
                          surfheader=fscanf(fid,'%f %f %f %f %f %f\n',6);
                          
                          normradius_new=sqrt(surfheader(1));
                          iota_new=surfheader(2);
                          G_new=surfheader(3)*NPeriods/2/pi*(4*pi*1e-7); %Tesla*meter
                          I_new=surfheader(4)/2/pi*(4*pi*1e-7);          %Tesla*meter
                          pPrimeHat_new=surfheader(5)*(4*pi*1e-7);       % p=pHat \bar{B}^2 / \mu_0
                          
                          fgetl(fid); %Skip units line
                          proceed=1;
                          modeind=0;
                          while proceed
                              tmp_str=fgetl(fid);
                              if length(tmp_str)==1
                                  if tmp_str==-1 %End of file has been reached
                                      proceed=0;
                                      end_of_file=1;
                                  end
                              elseif not(isempty(find(tmp_str=='s'))) %Next flux surface has been reached
                                  proceed=0;
                              else
                                  tmp=sscanf(tmp_str,'%d %d %f %f %f %f',6);
                                  if abs(tmp(6))>min_Bmn_to_load
                                      modeind=modeind+1;
                                      %if modeind > max_no_of_modes %Unnecessary to check this in matlab
                                      %  error(' modeind > max_no_of_modes !')
                                      %end
                                      modesm_new(modeind)=tmp(1);
                                      modesn_new(modeind)=tmp(2);
                                      modesb_new(modeind)=tmp(6);
                                  end
                              end
                          end
                          no_of_modes_new=modeind;
                          modesm_new(no_of_modes_new+1:end)=NaN;
                          modesn_new(no_of_modes_new+1:end)=NaN;
                          modesb_new(no_of_modes_new+1:end)=NaN;
                      end
                      fclose(fid);
                  catch me
                      error('%s\n\nFile\n\t%s\ndoes not seem to be a valid .bc geometry file.\n',...
                          me.message, JGboozer_file)
                  end

                  [~,minind]=min([(normradius_old-normradius_wish)^2,...
                                  (normradius_new-normradius_wish)^2]);
                  if minind==1
                    BHarmonics_l = modesm_old(1:no_of_modes_old);
                    BHarmonics_n = modesn_old(1:no_of_modes_old);
                    BHarmonics_amplitudes = modesb_old(1:no_of_modes_old);
                    iota=iota_old;
                    GHat=G_old;
                    IHat=I_old;
                    pPrimeHat=pPrimeHat_old;
                    normradius=normradius_old;
                  else %minind=2
                    BHarmonics_l = modesm_new(1:no_of_modes_new);
                    BHarmonics_n = modesn_new(1:no_of_modes_new);
                    BHarmonics_amplitudes = modesb_new(1:no_of_modes_new);
                    iota=iota_new;
                    GHat=G_new;
                    IHat=I_new;
                    pPrimeHat=pPrimeHat_new;
                    normradius=normradius_new;
                  end
                  dGdpHat=(G_new-G_old)/(normradius_new^2-normradius_old^2)/pPrimeHat;
                  
                  disp(['The calculation is performed for the normalised radius ',num2str(normradius)])

                  m0inds=find(BHarmonics_l==0);
                  n0m0inds=find(BHarmonics_n(m0inds)==0);
                  if isempty(n0m0inds)
                    error(' B00 component is missing!')
                  end
                  nm00ind=m0inds(n0m0inds);
                  B0OverBBar=BHarmonics_amplitudes(nm00ind); % Assumes \bar{B} = 1 Tesla
                  BHarmonics_amplitudes=[BHarmonics_amplitudes(1:nm00ind-1),...
                                         BHarmonics_amplitudes(nm00ind+1:end)]...
                                        /B0OverBBar;
                  BHarmonics_l = [BHarmonics_l(1:nm00ind-1),...
                                  BHarmonics_l(nm00ind+1:end)];
                  BHarmonics_n = [BHarmonics_n(1:nm00ind-1),...
                                  BHarmonics_n(nm00ind+1:end)];
                  BHarmonics_parity = ones(1,length(BHarmonics_amplitudes));
             case 12
                  %Non-stellarator symmetric case
                  fid = fopen(JGboozer_file_NonStelSym);
                  if fid<0
                      error('Unable to open file %s\n',JGboozer_file_NonStelSym)
                  end
                  
                  try
                      tmp_str=fgetl(fid);
                      while strcmp(tmp_str(1:2),'CC');
                          tmp_str=fgetl(fid); %Skip comment line
                      end
                      header=fscanf(fid,'%d %d %d %d %f %f %f\n',7);
                      fgetl(fid);  %Skip variable name line
                      
                      NPeriods = header(4);
                      psiAHat  = header(5)/2/pi; %Convert the flux from Tm^2 to Tm^2/rad
                      a        = header(6);      %minor radius %m
                      
                      max_no_of_modes=1000;
                      modesm_new=NaN*zeros(1,max_no_of_modes);
                      modesn_new=NaN*zeros(1,max_no_of_modes);
                      modesb_new=NaN*zeros(1,max_no_of_modes);
                      normradius_new=-inf;
                      no_of_modes_new=NaN;
                      iota_new=NaN;
                      G_new=NaN;
                      I_new=NaN;
                      pPrimeHat_new=NaN;
                      end_of_file=0;
                      
                      while (normradius_new<normradius_wish) && not(end_of_file)
                          normradius_old=normradius_new;
                          no_of_modes_old=no_of_modes_new;
                          modesm_old=modesm_new;
                          modesn_old=modesn_new;
                          modesb_old=modesb_new;
                          iota_old=iota_new;
                          G_old=G_new;
                          I_old=I_new;
                          pPrimeHat_old=pPrimeHat_new;
                          
                          fgetl(fid);
                          surfheader=fscanf(fid,'%f %f %f %f %f %f\n',6);
                          
                          normradius_new=sqrt(surfheader(1)); %r/a=sqrt(psi/psi_a)
                          iota_new=surfheader(2);
                          G_new=surfheader(3)*NPeriods/2/pi*(4*pi*1e-7); %Tesla*meter
                          I_new=surfheader(4)/2/pi*(4*pi*1e-7);          %Tesla*meter
                          pPrimeHat_new=surfheader(5)*(4*pi*1e-7);       % p=pHat \bar{B}^2 / \mu_0
                          
                          fgetl(fid); %Skip units line
                          proceed=1;
                          modeind=0;
                          while proceed
                              tmp_str=fgetl(fid);
                              if length(tmp_str)==1
                                  if tmp_str==-1 %End of file has been reached
                                      proceed=0;
                                      end_of_file=1;
                                  end
                              elseif not(isempty(find(tmp_str=='s'))) %Next flux surface has been reached
                                  proceed=0;
                              else
                                  tmp=sscanf(tmp_str,'%d %d %f %f %f %f %f %f %f %f',10);
                                  if (abs(tmp(9))>min_Bmn_to_load) || (abs(tmp(10))>min_Bmn_to_load)
                                      modeind=modeind+1;
                                      modesm_new(modeind)=tmp(1);
                                      modesn_new(modeind)=tmp(2);
                                      modesb_new(modeind)=tmp(9); %Cosinus component
                                      
                                      modeind=modeind+1;
                                      modesm_new(modeind)=tmp(1);
                                      modesn_new(modeind)=tmp(2);
                                      modesb_new(modeind)=tmp(10); %Sinus component
                                  end
                              end
                          end
                          no_of_modes_new=modeind;
                          modesm_new(no_of_modes_new+1:end)=NaN;
                          modesn_new(no_of_modes_new+1:end)=NaN;
                          modesb_new(no_of_modes_new+1:end)=NaN;
                      end
                      fclose(fid);
                  catch me
                      error('%s\n\nFile\n\t%s\ndoes not seem to be a valid .bc geometry file.\n',...
                          me.message, JGboozer_file_NonStelSym)
                  end

                  [~,minind]=min([(normradius_old-normradius_wish)^2,...
                                  (normradius_new-normradius_wish)^2]);
                  if minind==1
                    BHarmonics_l = modesm_old(1:no_of_modes_old);
                    BHarmonics_n = modesn_old(1:no_of_modes_old);
                    BHarmonics_amplitudes = modesb_old(1:no_of_modes_old);
                    iota=iota_old;
                    GHat=G_old;
                    IHat=I_old;
                    pPrimeHat=pPrimeHat_old;
                    normradius=normradius_old;
                  else %minind=2
                    BHarmonics_l = modesm_new(1:no_of_modes_new);
                    BHarmonics_n = modesn_new(1:no_of_modes_new);
                    BHarmonics_amplitudes = modesb_new(1:no_of_modes_new);
                    iota=iota_new;
                    GHat=G_new;
                    IHat=I_new;
                    pPrimeHat=pPrimeHat_new;
                    normradius=normradius_new;
                  end
                  dGdpHat=(G_new-G_old)/(normradius_new^2-normradius_old^2)/pPrimeHat;
                  
                  disp(['The calculation is performed for radius ' ...
                        ,num2str(normradius*a),' m , r/a=',num2str(normradius)])
                  
                  m0inds=find(BHarmonics_l==0);
                  n0m0inds=find(BHarmonics_n(m0inds)==0);
                  if isempty(n0m0inds)
                    error(' B00 component is missing!')
                  end
                  nm00ind=m0inds(n0m0inds(1));
                  B0OverBBar=BHarmonics_amplitudes(nm00ind); %Assumes \bar{B}=1T
                  BHarmonics_amplitudes=[BHarmonics_amplitudes(1:nm00ind-1),...
                                         BHarmonics_amplitudes(nm00ind+2:end)]...
                                        /B0OverBBar;
                  BHarmonics_l = [BHarmonics_l(1:nm00ind-1),...
                                  BHarmonics_l(nm00ind+2:end)];
                  BHarmonics_n = [BHarmonics_n(1:nm00ind-1),...
                                  BHarmonics_n(nm00ind+2:end)];
                  BHarmonics_parity=((-1).^(0:length(BHarmonics_n)-1)+1)/2; %[1,0,1,0,1,0,1,0,...], i.e. cos,sin.cos,sin,...
                  
               otherwise
                  error('Invalid setting for geometryScheme')
            end
                
            NHarmonics = numel(BHarmonics_amplitudes);
            BHat = B0OverBBar * ones(Ntheta,Nzeta);
            dBHatdtheta = zeros(Ntheta,Nzeta);
            dBHatdzeta = zeros(Ntheta,Nzeta);
            for i=1:NHarmonics
              if BHarmonics_parity(i) %The cosine components of BHat
                BHat = BHat + B0OverBBar * BHarmonics_amplitudes(i) *...
                       cos(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
                dBHatdtheta = dBHatdtheta - B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_l(i) *...
                    sin(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
                dBHatdzeta = dBHatdzeta + B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_n(i) * NPeriods *...
                    sin(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods ...
                        * zeta2D);
              else  %The sine components of BHat
                BHat = BHat + B0OverBBar * BHarmonics_amplitudes(i) *...
                       sin(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
                dBHatdtheta = dBHatdtheta + B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_l(i) *...
                    cos(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
                dBHatdzeta = dBHatdzeta - B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_n(i) * NPeriods *...
                    cos(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods ...
                        * zeta2D);                  
              end       
            end
            % ---------------------------------------------------------------------------------------
            % Calculate parallel current u from harmonics of 1/B^2. Used in NTV calculation.
            % \nabla_\parallel u = (2/B^4) \nabla B \times \vector{B} \cdot \iota \nabla \psi 
            % ---------------------------------------------------------------------------------------
            uHat = zeros(Ntheta,Nzeta);
            duHatdtheta = zeros(Ntheta,Nzeta);
            duHatdzeta = zeros(Ntheta,Nzeta);
            hHat=1./(BHat.^2);
            if any(BHarmonics_parity==0) %sine components exist
              for m=0:floor(Ntheta/2)-1 %Nyquist max freq.
                for n=0:floor(Nzeta/2)-1
                  if not(m==0 && n==0)
                    %cos
                    hHatHarmonics_amplitude = 2/(Ntheta*Nzeta) *...
                      sum(sum(cos(m * theta2D  - n * NPeriods * zeta2D).*hHat));
                    uHatHarmonics_amplitude = ...
                        iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude;
                    uHat = uHat + uHatHarmonics_amplitude * cos(m * theta2D - n * NPeriods * zeta2D);
                    duHatdtheta = duHatdtheta ...
                        - uHatHarmonics_amplitude * m * sin(m * theta2D - n * NPeriods * zeta2D);
                    duHatdzeta = duHatdzeta ...
                      + uHatHarmonics_amplitude * n * NPeriods * sin(m * theta2D - n * NPeriods * zeta2D); 
                    
                    %sin
                    hHatHarmonics_amplitude = 2/(Ntheta*Nzeta) *...
                        sum(sum(sin(m * theta2D  - n * NPeriods * zeta2D).*hHat));
                    uHatHarmonics_amplitude = ...
                        iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude;
                    uHat = uHat + uHatHarmonics_amplitude * sin(m * theta2D - n * NPeriods * zeta2D);
                    duHatdtheta = duHatdtheta ...
                        + uHatHarmonics_amplitude * m * cos(m * theta2D - n * NPeriods * zeta2D);
                    duHatdzeta = duHatdzeta ...
                        - uHatHarmonics_amplitude * n * NPeriods * cos(m * theta2D - n * NPeriods * zeta2D);   
                  end
                end
              end
            else %only cosinus components
              for m=0:Ntheta-1 %cos-series max freq.
                for n=0:Nzeta-1
                  if not(m==0 && n==0)
                    hHatHarmonics_amplitude = 2/(Ntheta*Nzeta) *...
                      sum(sum(cos(m * theta2D  - n * NPeriods * zeta2D).*hHat));
                    uHatHarmonics_amplitude = ...
                        iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude;
                    uHat = uHat + uHatHarmonics_amplitude * cos(m * theta2D - n * NPeriods * zeta2D);
                    duHatdtheta = duHatdtheta ...
                        - uHatHarmonics_amplitude * m * sin(m * theta2D - n * NPeriods * zeta2D);
                    duHatdzeta = duHatdzeta ...
                      + uHatHarmonics_amplitude * n * NPeriods * sin(m * theta2D - n * NPeriods * zeta2D);   
                  end
                end              
              end
            end
            NTVkernel = 2/5 * ( ...
                dGdpHat ./ BHat .* (iota * dBHatdtheta + dBHatdzeta) + ...
                1/2 * (iota * (duHatdtheta + uHat * 2./BHat .* dBHatdtheta) ...
                          + duHatdzeta + uHat * 2./BHat .* dBHatdzeta) );
        end  
    end
end
        



