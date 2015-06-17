function sfincs_createGrids()

global Ntheta Nzeta Nx Nxi NL NxPotentialsPerVth xMax solverTolerance xGrid_k
global constraintScheme collisionOperator NPeriods pointAtX0
global thetaDerivativeScheme zetaDerivativeScheme forceOddNthetaAndNzeta
global theta ddtheta thetaWeights ddtheta_preconditioner preconditioner_theta
global zeta ddzeta zetaWeights ddzeta_preconditioner preconditioner_zeta
global x xWeights ddx d2dx2 ddx_preconditioner xGridScheme xPotentialsGridScheme
global theta2D zeta2D RHSMode preconditioner_x xMaxPotentials
global xInterpolationScheme xPotentialsInterpolationScheme NxPotentials
global xPotentials ddxPotentials d2dx2Potentials interpolateXToXPotentials
global indexVars Nspecies includePhi1 transportMatrix
global BLOCK_F BLOCK_QN BLOCK_PHI1_CONSTRAINT BLOCK_DENSITY_CONSTRAINT BLOCK_PRESSURE_CONSTRAINT BLOCK_F_CONSTRAINT

% *************************************************************************
% Do a few sundry initialization tasks:
% *************************************************************************

if constraintScheme < 0
    if collisionOperator == 0
        constraintScheme = 1;
    else
        constraintScheme = 2;
    end
end 
    
if forceOddNthetaAndNzeta
    if mod(Ntheta,2)==0
        Ntheta=Ntheta+1;
    end
    if mod(Nzeta,2)==0
        Ntheta=Ntheta+1;
    end
end

indexVars = struct();
indexVars.Ntheta = Ntheta;
indexVars.Nzeta = Nzeta;
indexVars.Nspecies = Nspecies;
indexVars.Nx = Nx;
indexVars.Nxi = Nxi;
indexVars.includePhi1 = includePhi1;
indexVars.BLOCK_F = BLOCK_F;
indexVars.BLOCK_QN = BLOCK_QN;
indexVars.BLOCK_PHI1_CONSTRAINT = BLOCK_PHI1_CONSTRAINT;
indexVars.BLOCK_DENSITY_CONSTRAINT = BLOCK_DENSITY_CONSTRAINT;
indexVars.BLOCK_PRESSURE_CONSTRAINT = BLOCK_PRESSURE_CONSTRAINT;
indexVars.BLOCK_F_CONSTRAINT = BLOCK_F_CONSTRAINT;

if RHSMode==2
    transportMatrix = zeros(3);
elseif RHSMode==3
    transportMatrix = zeros(2);
end

fprintf('---- Numerical parameters: ----\n')
fprintf('            Ntheta = %d\n',Ntheta)
fprintf('             Nzeta = %d\n',Nzeta)
fprintf('               Nxi = %d\n',Nxi)
fprintf('                NL = %d\n',NL)
fprintf('                Nx = %d\n',Nx)
if xGridScheme<5
    fprintf('NxPotentialsPerVth = %d\n',NxPotentialsPerVth)
    fprintf('              xMax = %d\n',xMax)
end
fprintf('   solverTolerance = %d\n',solverTolerance)

sfincs_computeMatrixSize()

% *************************************************************************
% Generate abscissae, quadrature weights, and derivative matrix for theta grid.
% *************************************************************************

switch thetaDerivativeScheme
    case 0
        % Spectral uniform
        scheme = 20;
    case 1
        % Uniform periodic finite differences with 3-point stencil
        scheme = 0;
    case 2
        % Uniform periodic finite differences with 5-point stencil
        scheme = 10;
    otherwise
        error('Error! Invalid thetaDerivativeScheme')
end
[theta, thetaWeights, ddtheta, ~] = sfincs_uniformDiffMatrices(Ntheta, 0, 2*pi, scheme);

switch preconditioner_theta
    case 0
        ddtheta_preconditioner = ddtheta;
    case 1
        % Uniform periodic finite differences with 3-point stencil
        scheme = 0;
        [~, ~, ddtheta_preconditioner, ~] = sfincs_uniformDiffMatrices(Ntheta, 0, 2*pi, scheme);
    case 2
        ddtheta_preconditioner = zeros(size(ddtheta));
    case 3
        ddtheta_preconditioner = eye(Ntheta);
end
    


% *************************************************************************
% Generate abscissae, quadrature weights, and derivative matrix for zeta grid.
% *************************************************************************
zetaMax = 2*pi/NPeriods;

if Nzeta==1
    zeta=0;
    zetaWeights=2*pi;
    ddzeta=0;
    ddzeta_preconditioner=0;
else
    switch zetaDerivativeScheme
        case 0
            % Spectral uniform
            scheme = 20;
        case 1
            % Uniform periodic finite differences with 3-point stencil
            scheme = 0;
        case 2
            % Uniform periodic finite differences with 5-point stencil
            scheme = 10;
        otherwise
            error('Error! Invalid thetaDerivativeScheme')
    end
    [zeta, zetaWeights, ddzeta, ~] = sfincs_uniformDiffMatrices(Nzeta, 0, zetaMax, scheme);
    zetaWeights = zetaWeights * NPeriods;
    
    switch preconditioner_zeta
        case 0
            ddzeta_preconditioner = ddzeta;
        case 1
            % Uniform periodic finite differences with 3-point stencil
            scheme = 0;
            [~, ~, ddzeta_preconditioner, ~] = sfincs_uniformDiffMatrices(Nzeta, 0, zetaMax, scheme);
        case 2
            ddzeta_preconditioner = zeros(size(ddzeta));
        case 3
            ddzeta_preconditioner = eye(Nzeta);
    end
end

% *************************************************************************
% Generate abscissae, quadrature weights, and derivative matrices for
% the energy (x) grid used to represent the distribution function.
% *************************************************************************

switch xGridScheme
    case {1,2,5,6}
        % For these values of xGridScheme, xInterpolationScheme does not
        % matter.
        xInterpolationScheme = -1;
    case 3
        xInterpolationScheme = 1;
    case 4
        xInterpolationScheme = 2;
    otherwise
        error('Invalid xGridScheme')
end

switch xPotentialsGridScheme
    case {1,3}
        xPotentialsInterpolationScheme = 1;
    case {2,4}
        xPotentialsInterpolationScheme = 2;
    otherwise
        error('Invalid xPotentialsInterpolationScheme')
end

% Set pointAtX0, x, xWeights, ddx, d2dx2
if RHSMode==3
    % Monoenergetic calculation
    pointAtX0 = false;
    x = 1;
    xWeights = exp(1);
    ddx = 0;
    d2dx2 = 0;
    
else
    % RHSMode = 1 or 2.
    switch xGridScheme
        case {1,2,5,6}
            % Pseudospectral x grid
            
            if xGridScheme==1 || xGridScheme==5
                pointAtX0 = false;
            else
                if xGrid_k ~= 0
                    error('If xGridScheme is 2 or 6, then xGrid_k must be 0.')
                end
                pointAtX0 = true;
            end
            [x, xWeights, polynomials_a, polynomials_b, polynomials_c] = sfincs_GaussWeightsAndAbscissae(Nx, @sfincs_xWeight, 0, Inf, pointAtX0);
            xWeights = xWeights./sfincs_xWeight(x);
            weightFactors=[xGrid_k./x-2*x, xGrid_k*(xGrid_k-1)./(x.*x)-2*(2*xGrid_k+1)+4*x.*x]';
            differentiationMatrices = sfincs_poldif(x,sfincs_xWeight(x),weightFactors);
            ddx=differentiationMatrices(:,:,1);
            d2dx2=differentiationMatrices(:,:,2);

            if xPotentialsGridScheme == 3 || xPotentialsGridScheme == 4
                error('When xPotentialsGridScheme == 3 or 4, you must set xGridScheme = 3 or 4.')
            end
            
        case {3,4}
            % Uniform grid on [0, xMax]
            
            pointAtX0 = true;
            scheme = 12;
            [x_plus1, xWeights_plus1, ddx_plus1, d2dx2_plus1] = sfincs_uniformDiffMatrices(Nx+1, 0, xMax, scheme);
            % Discard the last point (the one at xMax)
            x = x_plus1(1:Nx);
            xWeights = xWeights_plus1(1:Nx);
            ddx = ddx_plus1(1:Nx, 1:Nx);
            d2dx2 = d2dx2_plus1(1:Nx, 1:Nx);
            
        otherwise
            error('Invalid xGridScheme')
    end
end

xMaxPotentials = max([xMax, max(x)]);
if xPotentialsGridScheme == 3 || xPotentialsGridScheme == 4
    NxPotentials = Nx+1;
else
    NxPotentials = ceil(NxPotentialsPerVth * xMaxPotentials);
end

%{
x
xWeights
ddx
%}

% Make the energy grid and differentiation matrices for the
% Rosenbluth potentials, and the interpolation matrix from x to xPotentials
if RHSMode == 3
    % Monoenergetic
    xPotentials = 0;
    ddxPotentials = 0;
    d2dx2Potentials = 0;
    interpolateXToXPotentials = 0;
    
else
    % Normal (not monoenergetic)
    scheme = 12;
    [xPotentials, ~, ddxPotentials, d2dx2Potentials] = sfincs_uniformDiffMatrices(NxPotentials, 0, xMaxPotentials, scheme);
    switch xGridScheme
        case {1,2,5,6}
            % Pseudospectral x grid
            interpolateXToXPotentials = sfincs_polynomialInterpolationMatrix(x,xPotentials,sfincs_xWeight(x),sfincs_xWeight(xPotentials));
        case {3,4}
            error('Not implemented yet')
        otherwise
            error('Invalid xGridScheme')
    end
end

switch preconditioner_x
    case 0
        % No simplification
        ddx_preconditioner = ddx;
    case 1
        % Keep only diagonal
        ddx_preconditioner = diag(diag(ddx));
    case 2
        ddx_preconditioner = triu(ddx);
    case 3
        mask = eye(Nx) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);
        ddx_preconditioner = ddx .* mask;
    case 4
        mask = eye(Nx) + diag(ones(Nx-1,1),1);
        ddx_preconditioner = ddx .* mask;
    otherwise
        error('Invalid preconditioner_x')
end

if RHSMode ~= 3 && (xGridScheme==5 || xGridScheme==6)
    sfincs_RosenbluthPotentialResponse(polynomials_a, polynomials_b, polynomials_c)
end

% *************************************************************************
% Compute the magnetic field on the grid.
% Also compute a few quantities related to the magnetic field
% *************************************************************************

[zeta2D, theta2D] = meshgrid(zeta,theta);

sfincs_computeBHat()
sfincs_computeBIntegrals()

fprintf('Done creating grids.\n')

end