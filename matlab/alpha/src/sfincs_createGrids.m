function sfincs_createGrids()

global Nalpha Nzeta Nx Nxi NL NxPotentialsPerVth xMax solverTolerance xGrid_k
global constraintScheme collisionOperator NPeriods pointAtX0
global alphaDerivativeScheme zetaDerivativeScheme forceOddNalphaAndNzeta
global alpha ddalpha alphaWeights ddalpha_preconditioner preconditioner_alpha
global zeta ddzeta zetaWeights ddzeta_preconditioner preconditioner_zeta
global x xWeights ddx d2dx2 ddx_preconditioner xGridScheme xPotentialsGridScheme
global alpha2D zeta2D RHSMode preconditioner_x xMaxPotentials
global xInterpolationScheme xPotentialsInterpolationScheme NxPotentials
global xPotentials ddxPotentials d2dx2Potentials interpolateXToXPotentials
global indexVars Nspecies includePhi1 transportMatrix
global BLOCK_F BLOCK_QN BLOCK_PHI1_CONSTRAINT BLOCK_DENSITY_CONSTRAINT BLOCK_PRESSURE_CONSTRAINT BLOCK_F_CONSTRAINT
global zeta_to_impose_DKE zetaMax

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
    
if forceOddNalphaAndNzeta
    if mod(Nalpha,2)==0
        Nalpha=Nalpha+1;
    end
    if mod(Nzeta,2)==0
        Nzeta=Nzeta+1;
    end
end

indexVars = struct();
indexVars.Nalpha = Nalpha;
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
fprintf('            Nalpha = %d\n',Nalpha)
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
% Generate abscissae, quadrature weights, and derivative matrix for alpha grid.
% *************************************************************************

switch alphaDerivativeScheme
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
        error('Error! Invalid alphaDerivativeScheme')
end
[alpha, alphaWeights, ddalpha, ~] = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, scheme);

switch preconditioner_alpha
    case 0
        ddalpha_preconditioner = ddalpha;
    case 1
        % Uniform periodic finite differences with 3-point stencil
        scheme = 0;
        [~, ~, ddalpha_preconditioner, ~] = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, scheme);
    case 2
        ddalpha_preconditioner = zeros(size(ddalpha));
    otherwise
        error('Invalid preconditioner_alpha')
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
        case 1
            zeta_stencil_size=3;
            zeta_to_impose_DKE = 2:(Nzeta-1);
            Delta_zeta = (2*pi)/(NPeriods*(Nzeta-2));
            zeta_scheme=2;
            [zeta, ~, ddzeta, ~] = sfincs_uniformDiffMatrices(Nzeta, -Delta_zeta, zetaMax, zeta_scheme);
            assert(abs(zeta(2)-zeta(1)-Delta_zeta)<1e-12)
            zetaWeights=ones(size(zeta));
            zetaWeights(1)  =0;
            zetaWeights(end)=0;
            zetaWeights = zetaWeights * Delta_zeta * NPeriods;
        case 2
            zeta_stencil_size=5;
            zeta_to_impose_DKE = 3:(Nzeta-2);
            Delta_zeta = (2*pi)/(NPeriods*(Nzeta-4));
            zeta_scheme=12;
            [zeta, ~, ddzeta, ~] = sfincs_uniformDiffMatrices(Nzeta, -2*Delta_zeta, zetaMax+Delta_zeta, zeta_scheme);
            assert(abs(zeta(2)-zeta(1)-Delta_zeta)<1e-12)
            zetaWeights=ones(size(zeta));
            zetaWeights(1:2)      =0;
            zetaWeights(end-1:end)=0;
            zetaWeights = zetaWeights * Delta_zeta * NPeriods;
        otherwise
            error('Invalid zetaDerivativeScheme')
    end

    
    
   
    
    switch preconditioner_zeta
        case 0
            ddzeta_preconditioner = ddzeta;
        case 1
            % Uniform periodic finite differences with 3-point stencil
            if zetaDerivativeScheme==1
                ddzeta_preconditioner = ddzeta;
            else
                zeta_scheme = 2;
                [~, ~, ddzeta_preconditioner, ~] = sfincs_uniformDiffMatrices(Nzeta, -2*Delta_zeta, zetaMax+Delta_zeta, zeta_scheme);
            end
        case 2
            ddzeta_preconditioner = zeros(size(ddzeta));
        otherwise
            error('Invalid preconditioner_zeta')
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
            if pointAtX0
                weightFactors(:,1) = [0, -2];
            end
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

[zeta2D, alpha2D] = meshgrid(zeta,alpha);

sfincs_computeBHat()
sfincs_computeBIntegrals()

fprintf('Done creating grids.\n')

end
