function sfincs_createGrids()

global Nalpha Nzeta Nx Nxi NL NxPotentialsPerVth xMax solverTolerance xGrid_k
global constraintScheme collisionOperator NPeriods pointAtX0
global forceOddNalphaAndNzeta
global streaming_theta_derivative_option preconditioner_streaming_theta_derivative_option
global alpha streaming_ddtheta_sum streaming_ddtheta_difference streaming_ddtheta_sum_preconditioner streaming_ddtheta_difference_preconditioner alphaWeights
global streaming_zeta_derivative_option preconditioner_streaming_zeta_derivative_option
global zeta streaming_ddzeta_sum streaming_ddzeta_difference streaming_ddzeta_sum_preconditioner streaming_ddzeta_difference_preconditioner zetaWeights
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

switch streaming_theta_derivative_option
    case 2
        fprintf('Streaming d/dtheta discretized using centered differences: 1 point on each side.\n')
        derivative_option_plus = 2;
        derivative_option_minus = derivative_option_plus;
    case 3
        fprintf('Streaming d/dtheta discretized using centered differences: 2 points on each side.\n')
        derivative_option_plus = 12;
        derivative_option_minus = derivative_option_plus;
    case 4
        fprintf('Streaming d/dtheta discretized using upwinded differences: 0 points on one side, 1 point on the other.\n')
        derivative_option_plus  = 32;
        derivative_option_minus = 42;
    case 5
        fprintf('Streaming d/dtheta discretized using upwinded differences: 0 points on one side, 2 points on the other.\n')
        derivative_option_plus  = 52;
        derivative_option_minus = 62;
    case 6
        fprintf('Streaming d/dtheta discretized using upwinded differences: 1 point on one side, 2 points on the other.\n')
        derivative_option_plus  = 82;
        derivative_option_minus = 92;
    case 7
        fprintf('Streaming d/dtheta discretized using upwinded differences: 1 point on one side, 3 points on the other.\n')
        derivative_option_plus  = 102;
        derivative_option_minus = 112;
    case 8
        fprintf('Streaming d/dtheta discretized using upwinded differences: 2 points on one side, 3 points on the other.\n')
        derivative_option_plus  = 122;
        derivative_option_minus = 132;
    otherwise
        error('Invalid streaming_theta_derivative_option: %d',streaming_theta_derivative_option)
end
quadrature_option = 0;
[alpha, alphaWeights, streaming_ddtheta_plus, ~] = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, derivative_option_plus,  quadrature_option);
[~, ~, streaming_ddtheta_minus, ~] = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, derivative_option_minus, quadrature_option);

streaming_ddtheta_sum        = 0.5*(streaming_ddtheta_plus+streaming_ddtheta_minus);
streaming_ddtheta_difference = 0.5*(streaming_ddtheta_plus-streaming_ddtheta_minus);

call_uniform_diff_matrices = true;
switch abs(preconditioner_streaming_theta_derivative_option)
    case 0
        fprintf('Streaming d/dtheta term is completely dropped in the preconditioner.\n')
        call_uniform_diff_matrices = false;
    case 100
        fprintf('Streaming d/dtheta term is the same in the preconditioner.\n')
        call_uniform_diff_matrices = false;
        streaming_ddtheta_plus_preconditioner = streaming_ddtheta_plus;
        streaming_ddtheta_minus_preconditioner = streaming_ddtheta_minus;
    case 2
        fprintf('Preconditioner streaming d/dtheta discretized using centered differences: 1 point on each side.\n')
        derivative_option_plus = 2;
        derivative_option_minus = derivative_option_plus;
    case 3
        fprintf('Preconditioner streaming d/dtheta discretized using centered differences: 2 points on each side.\n')
        derivative_option_plus = 12;
        derivative_option_minus = derivative_option_plus;
    case 4
        fprintf('Preconditioner streaming d/dtheta discretized using upwinded differences: 0 points on one side, 1 point on the other.\n')
        derivative_option_plus  = 32;
        derivative_option_minus = 42;
    case 5
        fprintf('Preconditioner streaming d/dtheta discretized using upwinded differences: 0 points on one side, 2 points on the other.\n')
        derivative_option_plus  = 52;
        derivative_option_minus = 62;
    case 6
        fprintf('Preconditioner streaming d/dtheta discretized using upwinded differences: 1 point on one side, 2 points on the other.\n')
        derivative_option_plus  = 82;
        derivative_option_minus = 92;
    case 7
        fprintf('Preconditioner streaming d/dtheta discretized using upwinded differences: 1 point on one side, 3 points on the other.\n')
        derivative_option_plus  = 102;
        derivative_option_minus = 112;
    case 8
        fprintf('Preconditioner streaming d/dtheta discretized using upwinded differences: 2 points on one side, 3 points on the other.\n')
        derivative_option_plus  = 122;
        derivative_option_minus = 132;
    otherwise
        error('Invalid preconditioner_streaming_theta_derivative_option: %d',preconditioner_streaming_theta_derivative_option)
end
if call_uniform_diff_matrices
    quadrature_option = 0;
    [~, ~, streaming_ddtheta_plus_preconditioner, ~]  = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, derivative_option_plus, quadrature_option);
    [~, ~, streaming_ddtheta_minus_preconditioner, ~] = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, derivative_option_minus, quadrature_option);
end

if discretizationParameters.theta_derivative_option<0
    fprintf('  But only the diagonal is kept.\n')
    streaming_ddtheta_plus_preconditioner = diag(diag(streaming_ddtheta_plus_preconditioner));
    streaming_ddtheta_minus_preconditioner = diag(diag(streaming_ddtheta_minus_preconditioner));
end

streaming_ddtheta_sum_preconditioner        = 0.5*(streaming_ddtheta_plus_preconditioner+streaming_ddtheta_minus_preconditioner);
streaming_ddtheta_difference_preconditioner = 0.5*(streaming_ddtheta_plus_preconditioner-streaming_ddtheta_minus_preconditioner);




% *************************************************************************
% Generate abscissae, quadrature weights, and derivative matrix for zeta grid.
% *************************************************************************
zetaMax = 2*pi/NPeriods;

if Nzeta==1
    zeta=0;
    zetaWeights=2*pi;
    buffer_zeta_points_on_each_side = 0;
    streaming_ddzeta_sum=0;
    streaming_ddzeta_difference=0;
    streaming_ddzeta_sum_preconditioner=0;
    streaming_ddzeta_difference_preconditioner=0;
else
    
    switch streaming_zeta_derivative_option
        case 2
            fprintf('Streaming d/dzeta discretized using centered differences: 1 point on each side.\n')
            derivative_option_plus = 2;
            derivative_option_minus = derivative_option_plus;
            buffer_zeta_points_on_each_side = 1;
        case 3
            fprintf('Streaming d/dzeta discretized using centered differences: 2 points on each side.\n')
            derivative_option_plus = 12;
            derivative_option_minus = derivative_option_plus;
            buffer_zeta_points_on_each_side = 2;
        case 4
            fprintf('Streaming d/dzeta discretized using upwinded differences: 0 points on one side, 1 point on the other.\n')
            derivative_option_plus  = 32;
            derivative_option_minus = 42;
            buffer_zeta_points_on_each_side = 1;
        case 5
            fprintf('Streaming d/dzeta discretized using upwinded differences: 0 points on one side, 2 points on the other.\n')
            derivative_option_plus  = 52;
            derivative_option_minus = 62;
            buffer_zeta_points_on_each_side = 2;
        case 6
            fprintf('Streaming d/dzeta discretized using upwinded differences: 1 point on one side, 2 points on the other.\n')
            derivative_option_plus  = 82;
            derivative_option_minus = 92;
            buffer_zeta_points_on_each_side = 2;
        case 7
            fprintf('Streaming d/dzeta discretized using upwinded differences: 1 point on one side, 3 points on the other.\n')
            derivative_option_plus  = 102;
            derivative_option_minus = 112;
            buffer_zeta_points_on_each_side = 3;
        case 8
            fprintf('Streaming d/dzeta discretized using upwinded differences: 2 points on one side, 3 points on the other.\n')
            derivative_option_plus  = 122;
            derivative_option_minus = 132;
            buffer_zeta_points_on_each_side = 3;
        otherwise
            error('Invalid streaming_zeta_derivative_option: %d',streaming_zeta_derivative_option)
    end
    Delta_zeta = (2*pi)/(NPeriods*(Nzeta-2*buffer_zeta_points_on_each_side));
    quadrature_option = 0;
    [zeta, ~, streaming_ddzeta_plus, ~]  = sfincs_uniformDiffMatrices(Nzeta, ...
        -buffer_zeta_points_on_each_side*Delta_zeta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta_zeta, derivative_option_plus, quadrature_option);
    [~   , ~, streaming_ddzeta_minus, ~] = sfincs_uniformDiffMatrices(Nzeta, ...
        -buffer_zeta_points_on_each_side*Delta_zeta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta_zeta, derivative_option_minus, quadrature_option);
    assert(abs(zeta(2)-zeta(1)-Delta_zeta)<1e-12)
    
    streaming_ddzeta_sum        = 0.5*(streaming_ddzeta_plus+streaming_ddzeta_minus);
    streaming_ddzeta_difference = 0.5*(streaming_ddzeta_plus-streaming_ddzeta_minus);
    
    call_uniform_diff_matrices = true;
    switch abs(preconditioner_streaming_zeta_derivative_option)
        case 0
            fprintf('Streaming d/dzeta term is completely dropped in the preconditioner.\n')
            call_uniform_diff_matrices = false;
        case 100
            fprintf('Streaming d/dzeta term is the same in the preconditioner.\n')
            call_uniform_diff_matrices = false;
            streaming_ddzeta_plus_preconditioner = streaming_ddzeta_plus;
            streaming_ddzeta_minus_preconditioner = streaming_ddzeta_minus;
        case 2
            fprintf('Preconditioner streaming d/dzeta discretized using centered differences: 1 point on each side.\n')
            derivative_option_plus = 2;
            derivative_option_minus = derivative_option_plus;
        case 3
            fprintf('Preconditioner streaming d/dzeta discretized using centered differences: 2 points on each side.\n')
            derivative_option_plus = 12;
            derivative_option_minus = derivative_option_plus;
        case 4
            fprintf('Preconditioner streaming d/dzeta discretized using upwinded differences: 0 points on one side, 1 point on the other.\n')
            derivative_option_plus  = 32;
            derivative_option_minus = 42;
        case 5
            fprintf('Preconditioner streaming d/dzeta discretized using upwinded differences: 0 points on one side, 2 points on the other.\n')
            derivative_option_plus  = 52;
            derivative_option_minus = 62;
        case 6
            fprintf('Preconditioner streaming d/dzeta discretized using upwinded differences: 1 point on one side, 2 points on the other.\n')
            derivative_option_plus  = 82;
            derivative_option_minus = 92;
        case 7
            fprintf('Preconditioner streaming d/dzeta discretized using upwinded differences: 1 point on one side, 3 points on the other.\n')
            derivative_option_plus  = 102;
            derivative_option_minus = 112;
        case 8
            fprintf('Preconditioner streaming d/dzeta discretized using upwinded differences: 2 points on one side, 3 points on the other.\n')
            derivative_option_plus  = 122;
            derivative_option_minus = 132;
        otherwise
            error('Invalid preconditioner_streaming_zeta_derivative_option: %d',preconditioner_streaming_zeta_derivative_option)
    end
    if call_uniform_diff_matrices
        quadrature_option = 0;
        [~, ~, streaming_ddzeta_plus_preconditioner, ~]  = sfincs_uniformDiffMatrices(Nzeta, ...
            -buffer_zeta_points_on_each_side*Delta_zeta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta_zeta, derivative_option_plus, quadrature_option);
        [~, ~, streaming_ddzeta_minus_preconditioner, ~] = sfincs_uniformDiffMatrices(Nzeta, ...
            -buffer_zeta_points_on_each_side*Delta_zeta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta_zeta, derivative_option_minus, quadrature_option);
    end
    
    if discretizationParameters.zeta_derivative_option<0
        fprintf('  But only the diagonal is kept.\n')
        streaming_ddzeta_plus_preconditioner = diag(diag(streaming_ddzeta_plus_preconditioner));
        streaming_ddzeta_minus_preconditioner = diag(diag(streaming_ddzeta_minus_preconditioner));
    end
    
    streaming_ddzeta_sum_preconditioner        = 0.5*(streaming_ddzeta_plus_preconditioner+streaming_ddzeta_minus_preconditioner);
    streaming_ddzeta_difference_preconditioner = 0.5*(streaming_ddzeta_plus_preconditioner-streaming_ddzeta_minus_preconditioner);
    
    
    
    zetaWeights=ones(size(zeta));
    zetaWeights(1:buffer_zeta_points_on_each_side)         = 0;
    zetaWeights(end-buffer_zeta_points_on_each_side+1:end) = 0;
    zetaWeights = zetaWeights * Delta * geometryParameters.Nperiods;
    assert(abs(sum(zetaWeights)-2*pi) < 1e-12)
    
    zeta_to_impose_DKE = (buffer_zeta_points_on_each_side+1):(Nzeta-buffer_zeta_points_on_each_side);

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
