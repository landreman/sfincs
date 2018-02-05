function sfincs_createGrids()

global Nalpha Nzeta Nx Nxi NL NxPotentialsPerVth xMax solverTolerance xGrid_k
global constraintScheme collisionOperator NPeriods pointAtX0
global forceOddNalphaAndNzeta buffer_zeta_points_on_each_side
global alpha alphaWeights zeta zetaWeights xi xiWeights xi_quadrature_option
global x xWeights ddx d2dx2 ddx_preconditioner xGridScheme xPotentialsGridScheme
global alpha2D zeta2D RHSMode preconditioner_x xMaxPotentials
global xInterpolationScheme xPotentialsInterpolationScheme NxPotentials
global xPotentials ddxPotentials d2dx2Potentials interpolateXToXPotentials
global indexVars Nspecies includePhi1 transportMatrix
global BLOCK_F BLOCK_QN BLOCK_PHI1_CONSTRAINT BLOCK_DENSITY_CONSTRAINT BLOCK_PRESSURE_CONSTRAINT BLOCK_F_CONSTRAINT
global zeta_to_impose_DKE zetaMax

global alpha_derivative_option preconditioner_alpha_derivative_option
global ddalpha_plus ddalpha_minus ddalpha_plus_preconditioner ddalpha_minus_preconditioner

global zeta_derivative_option preconditioner_zeta_derivative_option
global ddzeta_plus ddzeta_minus ddzeta_plus_preconditioner ddzeta_minus_preconditioner

global xi_derivative_option preconditioner_xi_derivative_option
global ddxi_plus ddxi_minus ddxi_plus_preconditioner ddxi_minus_preconditioner

global pitch_angle_scattering_option preconditioner_pitch_angle_scattering_option
global pitch_angle_scattering_operator pitch_angle_scattering_operator_preconditioner

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

switch alpha_derivative_option
    case 1
        fprintf('d/dalpha derivatives discretized using Fourier pseudospectral method.\n')
        derivative_option_plus = 20;
        derivative_option_minus = derivative_option_plus;
    case 2
        fprintf('d/dalpha derivatives discretized using centered finite differences: 1 point on either side.\n')
        derivative_option_plus = 0;
        derivative_option_minus = derivative_option_plus;
    case 3
        fprintf('d/dalpha derivatives discretized using centered finite differences: 2 points on either side.\n')
        derivative_option_plus = 10;
        derivative_option_minus = derivative_option_plus;
    case 4
        fprintf('d/dalpha derivatives discretized using upwinded finite differences: 0 points on one side, 1 point on the other side.\n')
        derivative_option_plus  = 30;
        derivative_option_minus = 40;
    case 5
        fprintf('d/dalpha derivatives discretized using upwinded finite differences: 0 points on one side, 2 point on the other side.\n')
        derivative_option_plus  = 50;
        derivative_option_minus = 60;
    case 6
        fprintf('d/dalpha derivatives discretized using upwinded finite differences: 1 points on one side, 2 point on the other side.\n')
        derivative_option_plus  = 80;
        derivative_option_minus = 90;
    case 7
        fprintf('d/dalpha derivatives discretized using upwinded finite differences: 1 point on one side, 3 point on the other side.\n')
        derivative_option_plus  = 100;
        derivative_option_minus = 110;
    case 8
        fprintf('d/dalpha derivatives discretized using upwinded finite differences: 2 points on one side, 3 points on the other side.\n')
        derivative_option_plus  = 120;
        derivative_option_minus = 130;
    otherwise
        error('Error! Invalid alpha_derivative_option')
end
quadrature_option = 0;
[alpha, alphaWeights, ddalpha_plus, ~]  = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, derivative_option_plus,  quadrature_option);
[~, ~, ddalpha_minus, ~] = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, derivative_option_minus, quadrature_option);

call_uniform_diff_matrices = true;
switch abs(preconditioner_alpha_derivative_option)
    case 0
        call_uniform_diff_matrices = false;
        fprintf('d/dalpha derivatives are completely dropped in the preconditioner.\n')
        ddalpha_plus_preconditioner = zeros(Nalpha);
        ddalpha_minus_preconditioner = zeros(Nalpha);
    case 100
        call_uniform_diff_matrices = false;
        fprintf('d/dalpha derivatives are the same in the preconditioner as in the main matrix.\n')
        ddalpha_plus_preconditioner = ddalpha_plus;
        ddalpha_minus_preconditioner = ddalpha_minus;
    case 1
        fprintf('Preconditioner d/dalpha derivatives discretized using Fourier pseudospectral method.\n')
        derivative_option_plus = 20;
        derivative_option_minus = derivative_option_plus;
    case 2
        fprintf('Preconditioner d/dalpha derivatives discretized using centered finite differences: 1 point on either side.\n')
        derivative_option_plus = 0;
        derivative_option_minus = derivative_option_plus;
    case 3
        fprintf('Preconditioner d/dalpha derivatives discretized using centered finite differences: 2 points on either side.\n')
        derivative_option_plus = 10;
        derivative_option_minus = derivative_option_plus;
    case 4
        fprintf('Preconditioner d/dalpha derivatives discretized using upwinded finite differences: 0 points on one side, 1 point on the other side.\n')
        derivative_option_plus  = 30;
        derivative_option_minus = 40;
    case 5
        fprintf('Preconditioner d/dalpha derivatives discretized using upwinded finite differences: 0 points on one side, 2 point on the other side.\n')
        derivative_option_plus  = 50;
        derivative_option_minus = 60;
    case 6
        fprintf('Preconditioner d/dalpha derivatives discretized using upwinded finite differences: 1 points on one side, 2 point on the other side.\n')
        derivative_option_plus  = 80;
        derivative_option_minus = 90;
    case 7
        fprintf('Preconditioner d/dalpha derivatives discretized using upwinded finite differences: 1 point on one side, 3 point on the other side.\n')
        derivative_option_plus  = 100;
        derivative_option_minus = 110;
    case 8
        fprintf('Preconditioner d/dalpha derivatives discretized using upwinded finite differences: 2 points on one side, 3 points on the other side.\n')
        derivative_option_plus  = 120;
        derivative_option_minus = 130;
    otherwise
        error('Error! Invalid preconditioner_alpha_derivative_option')
end
if call_uniform_diff_matrices
    quadrature_option = 0;
    [~, ~, ddalpha_plus_preconditioner, ~]  = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, derivative_option_plus,  quadrature_option);
    [~, ~, ddalpha_minus_preconditioner, ~] = sfincs_uniformDiffMatrices(Nalpha, 0, 2*pi, derivative_option_minus, quadrature_option);
end
if preconditioner_alpha_derivative_option<0
    fprintf('   But only the diagonal is kept.\n')
    ddalpha_plus_preconditioner = diag(diag(ddalpha_plus_preconditioner));
    ddalpha_minus_preconditioner = diag(diag(ddalpha_minus_preconditioner));
end


% *************************************************************************
% Generate abscissae, quadrature weights, and derivative matrix for zeta grid.
% *************************************************************************
zetaMax = 2*pi/NPeriods;

if Nzeta==1
    zeta=0;
    zetaWeights=2*pi;
    ddzeta_plus=0;
    ddzeta_minus=0;
    ddzeta_plus_preconditioner=0;
    ddzeta_minus_preconditioner=0;
    buffer_zeta_points_on_each_side = 0;
    zeta_to_impose_DKE = 1;
else
    switch zeta_derivative_option
        case 2
            fprintf('d/dzeta derivatives discretized using centered finite differences: 1 point on either side.\n')
            derivative_option_plus = 2;
            derivative_option_minus = derivative_option_plus;
            buffer_zeta_points_on_each_side = 1;
        case 3
            fprintf('d/dzeta derivatives discretized using centered finite differences: 2 points on either side.\n')
            derivative_option_plus = 12;
            derivative_option_minus = derivative_option_plus;
            buffer_zeta_points_on_each_side = 2;
        case 4
            fprintf('d/dzeta derivatives discretized using upwinded finite differences: 0 points on one side, 1 point on the other side.\n')
            derivative_option_plus  = 32;
            derivative_option_minus = 42;
            buffer_zeta_points_on_each_side = 1;
        case 5
            fprintf('d/dzeta derivatives discretized using upwinded finite differences: 0 points on one side, 2 point on the other side.\n')
            derivative_option_plus  = 52;
            derivative_option_minus = 62;
            buffer_zeta_points_on_each_side = 2;
        case 6
            fprintf('d/dzeta derivatives discretized using upwinded finite differences: 1 points on one side, 2 point on the other side.\n')
            derivative_option_plus  = 82;
            derivative_option_minus = 92;
            buffer_zeta_points_on_each_side = 2;
        case 7
            fprintf('d/dzeta derivatives discretized using upwinded finite differences: 1 point on one side, 3 point on the other side.\n')
            derivative_option_plus  = 102;
            derivative_option_minus = 112;
            buffer_zeta_points_on_each_side = 3;
        case 8
            fprintf('d/dzeta derivatives discretized using upwinded finite differences: 2 points on one side, 3 points on the other side.\n')
            derivative_option_plus  = 122;
            derivative_option_minus = 132;
            buffer_zeta_points_on_each_side = 3;
        otherwise
            error('Error! Invalid zeta_derivative_option')
    end
    Delta_zeta = (2*pi)/(NPeriods*(Nzeta-2*buffer_zeta_points_on_each_side));
    quadrature_option = 0;
    [zeta, ~, ddzeta_plus, ~]  = sfincs_uniformDiffMatrices(Nzeta, -buffer_zeta_points_on_each_side*Delta_zeta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta_zeta, derivative_option_plus,  quadrature_option);
    [~, ~, ddzeta_minus, ~] = sfincs_uniformDiffMatrices(Nzeta, -buffer_zeta_points_on_each_side*Delta_zeta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta_zeta, derivative_option_minus, quadrature_option);
    assert(abs(zeta(2)-zeta(1) - Delta_zeta) < 1e-12)
    
    zetaWeights=ones(Nzeta,1);
    zetaWeights(1:buffer_zeta_points_on_each_side)      =0;
    zetaWeights(end-buffer_zeta_points_on_each_side+1:end)=0;
    zetaWeights = zetaWeights * Delta_zeta * NPeriods;
    assert(abs(sum(zetaWeights)-2*pi)<1e-12)
    
    zeta_to_impose_DKE = (1+buffer_zeta_points_on_each_side):(Nzeta-buffer_zeta_points_on_each_side);
    
    call_uniform_diff_matrices = true;
    switch abs(preconditioner_zeta_derivative_option)
        case 0
            call_uniform_diff_matrices = false;
            fprintf('d/dzeta derivatives are completely dropped in the preconditioner.\n')
            ddzeta_plus_preconditioner = zeros(Nzeta);
            ddzeta_minus_preconditioner = zeros(Nzeta);
        case 100
            call_uniform_diff_matrices = false;
            fprintf('d/dzeta derivatives are the same in the preconditioner as in the main matrix.\n')
            ddzeta_plus_preconditioner = ddzeta_plus;
            ddzeta_minus_preconditioner = ddzeta_minus;
        case 2
            fprintf('Preconditioner d/dzeta derivatives discretized using centered finite differences: 1 point on either side.\n')
            derivative_option_plus = 2;
            derivative_option_minus = derivative_option_plus;
        case 3
            fprintf('Preconditioner d/dzeta derivatives discretized using centered finite differences: 2 points on either side.\n')
            derivative_option_plus = 12;
            derivative_option_minus = derivative_option_plus;
        case 4
            fprintf('Preconditioner d/dzeta derivatives discretized using upwinded finite differences: 0 points on one side, 1 point on the other side.\n')
            derivative_option_plus  = 32;
            derivative_option_minus = 42;
        case 5
            fprintf('Preconditioner d/dzeta derivatives discretized using upwinded finite differences: 0 points on one side, 2 point on the other side.\n')
            derivative_option_plus  = 52;
            derivative_option_minus = 62;
        case 6
            fprintf('Preconditioner d/dzeta derivatives discretized using upwinded finite differences: 1 points on one side, 2 point on the other side.\n')
            derivative_option_plus  = 82;
            derivative_option_minus = 92;
        case 7
            fprintf('Preconditioner d/dzeta derivatives discretized using upwinded finite differences: 1 point on one side, 3 point on the other side.\n')
            derivative_option_plus  = 102;
            derivative_option_minus = 112;
        case 8
            fprintf('Preconditioner d/dzeta derivatives discretized using upwinded finite differences: 2 points on one side, 3 points on the other side.\n')
            derivative_option_plus  = 122;
            derivative_option_minus = 132;
        otherwise
            error('Error! Invalid preconditioner_zeta_derivative_option')
    end
    if call_uniform_diff_matrices
        quadrature_option = 0;
        [~, ~, ddzeta_plus_preconditioner, ~]  = sfincs_uniformDiffMatrices(Nzeta, -buffer_zeta_points_on_each_side*Delta_zeta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta_zeta, derivative_option_plus,  quadrature_option);
        [~, ~, ddzeta_minus_preconditioner, ~] = sfincs_uniformDiffMatrices(Nzeta, -buffer_zeta_points_on_each_side*Delta_zeta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta_zeta, derivative_option_minus, quadrature_option);
    end
    if preconditioner_zeta_derivative_option<0
        fprintf('   But only the diagonal is kept.\n')
        ddzeta_plus_preconditioner = diag(diag(ddzeta_plus_preconditioner));
        ddzeta_minus_preconditioner = diag(diag(ddzeta_minus_preconditioner));
    end
    
    %{
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
    %}
    
    
end

% *************************************************************************
% Generate xi grid and pitch angle scattering operator
% *************************************************************************


switch pitch_angle_scattering_option
    case 2
        fprintf('Pitch angle scattering operator discretized using centered finite differences: 1 point on either side.\n')
        derivative_option = 2;
    case 3
        fprintf('Pitch angle scattering operator discretized using centered finite differences: 2 points on either side.\n')
        derivative_option = 12;
    otherwise
        error('Error! Invalid pitch_angle_scattering_option')
end
[xi, xiWeights, ddxi, d2dxi2]  = sfincs_uniformDiffMatrices(Nxi, -1, 1, derivative_option,  xi_quadrature_option);
pitch_angle_scattering_operator = 0.5*diag(1-xi.*xi)*d2dxi2 - diag(xi)*ddxi;
assert(abs(sum(xiWeights)-2)<1e-12)

call_uniform_diff_matrices = true;
switch abs(preconditioner_pitch_angle_scattering_option)
    case 0
        fprintf('Pitch angle scattering is completely dropped in the preconditioner.\n')
        call_uniform_diff_matrices = false;
        pitch_angle_scattering_operator_preconditioner = zeros(Nxi);
    case 100
        fprintf('Pitch angle scattering operator is the same in preconditioner as in the main matrix.\n')
        call_uniform_diff_matrices = false;
        pitch_angle_scattering_operator_preconditioner = pitch_angle_scattering_operator;
    case 2
        fprintf('Preconditioner pitch angle scattering operator discretized using centered finite differences: 1 point on either side.\n')
        derivative_option = 2;
    case 3
        fprintf('Preconditioner pitch angle scattering operator discretized using centered finite differences: 2 points on either side.\n')
        derivative_option = 12;
    otherwise
        error('Error! Invalid preconditioner_pitch_angle_scattering_option')
end
if call_uniform_diff_matrices
    [~, ~, ddxi, d2dxi2]  = sfincs_uniformDiffMatrices(Nxi, -1, 1, derivative_option,  xi_quadrature_option);
    pitch_angle_scattering_operator_preconditioner = 0.5*diag(1-xi.*xi)*d2dxi2 - diag(xi)*ddxi;
end

if preconditioner_pitch_angle_scattering_option<0
    fprintf('   But only the diagonal is kept.\n')
    pitch_angle_scattering_operator_preconditioner = diag(diag(pitch_angle_scattering_operator_preconditioner));
end

% *************************************************************************
% Generate differentiation matrices for the collisionless df/dxi term
% *************************************************************************

switch xi_derivative_option
    case 2
        fprintf('d/dxi derivatives discretized using centered finite differences: 1 point on either side.\n')
        derivative_option_plus = 2;
        derivative_option_minus = derivative_option_plus;
    case 3
        fprintf('d/dxi derivatives discretized using centered finite differences: 2 points on either side.\n')
        derivative_option_plus = 12;
        derivative_option_minus = derivative_option_plus;
    case 4
        fprintf('d/dxi derivatives discretized using upwinded finite differences: 0 points on one side, 1 point on the other side.\n')
        derivative_option_plus  = 32;
        derivative_option_minus = 42;
    case 5
        fprintf('d/dxi derivatives discretized using upwinded finite differences: 0 points on one side, 2 point on the other side.\n')
        derivative_option_plus  = 52;
        derivative_option_minus = 62;
    case 6
        fprintf('d/dxi derivatives discretized using upwinded finite differences: 1 points on one side, 2 point on the other side.\n')
        derivative_option_plus  = 82;
        derivative_option_minus = 92;
    case 7
        fprintf('d/dxi derivatives discretized using upwinded finite differences: 1 point on one side, 3 point on the other side.\n')
        derivative_option_plus  = 102;
        derivative_option_minus = 112;
    case 8
        fprintf('d/dxi derivatives discretized using upwinded finite differences: 2 points on one side, 3 points on the other side.\n')
        derivative_option_plus  = 122;
        derivative_option_minus = 132;
    otherwise
        error('Error! Invalid xi_derivative_option')
end
quadrature_option = 0;
[~, ~, ddxi_plus, ~]  = sfincs_uniformDiffMatrices(Nxi, -1, 1, derivative_option_plus,  quadrature_option);
[~, ~, ddxi_minus, ~] = sfincs_uniformDiffMatrices(Nxi, -1, 1, derivative_option_minus, quadrature_option);

call_uniform_diff_matrices = true;
switch abs(preconditioner_xi_derivative_option)
    case 0
        call_uniform_diff_matrices = false;
        fprintf('d/dxi derivatives are completely dropped in the preconditioner.\n')
        ddxi_plus_preconditioner = zeros(Nxi);
        ddxi_minus_preconditioner = zeros(Nxi);
    case 100
        call_uniform_diff_matrices = false;
        fprintf('d/dxi derivatives are the same in the preconditioner as in the main matrix.\n')
        ddxi_plus_preconditioner = ddxi_plus;
        ddxi_minus_preconditioner = ddxi_minus;
    case 2
        fprintf('Preconditioner d/dxi derivatives discretized using centered finite differences: 1 point on either side.\n')
        derivative_option_plus = 2;
        derivative_option_minus = derivative_option_plus;
    case 3
        fprintf('Preconditioner d/dxi derivatives discretized using centered finite differences: 2 points on either side.\n')
        derivative_option_plus = 12;
        derivative_option_minus = derivative_option_plus;
    case 4
        fprintf('Preconditioner d/dxi derivatives discretized using upwinded finite differences: 0 points on one side, 1 point on the other side.\n')
        derivative_option_plus  = 32;
        derivative_option_minus = 42;
    case 5
        fprintf('Preconditioner d/dxi derivatives discretized using upwinded finite differences: 0 points on one side, 2 point on the other side.\n')
        derivative_option_plus  = 52;
        derivative_option_minus = 62;
    case 6
        fprintf('Preconditioner d/dxi derivatives discretized using upwinded finite differences: 1 points on one side, 2 point on the other side.\n')
        derivative_option_plus  = 82;
        derivative_option_minus = 92;
    case 7
        fprintf('Preconditioner d/dxi derivatives discretized using upwinded finite differences: 1 point on one side, 3 point on the other side.\n')
        derivative_option_plus  = 102;
        derivative_option_minus = 112;
    case 8
        fprintf('Preconditioner d/dxi derivatives discretized using upwinded finite differences: 2 points on one side, 3 points on the other side.\n')
        derivative_option_plus  = 122;
        derivative_option_minus = 132;
    otherwise
        error('Error! Invalid preconditioner_xi_derivative_option')
end
if call_uniform_diff_matrices
    quadrature_option = 0;
    [~, ~, ddxi_plus_preconditioner, ~]  = sfincs_uniformDiffMatrices(Nxi, -1, 1, derivative_option_plus,  quadrature_option);
    [~, ~, ddxi_minus_preconditioner, ~] = sfincs_uniformDiffMatrices(Nxi, -1, 1, derivative_option_minus, quadrature_option);
end
if preconditioner_xi_derivative_option<0
    fprintf('   But only the diagonal is kept.\n')
    ddxi_plus_preconditioner = diag(diag(ddxi_plus_preconditioner));
    ddxi_minus_preconditioner = diag(diag(ddxi_minus_preconditioner));
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
    quadrature_option = 0;
    [xPotentials, ~, ddxPotentials, d2dx2Potentials] = sfincs_uniformDiffMatrices(NxPotentials, 0, xMaxPotentials, scheme, quadrature_option);
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
