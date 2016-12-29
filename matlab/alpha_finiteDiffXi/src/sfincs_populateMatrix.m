function matrix = sfincs_populateMatrix(whichMatrix)

global THats mHats Zs nHats RosenbluthPotentialTerms xGridScheme
global matrixSize Nalpha Nzeta Nxi Nx NL Nspecies pointAtX0
global ddalpha_plus ddalpha_minus ddalpha_plus_preconditioner ddalpha_minus_preconditioner
global ddzeta_plus ddzeta_minus ddzeta_plus_preconditioner ddzeta_minus_preconditioner
global ddxi_plus ddxi_minus ddxi_plus_preconditioner ddxi_minus_preconditioner
global pitch_angle_scattering_operator pitch_angle_scattering_operator_preconditioner
global ddx ddx_preconditioner x xWeights xMaxPotentials xi
global NxPotentials xPotentials ddxPotentials d2dx2Potentials interpolateXToXPotentials
global alphaWeights zetaWeights xiWeights
global preconditioner_species  
global collisionOperator constraintScheme
global BLOCK_F BLOCK_QN BLOCK_PHI1_CONSTRAINT BLOCK_DENSITY_CONSTRAINT BLOCK_PRESSURE_CONSTRAINT BLOCK_F_CONSTRAINT indexVars
global includePhi1 includePhi1InKineticEquation ExB_option includeXDotTerm includeElectricFieldTermInXiDot magneticDriftScheme
global Delta gamma nu_n Phi1Hat dPhi1Hatdalpha dPhi1Hatdzeta stateVector
global BDotCurlB FSABHat2 dPhiHatdpsiHat force0RadialCurrentInEquilibrium
global BHat dBHatdtheta dBHatdzeta dBHatdpsiHat GHat B0OverBBar VPrimeHat sqrt_g_sign
global DHat BHat_sub_psi BHat_sub_theta BHat_sub_zeta BHat_sup_theta BHat_sup_zeta
global dBHat_sub_psi_dtheta dBHat_sub_psi_dzeta
global dBHat_sub_theta_dpsiHat dBHat_sub_theta_dzeta
global dBHat_sub_zeta_dpsiHat dBHat_sub_zeta_dtheta
global adiabaticZ adiabaticNHat adiabaticTHat withAdiabatic quasineutralityOption
global dnHatdpsiHats dTHatdpsiHats reusePreconditioner
global zeta_to_impose_DKE zetaMax buffer_zeta_points_on_each_side
global alpha_interpolation_stencil preconditioner_alpha_interpolation_stencil alpha iota

populateMatrixTic = tic;

if pointAtX0
    ixMin = 2;
else
    ixMin = 1;
end

x2 = x.*x;
expx2 = exp(-x2);
sqrtpi = sqrt(pi);

switch whichMatrix
    case 0
        whichMatrixName = 'Jacobian preconditioner';
    case 1
        whichMatrixName = 'Jacobian';
    case 2
        whichMatrixName = 'Residual f0';
    case 3
        whichMatrixName = 'Residual f1';
    otherwise
        error('Invalid whichMatrix')
end

fprintf('Populating matrix: %s\n',whichMatrixName)

% To build the matrix as efficiently as possible, a reasonably
% accurate estimate of the number of nonzeros (nnz) is needed beforehand:
estimated_nnz = 1 * (Nx*Nx*Nspecies*Nspecies*Nxi*Nalpha*Nzeta ...
    + Nspecies*(nnz(ddalpha_plus)*Nx*(3*Nxi)*Nzeta + nnz(ddzeta_plus)*Nx*(3*Nxi)*Nalpha ...
    + Nx*(5*Nxi)*Nalpha*Nzeta + 3*Nx*Nx*Nxi*Nalpha*Nzeta ...
    + 2*2*Nx*1*Nalpha*Nzeta));

estimated_nnz = Nspecies*(Nalpha*5)*Nzeta*(Nxi*5)*Nx ...  %ddalpha terms
    + Nspecies*Nalpha*(Nzeta*5)*(Nxi*5)*Nx ...            %ddzeta terms
    + Nspecies*Nalpha*Nzeta*(Nxi*5)*Nx ...                %ddxi terms
    + Nspecies*Nalpha*Nzeta*(Nxi*5)*(Nx*Nx) ...           %ddx terms (collisionless)
    + Nspecies*Nalpha*Nzeta*(Nxi*5)*(Nx*Nx);              %collision terms

if constraintScheme==1
elseif constraintScheme==2
end

if includePhi1
    estimated_nnz = estimated_nnz + Nalpha*Nzeta*Nx ...  % quasineutrality equation
        + (Nalpha*5)*Nzeta*Nxi*Nx ...                    % dPhi1/dalpha
        + Nalpha*(Nzeta*5)*Nxi*Nx;                    % dPhi1/dzeta
end

sparseCreatorIndex=1;
sparseCreator_i=0;
sparseCreator_j=0;
sparseCreator_s=0;
resetSparseCreator()

% -----------------------------------------
% Add interpolations:
% -----------------------------------------

if whichMatrix==0
    stencil = preconditioner_alpha_interpolation_stencil;
else
    stencil = alpha_interpolation_stencil;
end

% First handle the points needed by the DKE to the left of zeta=0:
theta_left = alpha;

theta_right = alpha - iota * zetaMax;
interpolationMatrix_left = sfincs_periodicInterpolation(theta_left, theta_right, 2*pi, stencil);

theta_right = alpha + iota * zetaMax;
interpolationMatrix_right = sfincs_periodicInterpolation(theta_left, theta_right, 2*pi, stencil);

izetas_left = 1:buffer_zeta_points_on_each_side;
izetas_right = (Nzeta-buffer_zeta_points_on_each_side+1):Nzeta;
izeta_shift = Nzeta - 2*buffer_zeta_points_on_each_side;
%{
if zetaDerivativeScheme==1
    izetas_left = 1;
    izetas_right = Nzeta;
    izeta_shift = Nzeta-2;
elseif zetaDerivativeScheme==2
    izetas_left = [1,2];
    izetas_right = [Nzeta-1,Nzeta];
    izeta_shift = Nzeta-4;
else
    error('Should not get here')
end
%}
for ispecies = 1:Nspecies
    for ix = 1:Nx
        for ixi = 1:Nxi
            for izeta = izetas_left
                % Add 1's along the diagonal
                rowIndices = sfincs_indices(ispecies, ix, ixi, 1:Nalpha, izeta, BLOCK_F, indexVars);
                addToSparse(rowIndices, rowIndices, ones(size(rowIndices)));
                colIndices = sfincs_indices(ispecies, ix, ixi, 1:Nalpha, izeta+izeta_shift, BLOCK_F, indexVars);
                addSparseBlock(rowIndices,colIndices,-interpolationMatrix_left);
            end
            
            for izeta = izetas_right
                % Add 1's along the diagonal
                rowIndices = sfincs_indices(ispecies, ix, ixi, 1:Nalpha, izeta, BLOCK_F, indexVars);
                addToSparse(rowIndices, rowIndices, ones(size(rowIndices)));
                colIndices = sfincs_indices(ispecies, ix, ixi, 1:Nalpha, izeta-izeta_shift, BLOCK_F, indexVars);
                addSparseBlock(rowIndices,colIndices,-interpolationMatrix_right);
            end
        end
    end
end

% -----------------------------------------
% Add collisionless terms:
% -----------------------------------------

for ispecies = 1:Nspecies
    THat = THats(ispecies);
    nHat = nHats(ispecies);
    mHat = mHats(ispecies);
    Z = Zs(ispecies);
    sqrtT = sqrt(THat);
    sqrtm = sqrt(mHat);
    
    % -----------------------------------------
    % Add d/dzeta terms:
    % -----------------------------------------
    
    if whichMatrix ~= 2
        if whichMatrix==0
            ddzeta_plus_to_use  = ddzeta_plus_preconditioner;
            ddzeta_minus_to_use = ddzeta_minus_preconditioner;
        else
            ddzeta_plus_to_use  = ddzeta_plus;
            ddzeta_minus_to_use = ddzeta_minus;
        end
        switch ExB_option
            case 0
                Er_term = zeros(Nalpha,Nzeta);
            case 1
                Er_term = - (sqrt_g_sign*gamma*Delta*sqrtm*dPhiHatdpsiHat*BHat_sub_theta) ./ (2*sqrtT*BHat);
            case 2
                Er_term = - (sqrt_g_sign*gamma*Delta*sqrtm*dPhiHatdpsiHat*BHat_sub_theta.*BHat) / (2*sqrtT*FSABHat2);
            otherwise
                error('Invalid ExB_option')
        end
        for izeta_row = zeta_to_impose_DKE
            for ialpha = 1:Nalpha
                for ix = 1:Nx
                    for ixi = 1:Nxi
                        factor = sqrt_g_sign*x(ix)*xi(ixi) + Er_term(ialpha,izeta_row);
                        if factor>0
                            stuff_to_add = factor*ddzeta_plus_to_use(izeta_row,:);
                        else
                            stuff_to_add = factor*ddzeta_minus_to_use(izeta_row,:);
                        end

                        rowIndex = sfincs_indices(ispecies, ix, ixi, ialpha, izeta_row, BLOCK_F, indexVars);
                        colIndices = sfincs_indices(ispecies, ix, ixi, ialpha, 1:Nzeta, BLOCK_F, indexVars);

                        addSparseBlock(rowIndex, colIndices, stuff_to_add)
    
                    end
                end
            end
        end
    end
    
    % -----------------------------------------
    % Add d/dalpha terms:
    % -----------------------------------------
    
    if whichMatrix ~= 2
        if whichMatrix==0
            ddalpha_plus_to_use  = ddalpha_plus_preconditioner;
            ddalpha_minus_to_use = ddalpha_minus_preconditioner;
        else
            ddalpha_plus_to_use  = ddalpha_plus;
            ddalpha_minus_to_use = ddalpha_minus;
        end
        if Nzeta==1
            % Tokamak
            izeta = 1;
            switch ExB_option
                case 0
                    Er_term = zeros(Nalpha,Nzeta);
                case 1
                    Er_term = (sqrt_g_sign*gamma*Delta*sqrtm*dPhiHatdpsiHat*BHat_sub_zeta) ./ (2*sqrtT*BHat);
                case 2
                    Er_term = (sqrt_g_sign*gamma*Delta*sqrtm*dPhiHatdpsiHat*BHat_sub_zeta.*BHat) / (2*sqrtT*FSABHat2);
                otherwise
                    error('Invalid ExB_option')
            end
            for ialpha_row = 1:Nalpha
                for ix = 1:Nx
                    for ixi = 1:Nxi
                        factor = sqrt_g_sign*iota*x(ix)*xi(ixi) + Er_term(ialpha_row,izeta);
                        if factor>0
                            stuff_to_add = factor*ddalpha_plus_to_use(ialpha_row,:);
                        else
                            stuff_to_add = factor*ddalpha_minus_to_use(ialpha_row,:);
                        end

                        rowIndex = sfincs_indices(ispecies, ix, ixi, ialpha_row, izeta, BLOCK_F, indexVars);
                        colIndices = sfincs_indices(ispecies, ix, ixi, 1:Nalpha, izeta, BLOCK_F, indexVars);

                        addSparseBlock(rowIndex, colIndices, stuff_to_add)
                    end
                end
            end
        else
            % Nzeta>1, so we are in a stellarator.
            switch ExB_option
                case 0
                    Er_term = zeros(Nalpha,Nzeta);
                case 1
                    Er_term = (gamma*Delta*sqrtm*dPhiHatdpsiHat*BHat) ./ (2*sqrtT*abs(DHat));
                case 2
                    Er_term = (gamma*Delta*sqrtm*dPhiHatdpsiHat*(BHat.^3)) ./ (2*sqrtT*abs(DHat)*FSABHat2);
                otherwise
                    error('Invalid ExB_option')
            end
            for izeta = zeta_to_impose_DKE
                for ialpha_row = 1:Nalpha
                    factor = Er_term(ialpha_row,izeta);
                    if factor>0
                        factor2 = factor * ddalpha_plus_to_use(ialpha_row,:);
                    else
                        factor2 = factor * ddalpha_minus_to_use(ialpha_row,:);
                    end
                    for ialpha_col = 1:Nalpha
                        stuff_to_add = factor2(ialpha_col) * ones(Nxi,1);
                        for ix = 1:Nx
                            rowIndices = sfincs_indices(ispecies, ix, 1:Nxi, ialpha_row, izeta, BLOCK_F, indexVars);
                            colIndices = sfincs_indices(ispecies, ix, 1:Nxi, ialpha_col, izeta, BLOCK_F, indexVars);
                            addToSparse(rowIndices, colIndices, stuff_to_add)
                        end
                    end
                end
            end
        end
    end
    
    % -----------------------------------------
    % Add df/dxi terms:
    % -----------------------------------------
    
    if whichMatrix ~= 2
        if whichMatrix==0
            ddxi_plus_to_use  = ddxi_plus_preconditioner;
            ddxi_minus_to_use = ddxi_minus_preconditioner;
        else
            ddxi_plus_to_use  = ddxi_plus;
            ddxi_minus_to_use = ddxi_minus;
        end
        
        mirrorTermSpatialPart = -(sqrt_g_sign*iota*dBHatdtheta + dBHatdzeta) ./ (2*BHat);

        if includeElectricFieldTermInXiDot
            Er_term = (sqrt_g_sign*gamma*Delta*sqrtm*dPhiHatdpsiHat*(BHat_sub_zeta.*dBHatdtheta - BHat_sub_theta.*dBHatdzeta)) ...
                ./ (4*sqrtT*BHat.*BHat);
        else
            Er_term = zeros(Nalpha,Nzeta);
        end
        
        for ialpha = 1:Nalpha
            for izeta = zeta_to_impose_DKE
                for ix = 1:Nx
                    for ixi_row = 1:Nxi
                        factor = 1-xi(ixi_row)*xi(ixi_row);
                        factor2 = mirrorTermSpatialPart(ialpha,izeta)*x(ix)*factor + xi(ixi_row)*factor*Er_term(ialpha,izeta);
                        rowIndex = sfincs_indices(ispecies, ix, ixi_row, ialpha, izeta, BLOCK_F, indexVars);
                        colIndices = sfincs_indices(ispecies, ix, 1:Nxi, ialpha, izeta, BLOCK_F, indexVars);
                        if factor2>0
                            factor3 = factor2 * ddxi_plus_to_use(ixi_row,:);
                        else
                            factor3 = factor2 * ddxi_minus_to_use(ixi_row,:);
                        end

                        addSparseBlock(rowIndex, colIndices, factor3)
                    end
                end
            end
        end

    end
    
    % -----------------------------------------
    % Add the collisionless df1/dx term:
    % -----------------------------------------
    
    if whichMatrix ~= 2 && includeXDotTerm
        if whichMatrix==0
            x_part_of_x_dot = diag(x)*ddx_preconditioner;
        else
            x_part_of_x_dot = diag(x)*ddx;
        end
        
        factor = -(sqrt_g_sign*gamma*Delta*sqrtm*dPhiHatdpsiHat* (BHat_sub_theta.*dBHatdzeta - BHat_sub_zeta.*dBHatdtheta))./(4*sqrtT*(BHat.^2));
        for ialpha = 1:Nalpha
            for izeta = zeta_to_impose_DKE
                for ixi = 1:Nxi
                    indices = sfincs_indices(ispecies, 1:Nx, ixi, ialpha, izeta, BLOCK_F, indexVars);
                    addSparseBlock(indices, indices, factor(ialpha,izeta)*(1+xi(ixi)*xi(ixi))*x_part_of_x_dot)
                end
            end
        end
    end
    
end

switch (collisionOperator)
    case 0
        error('This collision operator is not yet implemented')
        % Linearized Fokker-Planck operator
        
        xWith0s = [0, xPotentials(2:(end-1))', 0];
        M21 = 4*pi*diag(xWith0s.^2) * interpolateXToXPotentials;
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
                % The subtraction in the next line causes a loss of some
                % digits at small x. Is there a better method?
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
                    regridSpecies(:,:,speciesA,speciesB) = sfincs_polynomialInterpolationMatrix(x,xb,sfincs_xWeight(x),sfincs_xWeight(xb));
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
                if whichMatrix > 0
                    % We're not making the preconditioner.
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
                        if xGridScheme>=5
                            CHat = M11 + squeeze(RosenbluthPotentialTerms(speciesA,speciesB,L+1,:,:));
                        else
                            M13 = M13IncludingX0(:,:,speciesA, speciesB, L+1);
                            M12 = M12IncludingX0(:,:,speciesA, speciesB, L+1);
                            
                            % Add Dirichlet or Neumann boundary condition for
                            % potentials at x=0:
                            if L~=0
                                M12(:,1) = 0;
                                M13(:,1) = 0;
                            end
                            
                            CHat = M11 -  (M12 - M13 * M33BackslashM32) * M22BackslashM21;
                        end
                    else
                        CHat = M11;
                    end
                    
                    % The lines below are invoked to make the preconditioner.
                    if whichMatrix == 0 && L >= preconditioner_x_min_L
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
                                error('Invalid preconditioner_x')
                        end
                        
                    end
                    
                    if pointAtX0
                        CHat(1,:) = 0;
                        if L ~= 0
                            CHat(:,1) = 0;
                        end
                    end
                    
                    % At this point, CHat holds the collision operator
                    % divided by \bar{nu}
                    
                    for ialpha = 1:Nalpha
                        for izeta = zeta_to_impose_DKE
                            rowIndices = sfincs_indices(speciesA, 1:Nx, L+1, ialpha, izeta, BLOCK_F, indexVars);
                            colIndices = sfincs_indices(speciesB, 1:Nx, L+1, ialpha, izeta, BLOCK_F, indexVars);
                            addSparseBlock(rowIndices, colIndices, -nu_n*CHat)
                        end
                    end
                    
                    
                end
            end
            
        end
        % End of multi-species Fokker-Planck collision
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
                % The subtraction in the next line causes a loss of some
                % digits at small x. Is there a better method?
                Psi = (erfs - 2/sqrtpi*xb .* expxb2) ./ (2*xb2);
                nuD(:,speciesA) = nuD(:,speciesA) + (speciesFactorTest * (erfs - Psi) ./ (x.^3));
            end
        end
        
        if whichMatrix==0
            pitch_angle_scattering_operator_to_use = pitch_angle_scattering_operator_preconditioner;
        else
            pitch_angle_scattering_operator_to_use = pitch_angle_scattering_operator;
        end
        
        for iSpecies = 1:Nspecies
            spatial_part = - (nu_n * BHat * sqrt(mHats(iSpecies))) ./ (abs(DHat) * sqrt(THats(iSpecies)));
            for ialpha = 1:Nalpha
                for izeta = zeta_to_impose_DKE
                    for ix = 1:Nx
                        indices = sfincs_indices(iSpecies, ix, 1:Nxi, ialpha, izeta, BLOCK_F, indexVars);
                        addSparseBlock(indices, indices, spatial_part(ialpha,izeta)*nuD(ix,iSpecies)*pitch_angle_scattering_operator_to_use)
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
% If there is a grid point at x=0, impose the appropriate boundary
% condition there.
% --------------------------------------------------

if whichMatrix ~= 2 && pointAtX0
    error('This version is not designed for pointAtX0')
    L = 0;
    if (whichMatrix==0 && L >= preconditioner_x_min_L)
        ddxToUse = ddx_preconditioner;
    else
        ddxToUse = ddx;
    end
    
    for ispecies = 1:Nspecies
        for ialpha = 1:Nalpha
            for izeta = zeta_to_impose_DKE
                % For L=0, force df/dx=0 at x=0 (regularity)
                L = 0;
                rowIndex = sfincs_indices(ispecies, 1, L+1, ialpha, izeta, BLOCK_F, indexVars);
                colIndices = sfincs_indices(ispecies, 1:Nx, L+1, ialpha, izeta, BLOCK_F, indexVars);
                addSparseBlock(rowIndex, colIndices, ddxToUse(1,:))
                
                % For L>0, set f=0 at x=0:
                indices = sfincs_indices(ispecies, 1, 2:Nxi, ialpha, izeta, BLOCK_F, indexVars);
                addToSparse(indices, indices, ones(size(indices)))
            end
        end
    end
end

% --------------------------------------------------
% Add density and pressure constraints.
% --------------------------------------------------

if whichMatrix ~= 2
    switch constraintScheme
        case 0
            % Do nothing.
            
        case 1
            for ialpha = 1:Nalpha
                for izeta = 1:Nzeta
                    factor = alphaWeights(ialpha)*zetaWeights(izeta)/(DHat(ialpha,izeta)*VPrimeHat);
                    for ispecies = 1:Nspecies
                        for ixi = 1:Nxi
                            colIndices = sfincs_indices(ispecies, 1:Nx, ixi, ialpha, izeta, BLOCK_F, indexVars);
                        
                            rowIndex = sfincs_indices(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT, indexVars);
                            addSparseBlock(rowIndex, colIndices, factor*xiWeights(ixi)*(x2.*xWeights)' )
                        
                            rowIndex = sfincs_indices(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT, indexVars);
                            addSparseBlock(rowIndex, colIndices, factor*xiWeights(ixi)*(x2.*x2.*xWeights)' )
                        end
                    end
                end
            end
            
        case 2
            for ialpha = 1:Nalpha
                for izeta = 1:Nzeta
                    factor = alphaWeights(ialpha)*zetaWeights(izeta)/(DHat(ialpha,izeta)*VPrimeHat);
                    for ispecies = 1:Nspecies
                        for ix = 1:Nx
                            rowIndex = sfincs_indices(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT, indexVars);
                            colIndices = sfincs_indices(ispecies, ix, 1:Nxi, ialpha, izeta, BLOCK_F, indexVars);
                            addSparseBlock(rowIndex, colIndices, factor * xiWeights')
                        end
                    end
                end
            end
            
        otherwise
            error('Invalid constraintScheme')
    end
end

% --------------------------------------------------
% Add sources.
% --------------------------------------------------

if whichMatrix ~= 2
    switch constraintScheme
        case 0
            % Do nothing
            
        case 1
            xPartOfSource1 = (1/(pi*sqrtpi))*(   -x2 + 5/2).*expx2;
            xPartOfSource2 = (1/(pi*sqrtpi))*(2/3*x2 -   1).*expx2;
            %xPartOfSource1 = (x2-5/2).*expx2;
            %xPartOfSource2 = (x2-3/2).*expx2;
            
            spatial_part_of_source = (B0OverBBar*BHat) ./ (GHat * DHat);
            
            for ialpha = 1:Nalpha
                for izeta = zeta_to_impose_DKE
                    stuff_to_add = spatial_part_of_source(ialpha,izeta)*ones(Nxi,1);
                    for ispecies = 1:Nspecies
                        for ix = ixMin:Nx
                            rowIndices = sfincs_indices(ispecies, ix, 1:Nxi, ialpha, izeta, BLOCK_F, indexVars);
                        
                            colIndex = sfincs_indices(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT, indexVars);
                            addSparseBlock(rowIndices, colIndex, xPartOfSource1(ix)*stuff_to_add)
                        
                            colIndex = sfincs_indices(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT, indexVars);
                            addSparseBlock(rowIndices, colIndex, xPartOfSource2(ix)*stuff_to_add)
                        end
                    end
                end
            end
            
        case 2
            spatial_part_of_source = (B0OverBBar*BHat) ./ (GHat * DHat);
            for ialpha = 1:Nalpha
                for izeta = zeta_to_impose_DKE
                    stuff_to_add = spatial_part_of_source(ialpha,izeta)*ones(Nxi,1);
                    for ispecies = 1:Nspecies
                        for ix = ixMin:Nx
                            colIndex = sfincs_indices(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT, indexVars);
                            rowIndices = sfincs_indices(ispecies, ix, 1:Nxi, ialpha, izeta, BLOCK_F, indexVars);
                            addSparseBlock(rowIndices, colIndex, stuff_to_add)
                        end
                    end
                end
            end
        otherwise
            error('Invalid constraintScheme')
    end
end

% --------------------------------------------------
% Add quasineutrality equation.
% --------------------------------------------------

% Part with column indices in BLOCK_F, which is also used for the residual:

if whichMatrix ~= 2 && includePhi1
    switch quasineutralityOption
        case 1
            ispecies_max = Nspecies;
        case 2
            ispecies_max = 1;
        otherwise
            error('Invalid quasineutralityOption')
    end
    L = 0;
    xPart = x2.*xWeights;
    speciesFactor = 4*pi*Zs .* ((THats./mHats).^(3/2));
    for ialpha = 1:Nalpha
        for izeta = 1:Nzeta
            rowIndex = sfincs_indices(1, 1, 1, ialpha, izeta, BLOCK_QN, indexVars);
            for ispecies = 1:ispecies_max
                colIndices = sfincs_indices(ispecies, 1:Nx, L+1, ialpha, izeta, BLOCK_F, indexVars);
                addSparseBlock(rowIndex, colIndices, xPart' *speciesFactor(ispecies))
            end
        end
    end
end

% For quasineutralityOption=1: Part of the Jacobian with column indices in BLOCK_QN, which is NOT used for the residual:
if includePhi1 && (quasineutralityOption==1) && (whichMatrix==0 || whichMatrix==1)
    stuffToAdd = zeros(Nalpha,Nzeta);
    for ispecies = 1:Nspecies
        stuffToAdd = stuffToAdd - gamma*Zs(ispecies)*Zs(ispecies)*nHats(ispecies)/THats(ispecies) ...
            *exp(-gamma*Zs(ispecies)/THats(ispecies)*Phi1Hat);
    end
    if withAdiabatic
        stuffToAdd = stuffToAdd - gamma*adiabaticZ*adiabaticZ*adiabaticNHat/adiabaticTHat ...
            *exp(-gamma*adiabaticZ/adiabaticTHat*Phi1Hat);
    end
    
    for ialpha = 1:Nalpha
        indices = sfincs_indices(1, 1, 1, ialpha, 1:Nzeta, BLOCK_QN, indexVars);
        addToSparse(indices, indices, stuffToAdd(ialpha,:))
    end
end

% For quasineutralityOption=2: Part of the Jacobian with column indices in BLOCK_QN, which IS used for the residual:
if includePhi1 && (quasineutralityOption==2) && (whichMatrix~=2)
    factor = -gamma*(Zs(1)*Zs(1)*nHats(1)/THats(1) ...
        +adiabaticZ*adiabaticZ*adiabaticNHat/adiabaticTHat);
    stuffToAdd = factor * ones(Nzeta,1);
    for ialpha = 1:Nalpha
        indices = sfincs_indices(1, 1, 1, ialpha, 1:Nzeta, BLOCK_QN, indexVars);
        addToSparse(indices, indices, stuffToAdd)
    end
end

% --------------------------------------------------
% Add Lagrange multiplier lambda
% --------------------------------------------------

if whichMatrix ~= 2 && includePhi1
    colIndex = sfincs_indices(1, 1, 1, 1, 1, BLOCK_PHI1_CONSTRAINT, indexVars);
    for ialpha = 1:Nalpha
        rowIndices = sfincs_indices(1, 1, 1, ialpha, 1:Nzeta, BLOCK_QN, indexVars);
        addSparseBlock(rowIndices, colIndex, ones(Nzeta,1))
    end
end

% --------------------------------------------------
% Add phi1 constraint.
% --------------------------------------------------

if whichMatrix ~= 2 && includePhi1
    for ialpha = 1:Nalpha
        colIndices = sfincs_indices(1, 1, 1, ialpha, 1:Nzeta, BLOCK_QN, indexVars);
        rowIndex = sfincs_indices(1, 1, 1, 1, 1, BLOCK_PHI1_CONSTRAINT, indexVars);
        addSparseBlock(rowIndex, colIndices, alphaWeights(ialpha)*(zetaWeights') ./ (DHat(ialpha,:)))
    end
end

% --------------------------------------------------
% End of adding entries to the matrix.
% --------------------------------------------------

fprintf('Time to contruct %s: %g seconds.\n',whichMatrixName,toc(populateMatrixTic))
tic
matrix = createSparse();
fprintf('Time to sparsify %s: %g seconds.\n',whichMatrixName,toc)
fprintf('This matrix has %d nonzeros. Fill fraction = %g. Original estimated nnz = %d\n',nnz(matrix), nnz(matrix)/(matrixSize*matrixSize),estimated_nnz)

index_lookup = zeros(matrixSize,5);
for ispecies = 1:Nspecies
    for ix = 1:Nx
        for ixi = 1:Nxi
            for ialpha = 1:Nalpha
                for izeta = 1:Nzeta
                    index = sfincs_indices(ispecies, ix, ixi, ialpha, izeta, BLOCK_F, indexVars);
                    index_lookup(index,1) = ispecies;
                    index_lookup(index,2) = ix;
                    index_lookup(index,3) = ixi;
                    index_lookup(index,4) = ialpha;
                    index_lookup(index,5) = izeta;
                end
            end
        end
    end
end
assignin('base','index_lookup',index_lookup)


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
        %fprintf('estimated nnz: %d   Actual value required: %d\n',estimated_nnz_original, sparseCreatorIndex)
        sparseMatrix = sparse(sparseCreator_i(1:(sparseCreatorIndex-1)), sparseCreator_j(1:(sparseCreatorIndex-1)), sparseCreator_s(1:(sparseCreatorIndex-1)), matrixSize, matrixSize);
        resetSparseCreator()
    end


end
