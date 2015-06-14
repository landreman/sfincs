function matrix = sfincs_populateMatrix(whichMatrix)

global THats mHats Zs nHats
global matrixSize Ntheta Nzeta Nxi Nx NL Nspecies pointAtX0
global ddtheta ddtheta_preconditioner ddzeta ddzeta_preconditioner
global ddx d2dx2 ddx_preconditioner
global xPotentials ddxPotentials d2dx2Potentials
global preconditioner_theta_min_L
global preconditioner_zeta_min_L
global preconditioner_species preconditioner_x preconditioner_x_min_L
global preconditioner_xi
global BLOCK_F BLOCK_QN BLOCK_PHI1_CONSTRAINT BLOCK_DENSITY_CONSTRAINT BLOCK_PRESSURE_CONSTRAINT BLOCK_F_CONSTRAINT
global includePhi1 nonlinear useDKESExBDrift includeXDotTerm includeElectricFieldTermInXiDot magneticDriftScheme
global Delta alpha nu_n dPhi1Hatdtheta dPhi1Hatdzeta

populateMatrixTic = tic;

if pointAtX0
    ixMin = 2;
else
    ixMin = 1;
end

x2 = x.*x;
expx2 = exp(-x2);

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
estimated_nnz = 1 * (Nx*Nx*Nspecies*Nspecies*Nxi*Ntheta*Nzeta ...
    + Nspecies*(nnz(ddtheta)*Nx*(3*Nxi)*Nzeta + nnz(ddzeta)*Nx*(3*Nxi)*Ntheta ...
    + Nx*(5*Nxi)*Ntheta*Nzeta + 3*Nx*Nx*Nxi*Ntheta*Nzeta ...
    + 2*2*Nx*1*Ntheta*Nzeta));

sparseCreatorIndex=1;
sparseCreator_i=0;
sparseCreator_j=0;
sparseCreator_s=0;
resetSparseCreator()

% -----------------------------------------
% Add collisionless terms:
% -----------------------------------------

for ispecies = 1:Nspecies
    THat = THats(ispecies);
    nHat = nHats(ispecies);
    mHat = mHats(ispecies);
    Z = Zs(ispecies);
    sqrtTHat = sqrt(THat);
    sqrtMass = sqrt(mHat);
    
    % -----------------------------------------
    % Add d/dtheta terms:
    % -----------------------------------------
    
    if whichMatrix ~= 2
        streamingTermSpatialPart = sqrtT/sqrtm*BHat_sup_theta./BHat;
        
        ExBTermSpatialPart = (alpha/2)*Delta*dPhiHatdpsiHat*DHat.*BHat_sub_zeta;
        if useDKESExBDrift
            ExBTermSpatialPart = ExBTermSpatialPart / FSABHat2;
        else
            ExBTermSpatialPart = ExBTermSpatialPart ./ (BHat.*BHat);
        end
        
        magneticDriftFactor = Delta*THat*DHat./(2*Z*(BHat.^3));
        if magneticDriftScheme>0
            magneticDriftSpatialPart1 = magneticDriftFactor .* (BHat_sub_zeta.*dBHatdpsiHat - BHat_sub_psi.*dBHatdzeta);
            magneticDriftSpatialPart2 = 2*BHat.*magneticDriftFactor .* (BHat_sub_psi_dzeta - BHat_sub_zeta_dpsiHat);
        else
            magneticDriftSpatialPart1 = zeros(Ntheta,Nzeta);
            magneticDriftSpatialPart2 = zeros(Ntheta,Nzeta);
        end
        if magneticDriftScheme==2
            magneticDriftSpatialPart1 = magneticDriftFactor .* BDotCurlB .* BHat_sup_theta ./ (BHat.*DHat);
        else
            magneticDriftSpatialPart3 = zeros(Ntheta,Nzeta);
        end
        
        if nonlinear
            nonlinearTermSpatialPart = -alpha*Delta*DHat.*BHat_sub_psi.*dPhi1Hatdzeta./(2*BHat.*BHat);
        else
            nonlinearTermSpatialPart = zeros(Ntheta,Nzeta);
        end
        
        for L=0:(Nxi-1)
            if whichMatrix==0 && L >= preconditioner_theta_min_L
                ddthetaToUse = ddtheta_preconditioner;
            else
                ddthetaToUse = ddtheta;
            end
                        
            for izeta = 1:Nzeta
                streamingTerm = diag(streamingTermSpatialPart(:,izeta))*ddthetaToUse;
                ExBTerm = diag(ExBTermSpatialPart(:,izeta))*ddthetaToUse;
                magneticDriftTerm1 = diag(magneticDriftSpatialPart1(:,izeta))*ddthetaToUse;
                magneticDriftTerm2 = diag(magneticDriftSpatialPart2(:,izeta))*ddthetaToUse;
                magneticDriftTerm3 = diag(magneticDriftSpatialPart3(:,izeta))*ddthetaToUse;
                nonlinearTerm = diag(nonlinearTermSpatialPart(:,izeta))*ddthetaToUse;
                
                for ix = ixMin:Nx
                    rowIndices = sfincs_indices(ispecies, ix, L+1, 1:Ntheta, izeta, BLOCK_F);

                    % Diagonal in L
                    colIndices = rowIndices;
                    addSparseBlock(rowIndices, colIndices, ExBTerm + nonlinearTerm ...
                        + x2(ix)*(magneticDriftTerm1+magneticDriftTerm2)*(L-1)*L/((2*L-3)*(2*L-1))...
                        + x2(ix)*magneticDriftTerm3*(-2)*L*(L+1)/((2*L+3)*(2*L-1)))
    
                    % Super-diagonal in L
                    ell = L + 1;
                    if (ell <= Nxi-1)
                        colIndices = sfincs_indices(ispecies, ix, ell+1, 1:Ntheta, izeta, BLOCK_F);
                        addSparseBlock(rowIndices, colIndices, streamingTerm*x(ix)*(L+1)/(2*L+3))
                    end
                    
                    % Sub-diagonal in L
                    ell = L - 1;
                    if (ell >= 0)
                        colIndices = sfincs_indices(ispecies, ix, ell+1, 1:Ntheta, izeta, BLOCK_F);
                        addSparseBlock(rowIndices, colIndices, streamingTerm*x(ix)*L/(2*L-1))
                    end
                    
                    if whichMatrix>0 || (preconditioner_xi==0)
                        % Super-super-diagonal in L
                        ell = L + 2;
                        if (ell <= Nxi-1)
                            colIndices = sfincs_indices(ispecies, ix, ell+1, 1:Ntheta, izeta, BLOCK_F);
                            addSparseBlock(rowIndices, colIndices, ...
                                x2(ix)*(magneticDriftTerm1+magneticDriftTerm2)*(L+2)*(L+1)/((2*L+5)*(2*L+3))...
                                + x2(ix)*magneticDriftTerm3*(-3)*(L+2)*(L+1)/((2*L+5)*(2*L+3)))
                        end
                        
                        % Sub-sub-diagonal in L
                        ell = L - 2;
                        if (ell >= 0)
                            colIndices = sfincs_indices(ispecies, ix, ell+1, 1:Ntheta, izeta, BLOCK_F);
                            addSparseBlock(rowIndices, colIndices, ...
                                x2(ix)*(magneticDriftTerm1+magneticDriftTerm2)*(L-1)*L/((2*L-3)*(2*L-1))...
                                + x2(ix)*magneticDriftTerm3*(-3)*(L-1)*L/((2*L-3)*(2*L-1)))
                        end
                    end
                end
            end
        end
    end
    
    % -----------------------------------------
    % Add d/dzeta terms:
    % -----------------------------------------
    
    if whichMatrix ~= 2
        streamingTermSpatialPart = sqrtT/sqrtm*BHat_sup_zeta./BHat;
        
        ExBTermSpatialPart = -(alpha/2)*Delta*dPhiHatdpsiHat*DHat.*BHat_sub_theta;
        if useDKESExBDrift
            ExBTermSpatialPart = ExBTermSpatialPart / FSABHat2;
        else
            ExBTermSpatialPart = ExBTermSpatialPart ./ (BHat.*BHat);
        end
        
        magneticDriftFactor = Delta*THat*DHat./(2*Z*(BHat.^3));
        if magneticDriftScheme>0
            magneticDriftSpatialPart1 = magneticDriftFactor .* (-BHat_sub_theta.*dBHatdpsiHat + BHat_sub_psi.*dBHatdtheta);
            magneticDriftSpatialPart2 = 2*BHat.*magneticDriftFactor .* (-BHat_sub_psi_dtheta + BHat_sub_theta_dpsiHat);
        else
            magneticDriftSpatialPart1 = zeros(Ntheta,Nzeta);
            magneticDriftSpatialPart2 = zeros(Ntheta,Nzeta);
        end
        if magneticDriftScheme==2
            magneticDriftSpatialPart1 = magneticDriftFactor .* BDotCurlB .* BHat_sup_zeta ./ (BHat.*DHat);
        else
            magneticDriftSpatialPart3 = zeros(Ntheta,Nzeta);
        end
        
        if nonlinear
            nonlinearTermSpatialPart = alpha*Delta*DHat.*BHat_sub_psi.*dPhi1Hatdtheta./(2*BHat.*BHat);
        else
            nonlinearTermSpatialPart = zeros(Ntheta,Nzeta);
        end
        
        for L=0:(Nxi-1)
            if whichMatrix==0 && L >= preconditioner_zeta_min_L
                ddzetaToUse = ddzeta_preconditioner;
            else
                ddzetaToUse = ddzeta;
            end
                        
            for itheta = 1:Ntheta
                streamingTerm = diag(streamingTermSpatialPart(itheta,:))*ddzetaToUse;
                ExBTerm = diag(ExBTermSpatialPart(itheta,:))*ddzetaToUse;
                magneticDriftTerm1 = diag(magneticDriftSpatialPart1(itheta,:))*ddzetaToUse;
                magneticDriftTerm2 = diag(magneticDriftSpatialPart2(itheta,:))*ddzetaToUse;
                magneticDriftTerm3 = diag(magneticDriftSpatialPart3(itheta,:))*ddzetaToUse;
                nonlinearTerm = diag(nonlinearTermSpatialPart(itheta,:))*ddzetaToUse;
                
                for ix = ixMin:Nx
                    rowIndices = sfincs_indices(ispecies, ix, L+1, itheta, 1:Nzeta, BLOCK_F);

                    % Diagonal in L
                    colIndices = rowIndices;
                    addSparseBlock(rowIndices, colIndices, ExBTerm + nonlinearTerm ...
                        + x2(ix)*(magneticDriftTerm1+magneticDriftTerm2)*(L-1)*L/((2*L-3)*(2*L-1))...
                        + x2(ix)*magneticDriftTerm3*(-2)*L*(L+1)/((2*L+3)*(2*L-1)))
    
                    % Super-diagonal in L
                    ell = L + 1;
                    if (ell <= Nxi-1)
                        colIndices = sfincs_indices(ispecies, ix, ell+1, itheta, 1:Nzeta, BLOCK_F);
                        addSparseBlock(rowIndices, colIndices, streamingTerm*x(ix)*(L+1)/(2*L+3))
                    end
                    
                    % Sub-diagonal in L
                    ell = L - 1;
                    if (ell >= 0)
                        colIndices = sfincs_indices(ispecies, ix, ell+1, itheta, 1:Nzeta, BLOCK_F);
                        addSparseBlock(rowIndices, colIndices, streamingTerm*x(ix)*L/(2*L-1))
                    end
                    
                    if whichMatrix>0 || (preconditioner_xi==0)
                        % Super-super-diagonal in L
                        ell = L + 2;
                        if (ell <= Nxi-1)
                            colIndices = sfincs_indices(ispecies, ix, ell+1, itheta, 1:Nzeta, BLOCK_F);
                            addSparseBlock(rowIndices, colIndices, ...
                                x2(ix)*(magneticDriftTerm1+magneticDriftTerm2)*(L+2)*(L+1)/((2*L+5)*(2*L+3))...
                                + x2(ix)*magneticDriftTerm3*(-3)*(L+2)*(L+1)/((2*L+5)*(2*L+3)))
                        end
                        
                        % Sub-sub-diagonal in L
                        ell = L - 2;
                        if (ell >= 0)
                            colIndices = sfincs_indices(ispecies, ix, ell+1, itheta, 1:Nzeta, BLOCK_F);
                            addSparseBlock(rowIndices, colIndices, ...
                                x2(ix)*(magneticDriftTerm1+magneticDriftTerm2)*(L-1)*L/((2*L-3)*(2*L-1))...
                                + x2(ix)*magneticDriftTerm3*(-3)*(L-1)*L/((2*L-3)*(2*L-1)))
                        end
                    end
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
    % Add the collisionless df1/dx term:
    % -----------------------------------------
    
    if whichMatrix ~= 2
        if includeXDotTerm
            factor = -alpha*Delta*DHat./(4*(BHat.^3));
            ErTermSpatialPart1 = factor .* (BHat_sub_theta.*dBHatdzeta - BHat_sub_zeta.*dBHatdtheta);
            if ~force0RadialCurrentInEquilibrium
                ErTermSpatialPart2 = factor.*(2*BHat).*(dBHat_sub_zeta_dtheta - dBHat_sub_theta_dzeta);
            else
                ErTermSpatialPart2 = zeros(Ntheta,Nzeta);
            end
        else
            ErTermSpatialPart1 = zeros(Ntheta,Nzeta);
            ErTermSpatialPart2 = zeros(Ntheta,Nzeta);
        end
        
        if nonlinear
            nonlinearTermSpatialPart = -alpha*Z./(2*sqrtT*sqrtm*BHat).*(BHat_sup_theta.*dPhi1Hatdtheta + BHat_sub_zeta.*dPhi1Hatdzeta);
        else
            nonlinearTermSpatialPart = zeros(Ntheta,Nzeta);
        end
        
        for L = 0:(Nxi-1)
            if whichMatrix==0 && L >= preconditioner_x_min_L
                ddxToUse = ddx_preconditioner;
            else
                ddxToUse = ddx;
            end
            if pointAtX0
                % Do not enforce the kinetic equation at x=0:
                ddxToUse(1,:)=0;
            end
            xddxToUse = diag(x)*ddxToUse;
            
            for itheta = 1:Ntheta
                for izeta = 1:Nzeta
                    rowIndices = sfincs_indices(ispecies, 1:Nx, L+1, itheta, izeta, BLOCK_F);

                    % Diagonal in L
                    colIndices = rowIndices;
                    addSparseBlock(rowIndices, colIndices, ...
                        xddxToUse*ErTermSpatialPart1(itheta,izeta)*2*(3*L*L+3*L-2)/((2*L+3)*(2*L-1))...
                        + ErTermSpatialPart2(itheta,izeta)*(2*L*L+2*L-1)/((2*L+3)*(2*L-1)))
                    
                    % Super-diagonal in L
                    ell = L + 1;
                    if (ell <= Nxi-1)
                        colIndices = sfincs_indices(ispecies, 1:Nx, ell+1, itheta, izeta, BLOCK_F);
                        addSparseBlock(rowIndices, colIndices, ddxToUse*nonlinearTermSpatialPart(itheta,izeta)*(L+1)/(2*L+3))
                    end
                    
                    % Sub-diagonal in L
                    ell = L - 1;
                    if (ell >= 0)
                        colIndices = sfincs_indices(ispecies, 1:Nx, ell+1, itheta, izeta, BLOCK_F);
                        addSparseBlock(rowIndices, colIndices, ddxToUse*nonlinearTermSpatialPart(itheta,izeta)*L/(2*L-1))
                    end

                    if whichMatrix>0 || (preconditioner_xi==0)
                        % Super-super-diagonal in L
                        ell = L + 2;
                        if (ell <= Nxi-1)
                            colIndices = sfincs_indices(ispecies, 1:Nx, ell+1, itheta, izeta, BLOCK_F);
                            addSparseBlock(rowIndices, colIndices, ...
                                xddxToUse*(ErTermSpatialPart1(itheta,izeta)+ErTermSpatialPart2(itheta,izeta))...
                                *(L+2)*(L+1)/((2*L+5)*(2*L+3)))
                        end
                        
                        % Sub-sub-diagonal in L
                        ell = L - 2;
                        if (ell >= 0)
                            colIndices = sfincs_indices(ispecies, 1:Nx, ell+1, itheta, izeta, BLOCK_F);
                            addSparseBlock(rowIndices, colIndices, ...
                                xddxToUse*(ErTermSpatialPart1(itheta,izeta)+ErTermSpatialPart2(itheta,izeta))...
                                *(L-1)*L/((2*L-3)*(2*L-1)))
                        end
                    end
                end
            end
        end
    end
    
    
    
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
    
    
    
end

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
                                error('Invalid preconditionerMethod_x')
                        end
                        
                    end
                    
                    % At this point, CHat holds the collision operator
                    % divided by \bar{nu}
                    
                    for itheta = 1:Ntheta
                        for izeta = 1:Nzeta
                            rowIndices = getIndices(speciesA, 1:Nx, L+1, itheta, izeta, 0);
                            colIndices = getIndices(speciesB, 1:Nx, L+1, itheta, izeta, 0);
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
                        addToSparse(indices, indices, -nu_n*CHat)
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
    L = 0;
    if (whichMatrix==0 && L >= preconditioner_x_min_L)
        ddxToUse = ddx_preconditioner;
    else
        ddxToUse = ddx;
    end
    
    for ispecies = 1:Nspecies
        for itheta = 1:Ntheta
            for izeta = 1:Nzeta
                % For L=0, force df/dx=0 at x=0 (regularity)
                L = 0;
                rowIndex = sfincs_indices(ispecies, 1, L+1, itheta, izeta, BLOCK_F);
                colIndices = sfincs_indices(ispecies, 1:Nx, L+1, itheta, izeta, BLOCK_F);
                addSparseBlock(rowIndex, colIndices, ddxToUse(1,:))
                
                % For L>0, set f=0 at x=0:
                indices = sfincs_indices(ispecies, 1, 2:Nxi, itheta, izeta, BLOCK_F);
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
            L=0;
            for ispecies = 1:Nspecies
                for itheta = 1:Ntheta
                    for izeta = 1:Nzeta
                        colIndices = sfincs_indices(ispecies, 1:Nx, L+1, itheta, izeta, BLOCK_F);
                        
                        rowIndex = sfincs_indices(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT);
                        addSparseBlock(rowIndex, colIndices, thetaWeights(itheta)*zetaWeights(izeta)*(x2.*xWeights)' / (DHat(itheta,izeta)))
                        
                        rowIndex = sfincs_indices(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT);
                        addSparseBlock(rowIndex, colIndices, thetaWeights(itheta)*zetaWeights(izeta)*(x2.*x2.*xWeights)' / (DHat(itheta,izeta)))
                    end
                end
            end
            
        case 2
            L=0;
            for ispecies = 1:Nspecies
                % I think this loop should go from 1 rather than from
                % ixMin. But I could be convinced otherwise.
                for ix = 1:Nx
                    rowIndex = sfincs_indices(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT);
                    for itheta = 1:Ntheta
                        colIndices = sfincs_indices(ispecies, ix, L+1, itheta, 1:Nzeta, BLOCK_F);
                        addSparseBlock(rowIndex, colIndices, thetaWeights(itheta)*zetaWeights./DHat(itheta,:))
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
            xPartOfSource1 = (x2-5/2).*expx2;
            xPartOfSource2 = (x2-3/2).*expx2;
            
            L=0;
            for ispecies = 1:Nspecies
                for ix = ixMin:Nx
                    for itheta = 1:Ntheta
                        rowIndices = sfincs_indices(ispecies, ix, L+1, itheta, 1:Nzeta, BLOCK_F);
                        
                        colIndex = sfincs_indices(ispecies, 1, 1, 1, 1, BLOCK_DENSITY_CONSTRAINT);
                        addSparseBlock(rowIndices, colIndex, xPartOfSource1(ix)*ones(Nzeta,1))
                        
                        colIndex = sfincs_indices(ispecies, 1, 1, 1, 1, BLOCK_PRESSURE_CONSTRAINT);
                        addSparseBlock(rowIndices, colIndex, xPartOfSource2(ix)*ones(Nzeta,1))
                    end
                end
            end
            
        case 2
            L=0;
            for ispecies = 1:Nspecies
                for ix = ixMin:Nx
                    colIndex = sfincs_indices(ispecies, ix, 1, 1, 1, BLOCK_F_CONSTRAINT);
                    for itheta = 1:Ntheta
                        rowIndices = sfincs_indices(ispecies, ix, L+1, itheta, 1:Nzeta, BLOCK_F);
                        addSparseBlock(rowIndices, colIndex, ones(Nzeta,1))
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

if whichMatrix ~= 2 && includePhi1
    L = 0;
    xPart = x2.*xWeights;
    speciesFactor = Zs .* ((THats./mHats).^(3/2));
    for itheta = 1:Ntheta
        for izeta = 1:Nzeta
            rowIndex = sfincs_indices(1, 1, 1, itheta, izeta, BLOCK_QN);
            for ispecies = 1:Nspecies
                colIndices = sfincs_indices(ispecies, 1:Nx, L+1, itheta, izeta, BLOCK_F);
                addSparseBlock(rowIndex, colIndices, xPart*speciesFactor(ispecies))
            end
        end
    end
end

% --------------------------------------------------
% Add Lagrange multiplier lambda
% --------------------------------------------------

if whichMatrix ~= 2 && includePhi1
    colIndex = sfincs_indices(1, 1, 1, 1, 1, BLOCK_PHI1_CONSTRAINT);
    for itheta = 1:Ntheta
        rowIndices = sfincs_indices(1, 1, 1, itheta, 1:Nzeta, BLOCK_QN);
        addSparseBlock(rowIndices, colIndex, ones(Nzeta,1))
    end
end

% --------------------------------------------------
% Add phi1 constraint.
% --------------------------------------------------

if whichMatrix ~= 2 && includePhi1
    for itheta = 1:Ntheta
        colIndices = sfincs_indices(1, 1, 1, itheta, 1:Nzeta, BLOCK_QN);
        rowIndex = sfincs_indices(1, 1, 1, 1, 1, BLOCK_PHI1_CONSTRAINT);
        addSparseBlock(rowIndex, colIndices, thetaWeights(itheta)*zetaWeights ./ (DHat(itheta,:)))
    end
end

% --------------------------------------------------
% End of adding entries to the matrix.
% --------------------------------------------------

fprintf('Time to contruct %s: %g seconds.\n',whichMatrixName,toc(populateMatrixTic))
tic
matrix = createSparse();
fprintf('Time to sparsify %s: %g seconds.\n',whichMatrixName,toc)
fprintf('This matrix has %d nonzeros. Fill fraction = %g\n',nnz(matrix), nnz(matrix)/(matrixSize*matrixSize))


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


end