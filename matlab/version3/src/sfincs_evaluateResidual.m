function residual = sfincs_evaluateResidual()

global stateVector matrixSize RHSMode pointAtX0 dPhiHatdpsiHat includeTemperatureEquilibrationTerm f0
global x Nspecies Ntheta Nzeta Nx Delta alpha BLOCK_F indexVars
global Zs THats mHats nHats dnHatdpsiHats dTHatdpsiHats EParallelHat
global BHat DHat FSABHat2 BHat_sub_zeta BHat_sub_theta
global dBHatdtheta dBHatdzeta

fprintf('Evaluating residual.\n')

if norm(stateVector)>1e-100
    % Part of the residual comes from multiplying the state vector by a big
    % matrix.
    whichMatrix = 3;
    residualMatrix = sfincs_populateMatrix(whichMatrix);
    residual = residualMatrix * stateVector;
else
    % There is no need to assemble that matrix since the state vector is 0.
    residual = zeros(matrixSize,1);
    fprintf('The state vector is 0 so I will skip building the first matrix when evaluating the residual.\n')
end

if includeTemperatureEquilibrationTerm
    whichMatrix = 2;
    residualMatrix = sfincs_populateMatrix(whichMatrix);
    residual = residual + residualMatrix * f0;
end

% Next, evaluate the remaining inhomogeneous terms.

if RHSMode == 1
    dPhiHatdpsiHatToUseInRHS = dPhiHatdpsiHat;
else
    dPhiHatdpsiHatToUseInRHS = 0;
end

if pointAtX0
    ixMin = 2;
else
    ixMin = 1;
end

x2 = x.*x;
expx2 = exp(-x2);
sqrtpi = sqrt(pi);

allZeta = 1:Nzeta;
rhs = zeros(matrixSize,1);
for ispecies = 1:Nspecies
    Z = Zs(ispecies);
    THat = THats(ispecies);
    mHat = mHats(ispecies);
    nHat = nHats(ispecies);
    sqrtTHat = sqrt(THat);
    sqrtmHat = sqrt(mHat);
    
    spatialFactor = Delta*nHat*mHat*sqrtmHat ...
        ./(2*pi*sqrtpi*Z*(BHat.^3)*sqrtTHat) ...
        .*(BHat_sub_zeta.*dBHatdtheta - BHat_sub_theta.*dBHatdzeta) ...
        .* DHat;
    
    for ix = ixMin:Nx
        xPartOfRHS = x2(ix)*expx2(ix)*( dnHatdpsiHats(ispecies)/nHat ...
            + alpha*Z/THat*dPhiHatdpsiHatToUseInRHS ...
            + (x2(ix) - 3/2)*dTHatdpsiHats(ispecies)/THat);
        
        inductiveFactor = alpha*Z*x(ix)*expx2(ix)*EParallelHat ...
            *nHat*mHat/(pi*sqrtpi*THat*THat*FSABHat2);
    
        for itheta = 1:Ntheta
            % Gradient terms:
            
            L = 0;
            indices = sfincs_indices(ispecies, ix, L+1, itheta, allZeta, BLOCK_F, indexVars);
            rhs(indices) =  (4/3)*xPartOfRHS*spatialFactor(itheta,:);
            
            L = 2;
            indices = sfincs_indices(ispecies, ix, L+1, itheta, allZeta, BLOCK_F, indexVars);
            rhs(indices) = (2/3)*xPartOfRHS*spatialFactor(itheta,:);
            
            % Inductive term:
            L = 1;
            indices = sfincs_indices(ispecies, ix, L+1, itheta, allZeta, BLOCK_F, indexVars);
            rhs(indices) = inductiveFactor * BHat(itheta,:);
            
        end
    end
end

residual = residual - rhs;

end