function residual = sfincs_evaluateResidual()

global stateVector matrixSize RHSMode pointAtX0 dPhiHatdpsiHat includeTemperatureEquilibrationTerm f0
global x Nspecies Nalpha Nzeta Nx Delta gamma BLOCK_F BLOCK_QN indexVars
global Zs THats mHats nHats dnHatdpsiHats dTHatdpsiHats EParallelHat
global BHat DHat FSABHat2 BHat_sub_zeta BHat_sub_theta
global dBHatdtheta dBHatdzeta Phi1Hat includePhi1 includePhi1InKineticEquation
global quasineutralityOption withAdiabatic adiabaticZ adiabaticNHat adiabaticTHat

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

% Add R_m, the part of the residual involving the radial magnetic drift.
% Also add the inductive parallel electric field term. 
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
    
    if includePhi1 && includePhi1InKineticEquation
        spatialFactor = spatialFactor .* exp(-(gamma*Z/THat)*Phi1Hat);
    end
    
    for ix = ixMin:Nx
        if includePhi1 && includePhi1InKineticEquation
            xAndSpatialPartOfRHS = x2(ix)*expx2(ix)*( dnHatdpsiHats(ispecies)/nHat ...
                + gamma*Z/THat*dPhiHatdpsiHatToUseInRHS ...
                + (x2(ix) - 3/2 + (gamma*Z/THat)*Phi1Hat)*dTHatdpsiHats(ispecies)/THat) ...
                .* spatialFactor;
        else
            xAndSpatialPartOfRHS = x2(ix)*expx2(ix)*( dnHatdpsiHats(ispecies)/nHat ...
                + gamma*Z/THat*dPhiHatdpsiHatToUseInRHS ...
                + (x2(ix) - 3/2)*dTHatdpsiHats(ispecies)/THat) ...
                * spatialFactor;
        end
        
        inductiveFactor = gamma*Z*x(ix)*expx2(ix)*EParallelHat ...
            *nHat*mHat/(pi*sqrtpi*THat*THat*FSABHat2);
    
        for ialpha = 1:Nalpha
            % Gradient terms:
            
            L = 0;
            indices = sfincs_indices(ispecies, ix, L+1, ialpha, allZeta, BLOCK_F, indexVars);
            rhs(indices) =  (4/3)*xAndSpatialPartOfRHS(ialpha,:);
            
            L = 2;
            indices = sfincs_indices(ispecies, ix, L+1, ialpha, allZeta, BLOCK_F, indexVars);
            rhs(indices) = (2/3)*xAndSpatialPartOfRHS(ialpha,:);
            % Inductive term:
            L = 1;
            indices = sfincs_indices(ispecies, ix, L+1, ialpha, allZeta, BLOCK_F, indexVars);
            rhs(indices) = inductiveFactor * BHat(ialpha,:);
            
        end
    end
end

% Add Z n_0 exp(-Z e Phi1/T) term in quasineutrality:
if includePhi1 && (quasineutralityOption==1)
    stuffToAdd = zeros(Nalpha,Nzeta);
    for ispecies = 1:Nspecies
        stuffToAdd = stuffToAdd + Zs(ispecies)*nHats(ispecies) ...
            *exp(-gamma*Zs(ispecies)/THats(ispecies)*Phi1Hat);
    end
    if withAdiabatic
        stuffToAdd = stuffToAdd + adiabaticZ*adiabaticNHat ...
            *exp(-gamma*adiabaticZ/adiabaticTHat*Phi1Hat);
    end
    
    for ialpha = 1:Nalpha
        indices = sfincs_indices(1, 1, 1, ialpha, 1:Nzeta, BLOCK_QN, indexVars);
        rhs(indices) = -stuffToAdd(ialpha,:);
    end
end

residual = residual - rhs;

end
