function sfincs_diagnostics()

global stateVector f0 constraintScheme sources force0RadialCurrentInEquilibrium
global Ntheta Nzeta Nspecies Nx thetaWeights zetaWeights RHSMode
global Phi1Hat dPhi1Hatdtheta dPhi1Hatdzeta ddtheta ddzeta includePhi1
global BLOCK_F BLOCK_QN BLOCK_PHI1_CONSTRAINT BLOCK_DENSITY_CONSTRAINT BLOCK_PRESSURE_CONSTRAINT BLOCK_F_CONSTRAINT
global BHat DHat BHat_sub_theta BHat_sub_zeta dBHatdtheta dBHatdzeta dBHat_sub_zeta_dtheta dBHat_sub_theta_dzeta
global THats Zs mHats nHats x xWeights Delta alpha
global B0OverBBar FSABHat2 VPrimeHat

global FSADensityPerturbation FSABFlow FSAPressurePerturbation
global particleFlux_vm0_psiHat particleFlux_vm_psiHat particleFlux_vE0_psiHat particleFlux_vE_psiHat
global momentumFlux_vm0_psiHat momentumFlux_vm_psiHat momentumFlux_vE0_psiHat momentumFlux_vE_psiHat
global heatFlux_vm0_psiHat heatFlux_vm_psiHat heatFlux_vE0_psiHat heatFlux_vE_psiHat
global jHat FSABjHat FSABjHatOverB0 FSABjHatOverRootFSAB2
global totalDensity totalPressure velocityUsingFSADensity velocityUsingTotalDensity MachUsingFSAThermalSpeed

fprintf('Computing diagnostics.\n')

if includePhi1
    for itheta = 1:Ntheta
        indices = sfincs_indices(1,1,1,itheta,1:Nzeta,BLOCK_QN);
        Phi1Hat(itheta,:) = stateVector(indices);
    end
    
    dPhi1Hatdtheta = ddtheta * Phi1Hat;
    dPhi1Hatdzeta = (ddzeta * (Phi1Hat'))';
    
    index = sfincs_indices(1,1,1,1,1,BLOCK_PHI1_CONSTRAINT);
    lambda = stateVector(index);
else
    lambda = 0;
end

switch constraintScheme
    case 0
        sources = [];
    case 1
        sources = zeros(Nspecies,2);
        for ispecies = 1:Nspecies
            sources(ispecies,1) = stateVector(sfincs_indices(ispecies,1,1,1,1,BLOCK_DENSITY_CONSTRAINT));
            sources(ispecies,2) = stateVector(sfincs_indices(ispecies,1,1,1,1,BLOCK_PRESSURE_CONSTRAINT));
        end
    case 2
        sources = zeros(Nspecies,Nx);
        for ispecies = 1:Nspecies
            sources(ispecies,:) = stateVector(sfincs_indices(ispecies,1:Nx,1,1,1,BLOCK_F_CONSTRAINT));
        end
    otherwise
        error('Invalid constraintScheme')
end

fullf = stateVector + f0;

jHat = zeros(Ntheta,Nzeta);
densityPerturbation = zeros(Nspecies,Ntheta,Nzeta);
totalDensity = zeros(Nspecies,Ntheta,Nzeta);
flow = zeros(Nspecies,Ntheta,Nzeta);
velocityUsingFSADensity = zeros(Nspecies,Ntheta,Nzeta);
velocityUsingTotalDensity = zeros(Nspecies,Ntheta,Nzeta);
MachUsingFSAThermalSpeed = zeros(Nspecies,Ntheta,Nzeta);
pressurePerturbation = zeros(Nspecies,Ntheta,Nzeta);
totalPressure = zeros(Nspecies,Ntheta,Nzeta);
particleFluxBeforeSurfaceIntegral_vm0 = zeros(Nspecies,Ntheta,Nzeta);
particleFluxBeforeSurfaceIntegral_vm = zeros(Nspecies,Ntheta,Nzeta);
particleFluxBeforeSurfaceIntegral_vE0 = zeros(Nspecies,Ntheta,Nzeta);
particleFluxBeforeSurfaceIntegral_vE = zeros(Nspecies,Ntheta,Nzeta);
momentumFluxBeforeSurfaceIntegral_vm0 = zeros(Nspecies,Ntheta,Nzeta);
momentumFluxBeforeSurfaceIntegral_vm = zeros(Nspecies,Ntheta,Nzeta);
momentumFluxBeforeSurfaceIntegral_vE0 = zeros(Nspecies,Ntheta,Nzeta);
momentumFluxBeforeSurfaceIntegral_vE = zeros(Nspecies,Ntheta,Nzeta);
heatFluxBeforeSurfaceIntegral_vm0 = zeros(Nspecies,Ntheta,Nzeta);
heatFluxBeforeSurfaceIntegral_vm = zeros(Nspecies,Ntheta,Nzeta);
heatFluxBeforeSurfaceIntegral_vE0 = zeros(Nspecies,Ntheta,Nzeta);
heatFluxBeforeSurfaceIntegral_vE = zeros(Nspecies,Ntheta,Nzeta);

FSADensityPerturbation = zeros(Nspecies,1);
FSABFlow = zeros(Nspecies,1);
FSAPressurePerturbation = zeros(Nspecies,1);
particleFlux_vm0_psiHat = zeros(Nspecies,1);
particleFlux_vm_psiHat = zeros(Nspecies,1);
particleFlux_vE0_psiHat = zeros(Nspecies,1);
particleFlux_vE_psiHat = zeros(Nspecies,1);
momentumFlux_vm0_psiHat = zeros(Nspecies,1);
momentumFlux_vm_psiHat = zeros(Nspecies,1);
momentumFlux_vE0_psiHat = zeros(Nspecies,1);
momentumFlux_vE_psiHat = zeros(Nspecies,1);
heatFlux_vm0_psiHat = zeros(Nspecies,1);
heatFlux_vm_psiHat = zeros(Nspecies,1);
heatFlux_vE0_psiHat = zeros(Nspecies,1);
heatFlux_vE_psiHat = zeros(Nspecies,1);

factor = (BHat_sub_theta.*dBHatdzeta - BHat_sub_zeta.*dBHatdtheta) ./ (BHat.^3);
if force0RadialCurrentInEquilibrium
    factor2 = zeros(Ntheta,Nzeta);
else
    factor2 = 2 * (dBHat_sub_zeta_dtheta-dBHat_sub_theta_dzeta) ./ (BHat.*BHat);
end
factor_vE = (BHat_sub_theta.*dPhi1Hatdzeta - BHat_sub_zeta.*dPhi1Hatdtheta) ./ (BHat.*BHat);

for ispecies = 1:Nspecies
    THat = THats(ispecies);
    nHat = nHats(ispecies);
    mHat = mHats(ispecies);
    Z = Zs(ispecies);
    sqrtT = sqrt(THat);
    sqrtm = sqrt(mHat);
    
    densityWeights = 4*pi*THat*sqrtT/(mHat*sqrtm)*(x.*x.*xWeights)';
    flowWeights = 4*pi*THat*THat/(3*mHat*mHat)*(x.*x.*x.*xWeights)';
    pressureWeights = 8*pi*THat*THat*sqrtT/(3*mHat*sqrtm)*(x.*x.*x.*x.*xWeights)';
    particleFluxWeights_vm = pi*Delta*THat*THat*sqrtT/(Z*VPrimeHat*mHat*sqrtm)*(x.*x.*x.*x.*xWeights)';
    particleFluxWeights_vE = 2*alpha*pi*Delta*THat*sqrtT/(VPrimeHat*mHat*sqrtm)*(x.*x.*xWeights)';
    momentumFluxWeights_vm = pi*Delta*THat*THat*THat/(Z*VPrimeHat*mHat)*(x.*x.*x.*x.*x.*xWeights)';
    momentumFluxWeights_vE = 2*alpha*pi*Delta*THat*THat/(VPrimeHat*mHat)*(x.*x.*x.*xWeights)';
    heatFluxWeights_vm = pi*Delta*THat*THat*THat*sqrtT/(2*Z*VPrimeHat*mHat*sqrtm)*(x.*x.*x.*x.*x.*x.*xWeights)';
    heatFluxWeights_vE = 2*alpha*pi*Delta*THat*THat*sqrtT/(2*VPrimeHat*mHat*sqrtm)*(x.*x.*x.*x.*xWeights)';
        
    % Do all the integrals over x and xi:
    for itheta = 1:Ntheta
        for izeta = 1:Nzeta
            
            L = 0;
            indices = sfincs_indices(ispecies, 1:Nx, L+1, itheta, izeta, BLOCK_F);
            densityPerturbation(ispecies, itheta, izeta) = densityWeights * stateVector(indices);
            pressurePerturbation(ispecies, itheta, izeta) = pressureWeights * stateVector(indices);
            particleFluxBeforeSurfaceIntegral_vm0(ispecies, itheta, izeta) = (factor(itheta,izeta)*8/3+factor2(itheta,izeta)*2/3)*particleFluxWeights_vm * f0(indices);
            particleFluxBeforeSurfaceIntegral_vm(ispecies, itheta, izeta) = (factor(itheta,izeta)*8/3+factor2(itheta,izeta)*2/3)*particleFluxWeights_vm * fullf(indices);
            particleFluxBeforeSurfaceIntegral_vE0(ispecies, itheta, izeta) = factor_vE(itheta,izeta) * particleFluxWeights_vE * f0(indices);
            particleFluxBeforeSurfaceIntegral_vE(ispecies, itheta, izeta) = factor_vE(itheta,izeta) * particleFluxWeights_vE * fullf(indices);
            heatFluxBeforeSurfaceIntegral_vm0(ispecies, itheta, izeta) = (factor(itheta,izeta)*8/3+factor2(itheta,izeta)*2/3)*heatFluxWeights_vm * f0(indices);
            heatFluxBeforeSurfaceIntegral_vm(ispecies, itheta, izeta) = (factor(itheta,izeta)*8/3+factor2(itheta,izeta)*2/3)*heatFluxWeights_vm * fullf(indices);
            heatFluxBeforeSurfaceIntegral_vE0(ispecies, itheta, izeta) = factor_vE(itheta,izeta) * heatFluxWeights_vE * f0(indices);
            heatFluxBeforeSurfaceIntegral_vE(ispecies, itheta, izeta) = factor_vE(itheta,izeta) * heatFluxWeights_vE * fullf(indices);
            
            L = 1;
            indices = sfincs_indices(ispecies, 1:Nx, L+1, itheta, izeta, BLOCK_F);
            flow(ispecies, itheta, izeta) = flowWeights * stateVector(indices);
            momentumFluxBeforeSurfaceIntegral_vm0(ispecies, itheta, izeta) = (factor(itheta,izeta)*16/15+factor2(itheta,izeta)*2/5)*BHat(itheta,izeta)*momentumFluxWeights_vm * f0(indices);
            momentumFluxBeforeSurfaceIntegral_vm(ispecies, itheta, izeta) = (factor(itheta,izeta)*16/15+factor2(itheta,izeta)*2/5)*BHat(itheta,izeta)*momentumFluxWeights_vm * fullf(indices);
            momentumFluxBeforeSurfaceIntegral_vE0(ispecies, itheta, izeta) = factor_vE(itheta,izeta)*2/3*BHat(itheta,izeta)*momentumFluxWeights_vE * f0(indices);
            momentumFluxBeforeSurfaceIntegral_vE(ispecies, itheta, izeta) = factor_vE(itheta,izeta)*2/3*BHat(itheta,izeta)*momentumFluxWeights_vE * fullf(indices);
            
            L = 2;
            indices = sfincs_indices(ispecies, 1:Nx, L+1, itheta, izeta, BLOCK_F);
            particleFluxBeforeSurfaceIntegral_vm0(ispecies, itheta, izeta) = particleFluxBeforeSurfaceIntegral_vm0(ispecies, itheta, izeta) + (factor(itheta,izeta)+factor2(itheta,izeta))*4/15*particleFluxWeights_vm * f0(indices);
            particleFluxBeforeSurfaceIntegral_vm(ispecies, itheta, izeta) = particleFluxBeforeSurfaceIntegral_vm(ispecies, itheta, izeta) + (factor(itheta,izeta)+factor2(itheta,izeta))*4/15*particleFluxWeights_vm * fullf(indices);
            heatFluxBeforeSurfaceIntegral_vm0(ispecies, itheta, izeta) = heatFluxBeforeSurfaceIntegral_vm0(ispecies, itheta, izeta) + (factor(itheta,izeta)+factor2(itheta,izeta))*4/15*heatFluxWeights_vm * f0(indices);
            heatFluxBeforeSurfaceIntegral_vm(ispecies, itheta, izeta) = heatFluxBeforeSurfaceIntegral_vm(ispecies, itheta, izeta) + (factor(itheta,izeta)+factor2(itheta,izeta))*4/15*heatFluxWeights_vm * fullf(indices);
            
            L = 3;
            indices = sfincs_indices(ispecies, 1:Nx, L+1, itheta, izeta, BLOCK_F);
            momentumFluxBeforeSurfaceIntegral_vm0(ispecies, itheta, izeta) = momentumFluxBeforeSurfaceIntegral_vm0(ispecies, itheta, izeta) + (factor(itheta,izeta)+factor2(itheta,izeta))*4/35*BHat(itheta,izeta)*momentumFluxWeights_vm * f0(indices);
            momentumFluxBeforeSurfaceIntegral_vm(ispecies, itheta, izeta) = momentumFluxBeforeSurfaceIntegral_vm(ispecies, itheta, izeta) + (factor(itheta,izeta)+factor2(itheta,izeta))*4/35*BHat(itheta,izeta)*momentumFluxWeights_vm * fullf(indices);
        end
    end
    
    % Done with the integrals over x and xi.
    % Now do the integrals over theta and zeta.
    
    FSADensityPerturbation(ispecies) = (thetaWeights') * (squeeze(densityPerturbation(ispecies,:,:)) ./ DHat) * zetaWeights;
    FSABFlow(ispecies) = (thetaWeights') * (squeeze(flow(ispecies,:,:)) .* BHat./ DHat) * zetaWeights;
    FSAPressurePerturbation(ispecies) = (thetaWeights') * (squeeze(pressurePerturbation(ispecies,:,:)) ./ DHat) * zetaWeights;
    particleFlux_vm0_psiHat(ispecies) = (thetaWeights') * (squeeze(particleFluxBeforeSurfaceIntegral_vm0(ispecies,:,:))) * zetaWeights;
    particleFlux_vm_psiHat(ispecies) = (thetaWeights') * (squeeze(particleFluxBeforeSurfaceIntegral_vm(ispecies,:,:))) * zetaWeights;
    particleFlux_vE0_psiHat(ispecies) = (thetaWeights') * (squeeze(particleFluxBeforeSurfaceIntegral_vE0(ispecies,:,:))) * zetaWeights;
    particleFlux_vE_psiHat(ispecies) = (thetaWeights') * (squeeze(particleFluxBeforeSurfaceIntegral_vE(ispecies,:,:))) * zetaWeights;
    momentumFlux_vm0_psiHat(ispecies) = (thetaWeights') * (squeeze(momentumFluxBeforeSurfaceIntegral_vm0(ispecies,:,:))) * zetaWeights;
    momentumFlux_vm_psiHat(ispecies) = (thetaWeights') * (squeeze(momentumFluxBeforeSurfaceIntegral_vm(ispecies,:,:))) * zetaWeights;
    momentumFlux_vE0_psiHat(ispecies) = (thetaWeights') * (squeeze(momentumFluxBeforeSurfaceIntegral_vE0(ispecies,:,:))) * zetaWeights;
    momentumFlux_vE_psiHat(ispecies) = (thetaWeights') * (squeeze(momentumFluxBeforeSurfaceIntegral_vE(ispecies,:,:))) * zetaWeights;
    heatFlux_vm0_psiHat(ispecies) = (thetaWeights') * (squeeze(heatFluxBeforeSurfaceIntegral_vm0(ispecies,:,:))) * zetaWeights;
    heatFlux_vm_psiHat(ispecies) = (thetaWeights') * (squeeze(heatFluxBeforeSurfaceIntegral_vm(ispecies,:,:))) * zetaWeights;
    heatFlux_vE0_psiHat(ispecies) = (thetaWeights') * (squeeze(heatFluxBeforeSurfaceIntegral_vE0(ispecies,:,:))) * zetaWeights;
    heatFlux_vE_psiHat(ispecies) = (thetaWeights') * (squeeze(heatFluxBeforeSurfaceIntegral_vE(ispecies,:,:))) * zetaWeights;
    
    jHat = jHat + Z*squeeze(flow(ispecies,:,:));
    totalDensity(ispecies,:,:) = densityPerturbation(ispecies,:,:) + nHat;
    totalPressure(ispecies,:,:) = densityPerturbation(ispecies,:,:) + nHat*THat;
    velocityUsingFSADensity(ispecies,:,:) = flow(ispecies,:,:)/nHat;
    velocityUsingTotalDensity(ispecies,:,:) = flow(ispecies,:,:)./totalDensity(ispecies,:,:);
    MachUsingFSAThermalSpeed(ispecies,:,:) = velocityUsingFSADensity(ispecies,:,:)*sqrtm/sqrtT;
end

FSADensityPerturbation = FSADensityPerturbation / VPrimeHat;
FSABFlow = FSABFlow / VPrimeHat;
FSAPressurePerturbation = FSAPressurePerturbation / VPrimeHat;
FSABjHat = Zs(:)' * FSABFlow;
FSABjHatOverB0 = FSABjHat / B0OverBBar;
FSABjHatOverRootFSAB2 = FSABjHat / sqrt(FSABHat2);


for ispecies = 1:Nspecies
    fprintf('Results for species %d:\n',ispecies)
    fprintf('   FSADensityPerturbation: %g\n',FSADensityPerturbation(ispecies))
    fprintf('                 FSABFlow: %g\n',FSABFlow(ispecies))
    fprintf('       max and min Mach #: %g, %g\n',max(max(MachUsingFSAThermalSpeed(ispecies,:,:))), min(min(MachUsingFSAThermalSpeed(ispecies,:,:))))
    fprintf('  FSAPressurePerturbation: %g\n',FSAPressurePerturbation(ispecies))
    fprintf('  particleFlux_vm0_psiHat: %g\n',particleFlux_vm0_psiHat(ispecies))
    fprintf('   particleFlux_vm_psiHat: %g\n',particleFlux_vm_psiHat(ispecies))
    if includePhi1
        fprintf('  particleFlux_vE0_psiHat: %g\n',particleFlux_vE0_psiHat(ispecies))
        fprintf('   particleFlux_vE_psiHat: %g\n',particleFlux_vE_psiHat(ispecies))
        fprintf('  particleFlux_vd1_psiHat: %g\n',particleFlux_vd1_psiHat(ispecies))
        fprintf('   particleFlux_vd_psiHat: %g\n',particleFlux_vd_psiHat(ispecies))
    end
    fprintf('  momentumFlux_vm0_psiHat: %g\n',momentumFlux_vm0_psiHat(ispecies))
    fprintf('   momentumFlux_vm_psiHat: %g\n',momentumFlux_vm_psiHat(ispecies))
    if includePhi1
        fprintf('  momentumFlux_vE0_psiHat: %g\n',momentumFlux_vE0_psiHat(ispecies))
        fprintf('   momentumFlux_vE_psiHat: %g\n',momentumFlux_vE_psiHat(ispecies))
        fprintf('  momentumFlux_vd1_psiHat: %g\n',momentumFlux_vd1_psiHat(ispecies))
        fprintf('   momentumFlux_vd_psiHat: %g\n',momentumFlux_vd_psiHat(ispecies))
    end
    fprintf('      heatFlux_vm0_psiHat: %g\n',heatFlux_vm0_psiHat(ispecies))
    fprintf('       heatFlux_vm_psiHat: %g\n',heatFlux_vm_psiHat(ispecies))
    if includePhi1
        fprintf('      heatFlux_vE0_psiHat: %g\n',heatFlux_vE0_psiHat(ispecies))
        fprintf('       heatFlux_vE_psiHat: %g\n',heatFlux_vE_psiHat(ispecies))
        fprintf('      heatFlux_vd1_psiHat: %g\n',heatFlux_vd1_psiHat(ispecies))
        fprintf('       heatFlux_vd_psiHat: %g\n',heatFlux_vd_psiHat(ispecies))
    end
    if constraintScheme==1
        fprintf('          particle source: %g\n',sources(ispecies,1))
        fprintf('              heat source: %g\n',sources(ispecies,2))
    elseif constraintScheme==2
        fprintf('          max/min sources: %g ... %g\n',max(sources(ispecies,:)),min(sources(ispecies,:)))
    end
end
fprintf('   FSABjHat (j_bootstrap): %g\n',FSABjHat)
if includePhi1
    fprintf('                   lambda: %g\n',lambda)
end

if RHSMode>1
    fprintf('Transport matrix:\n')
    transportMatrix
end

end