function sfincs_diagnostics()

global stateVector f0 constraintScheme sources force0RadialCurrentInEquilibrium
global Nalpha Nzeta Nspecies Nx Nxi alphaWeights zetaWeights RHSMode
global Phi1Hat dPhi1Hatdalpha dPhi1Hatdzeta ddalpha ddzeta includePhi1
global BLOCK_F BLOCK_QN BLOCK_PHI1_CONSTRAINT BLOCK_DENSITY_CONSTRAINT BLOCK_PRESSURE_CONSTRAINT BLOCK_F_CONSTRAINT indexVars
global BHat DHat BHat_sub_theta BHat_sub_zeta dBHatdtheta dBHatdzeta dBHat_sub_zeta_dtheta dBHat_sub_theta_dzeta
global THats Zs mHats nHats x xWeights xi xiWeights Delta gamma
global B0OverBBar FSABHat2 VPrimeHat
global GHat iota IHat transportMatrix whichRHS

global FSADensityPerturbation FSABFlow FSAPressurePerturbation
global particleFlux_vm0_psiHat particleFlux_vm_psiHat particleFlux_vE0_psiHat particleFlux_vE_psiHat particleFlux_vd_psiHat particleFlux_vd1_psiHat
global momentumFlux_vm0_psiHat momentumFlux_vm_psiHat momentumFlux_vE0_psiHat momentumFlux_vE_psiHat momentumFlux_vd_psiHat momentumFlux_vd1_psiHat
global heatFlux_vm0_psiHat heatFlux_vm_psiHat heatFlux_vE0_psiHat heatFlux_vE_psiHat heatFlux_vd_psiHat heatFlux_vd1_psiHat heatFlux_withoutPhi1_psiHat
global jHat FSABjHat FSABjHatOverB0 FSABjHatOverRootFSAB2
global totalDensity totalPressure velocityUsingFSADensity velocityUsingTotalDensity MachUsingFSAThermalSpeed

fprintf('Computing diagnostics.\n')

if includePhi1
    for ialpha = 1:Nalpha
        indices = sfincs_indices(1,1,1,ialpha,1:Nzeta,BLOCK_QN, indexVars);
        Phi1Hat(ialpha,:) = stateVector(indices);
    end
    
    dPhi1Hatdalpha = ddalpha * Phi1Hat;
    dPhi1Hatdzeta = (ddzeta * (Phi1Hat'))';
    
    index = sfincs_indices(1,1,1,1,1,BLOCK_PHI1_CONSTRAINT, indexVars);
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
            sources(ispecies,1) = stateVector(sfincs_indices(ispecies,1,1,1,1,BLOCK_DENSITY_CONSTRAINT, indexVars));
            sources(ispecies,2) = stateVector(sfincs_indices(ispecies,1,1,1,1,BLOCK_PRESSURE_CONSTRAINT, indexVars));
        end
    case 2
        sources = zeros(Nspecies,Nx);
        for ispecies = 1:Nspecies
            sources(ispecies,:) = stateVector(sfincs_indices(ispecies,1:Nx,1,1,1,BLOCK_F_CONSTRAINT, indexVars));
        end
    otherwise
        error('Invalid constraintScheme')
end

fullf = stateVector + f0;

jHat = zeros(Nalpha,Nzeta);
densityPerturbation = zeros(Nspecies,Nalpha,Nzeta);
totalDensity = zeros(Nspecies,Nalpha,Nzeta);
flow = zeros(Nspecies,Nalpha,Nzeta);
velocityUsingFSADensity = zeros(Nspecies,Nalpha,Nzeta);
velocityUsingTotalDensity = zeros(Nspecies,Nalpha,Nzeta);
MachUsingFSAThermalSpeed = zeros(Nspecies,Nalpha,Nzeta);
pressurePerturbation = zeros(Nspecies,Nalpha,Nzeta);
totalPressure = zeros(Nspecies,Nalpha,Nzeta);
particleFluxBeforeSurfaceIntegral_vm0 = zeros(Nspecies,Nalpha,Nzeta);
particleFluxBeforeSurfaceIntegral_vm = zeros(Nspecies,Nalpha,Nzeta);
particleFluxBeforeSurfaceIntegral_vE0 = zeros(Nspecies,Nalpha,Nzeta);
particleFluxBeforeSurfaceIntegral_vE = zeros(Nspecies,Nalpha,Nzeta);
momentumFluxBeforeSurfaceIntegral_vm0 = zeros(Nspecies,Nalpha,Nzeta);
momentumFluxBeforeSurfaceIntegral_vm = zeros(Nspecies,Nalpha,Nzeta);
momentumFluxBeforeSurfaceIntegral_vE0 = zeros(Nspecies,Nalpha,Nzeta);
momentumFluxBeforeSurfaceIntegral_vE = zeros(Nspecies,Nalpha,Nzeta);
heatFluxBeforeSurfaceIntegral_vm0 = zeros(Nspecies,Nalpha,Nzeta);
heatFluxBeforeSurfaceIntegral_vm = zeros(Nspecies,Nalpha,Nzeta);
heatFluxBeforeSurfaceIntegral_vE0 = zeros(Nspecies,Nalpha,Nzeta);
heatFluxBeforeSurfaceIntegral_vE = zeros(Nspecies,Nalpha,Nzeta);

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
    factor2 = zeros(Nalpha,Nzeta);
else
    factor2 = 2 * (dBHat_sub_zeta_dtheta-dBHat_sub_theta_dzeta) ./ (BHat.*BHat);
end
factor_vE = (BHat_sub_theta.*dPhi1Hatdzeta - BHat_sub_zeta.*dPhi1Hatdalpha) ./ (BHat.*BHat);

for ispecies = 1:Nspecies
    THat = THats(ispecies);
    nHat = nHats(ispecies);
    mHat = mHats(ispecies);
    Z = Zs(ispecies);
    sqrtT = sqrt(THat);
    sqrtm = sqrt(mHat);
    
    densityWeights = 2*pi*nHat*(x.*x.*xWeights)';
    flowWeights = 2*pi*nHat*sqrtT/(sqrtm)*(x.*x.*x.*xWeights)';
    pressureWeights = 4*pi*THat*nHat/(3)*(x.*x.*x.*x.*xWeights)';
    particleFluxWeights_vm = pi*Delta*THat/(Z*VPrimeHat)*(x.*x.*x.*x.*xWeights)';
    particleFluxWeights_vE = 2*gamma*pi*Delta/(VPrimeHat)*(x.*x.*xWeights)';
    momentumFluxWeights_vm = pi*Delta*THat*sqrtT*sqrtm/(Z*VPrimeHat)*(x.*x.*x.*x.*x.*xWeights)';
    momentumFluxWeights_vE = 2*gamma*pi*Delta*sqrtT*sqrtm/(VPrimeHat)*(x.*x.*x.*xWeights)';
    heatFluxWeights_vm = pi*Delta*THat*THat/(2*Z*VPrimeHat)*(x.*x.*x.*x.*x.*x.*xWeights)';
    heatFluxWeights_vE = 2*gamma*pi*Delta*THat/(2*VPrimeHat)*(x.*x.*x.*x.*xWeights)';
        
    % Do all the integrals over x and xi:
    for ialpha = 1:Nalpha
        for izeta = 1:Nzeta
            for ixi = 1:Nxi
                indices = sfincs_indices(ispecies, 1:Nx, ixi, ialpha, izeta, BLOCK_F, indexVars);
                
                densityPerturbation(ispecies, ialpha, izeta) = densityPerturbation(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*densityWeights * stateVector(indices);
                
                flow(ispecies, ialpha, izeta) = flow(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*xi(ixi) * flowWeights * stateVector(indices);
                
                pressurePerturbation(ispecies, ialpha, izeta) = pressurePerturbation(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*pressureWeights * stateVector(indices);
                
                particleFluxBeforeSurfaceIntegral_vm0(ispecies, ialpha, izeta) = particleFluxBeforeSurfaceIntegral_vm0(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*(factor(ialpha,izeta)*(1+xi(ixi)*xi(ixi))+factor2(ialpha,izeta)*xi(ixi)*xi(ixi))*particleFluxWeights_vm * f0(indices);
                
                particleFluxBeforeSurfaceIntegral_vm(ispecies, ialpha, izeta) = particleFluxBeforeSurfaceIntegral_vm(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*(factor(ialpha,izeta)*(1+xi(ixi)*xi(ixi))+factor2(ialpha,izeta)*xi(ixi)*xi(ixi))*particleFluxWeights_vm * fullf(indices);
                
                particleFluxBeforeSurfaceIntegral_vE0(ispecies, ialpha, izeta) = particleFluxBeforeSurfaceIntegral_vE0(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*factor_vE(ialpha,izeta) * particleFluxWeights_vE * f0(indices);
                
                particleFluxBeforeSurfaceIntegral_vE(ispecies, ialpha, izeta) = particleFluxBeforeSurfaceIntegral_vE(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*factor_vE(ialpha,izeta) * particleFluxWeights_vE * fullf(indices);
                
                heatFluxBeforeSurfaceIntegral_vm0(ispecies, ialpha, izeta) = heatFluxBeforeSurfaceIntegral_vm0(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*(factor(ialpha,izeta)*(1+xi(ixi)*xi(ixi))+factor2(ialpha,izeta)*xi(ixi)*xi(ixi))*heatFluxWeights_vm * f0(indices);
                
                heatFluxBeforeSurfaceIntegral_vm(ispecies, ialpha, izeta) = heatFluxBeforeSurfaceIntegral_vm(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*(factor(ialpha,izeta)*(1+xi(ixi)*xi(ixi))+factor2(ialpha,izeta)*xi(ixi)*xi(ixi))*heatFluxWeights_vm * fullf(indices);
                
                heatFluxBeforeSurfaceIntegral_vE0(ispecies, ialpha, izeta) = heatFluxBeforeSurfaceIntegral_vE0(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*factor_vE(ialpha,izeta) * heatFluxWeights_vE * f0(indices);
                
                heatFluxBeforeSurfaceIntegral_vE(ispecies, ialpha, izeta) = heatFluxBeforeSurfaceIntegral_vE(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*factor_vE(ialpha,izeta) * heatFluxWeights_vE * fullf(indices);
                
                momentumFluxBeforeSurfaceIntegral_vm0(ispecies, ialpha, izeta) = momentumFluxBeforeSurfaceIntegral_vm0(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*xi(ixi)*(factor(ialpha,izeta)*(1+xi(ixi)*xi(ixi))+factor2(ialpha,izeta)*xi(ixi)*xi(ixi))*BHat(ialpha,izeta)*momentumFluxWeights_vm * f0(indices);
                
                momentumFluxBeforeSurfaceIntegral_vm(ispecies, ialpha, izeta) = momentumFluxBeforeSurfaceIntegral_vm(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*xi(ixi)*(factor(ialpha,izeta)*(1+xi(ixi)*xi(ixi))+factor2(ialpha,izeta)*xi(ixi)*xi(ixi))*BHat(ialpha,izeta)*momentumFluxWeights_vm * fullf(indices);
                
                momentumFluxBeforeSurfaceIntegral_vE0(ispecies, ialpha, izeta) = momentumFluxBeforeSurfaceIntegral_vE0(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*xi(ixi)*factor_vE(ialpha,izeta)*BHat(ialpha,izeta)*momentumFluxWeights_vE * f0(indices);
                
                momentumFluxBeforeSurfaceIntegral_vE(ispecies, ialpha, izeta) = momentumFluxBeforeSurfaceIntegral_vE(ispecies, ialpha, izeta) + ...
                    xiWeights(ixi)*xi(ixi)*factor_vE(ialpha,izeta)*BHat(ialpha,izeta)*momentumFluxWeights_vE * fullf(indices);
            end
        end
    end
    
    % Done with the integrals over x and xi.
    % Now do the integrals over alpha and zeta.
   
    % We need to do this next step because Matlab thinks
    % squeeze(densityPerturbation(ispecies,:,:)) and DHat have different
    % dimensions when Nzeta=1.
    if Nzeta==1
        DHat_x = DHat';
        BHat_x = BHat';
    else
        DHat_x = DHat;
        BHat_x = BHat;
    end
    
    FSADensityPerturbation(ispecies) = dot(alphaWeights, (squeeze(densityPerturbation(ispecies,:,:)) ./ DHat_x) * zetaWeights);
    FSABFlow(ispecies) = dot(alphaWeights, (squeeze(flow(ispecies,:,:)) .* BHat_x./ DHat_x) * zetaWeights);
    FSAPressurePerturbation(ispecies) = dot(alphaWeights, (squeeze(pressurePerturbation(ispecies,:,:)) ./ DHat_x) * zetaWeights);
    particleFlux_vm0_psiHat(ispecies) = dot(alphaWeights, (squeeze(particleFluxBeforeSurfaceIntegral_vm0(ispecies,:,:))) * zetaWeights);
    particleFlux_vm_psiHat(ispecies) = dot(alphaWeights, (squeeze(particleFluxBeforeSurfaceIntegral_vm(ispecies,:,:))) * zetaWeights);
    particleFlux_vE0_psiHat(ispecies) = dot(alphaWeights, (squeeze(particleFluxBeforeSurfaceIntegral_vE0(ispecies,:,:))) * zetaWeights);
    particleFlux_vE_psiHat(ispecies) = dot(alphaWeights, (squeeze(particleFluxBeforeSurfaceIntegral_vE(ispecies,:,:))) * zetaWeights);
    momentumFlux_vm0_psiHat(ispecies) = dot(alphaWeights, (squeeze(momentumFluxBeforeSurfaceIntegral_vm0(ispecies,:,:))) * zetaWeights);
    momentumFlux_vm_psiHat(ispecies) = dot(alphaWeights, (squeeze(momentumFluxBeforeSurfaceIntegral_vm(ispecies,:,:))) * zetaWeights);
    momentumFlux_vE0_psiHat(ispecies) = dot(alphaWeights, (squeeze(momentumFluxBeforeSurfaceIntegral_vE0(ispecies,:,:))) * zetaWeights);
    momentumFlux_vE_psiHat(ispecies) = dot(alphaWeights, (squeeze(momentumFluxBeforeSurfaceIntegral_vE(ispecies,:,:))) * zetaWeights);
    heatFlux_vm0_psiHat(ispecies) = dot(alphaWeights, (squeeze(heatFluxBeforeSurfaceIntegral_vm0(ispecies,:,:))) * zetaWeights);
    heatFlux_vm_psiHat(ispecies) = dot(alphaWeights, (squeeze(heatFluxBeforeSurfaceIntegral_vm(ispecies,:,:))) * zetaWeights);
    heatFlux_vE0_psiHat(ispecies) = dot(alphaWeights, (squeeze(heatFluxBeforeSurfaceIntegral_vE0(ispecies,:,:))) * zetaWeights);
    heatFlux_vE_psiHat(ispecies) = dot(alphaWeights, (squeeze(heatFluxBeforeSurfaceIntegral_vE(ispecies,:,:))) * zetaWeights);
    
    jHat = jHat + reshape(Z*squeeze(flow(ispecies,:,:)),[Nalpha,Nzeta]);
    totalDensity(ispecies,:,:) = densityPerturbation(ispecies,:,:) + nHat;
    totalPressure(ispecies,:,:) = densityPerturbation(ispecies,:,:) + nHat*THat;
    velocityUsingFSADensity(ispecies,:,:) = flow(ispecies,:,:)/nHat;
    velocityUsingTotalDensity(ispecies,:,:) = flow(ispecies,:,:)./totalDensity(ispecies,:,:);
    MachUsingFSAThermalSpeed(ispecies,:,:) = velocityUsingFSADensity(ispecies,:,:)*sqrtm/sqrtT;
end

particleFlux_vd_psiHat = particleFlux_vm_psiHat + particleFlux_vE_psiHat;
momentumFlux_vd_psiHat = momentumFlux_vm_psiHat + momentumFlux_vE_psiHat;
heatFlux_vd_psiHat = heatFlux_vm_psiHat + heatFlux_vE_psiHat;

particleFlux_vd1_psiHat = particleFlux_vm_psiHat + particleFlux_vE0_psiHat;
momentumFlux_vd1_psiHat = momentumFlux_vm_psiHat + momentumFlux_vE0_psiHat;
heatFlux_vd1_psiHat = heatFlux_vm_psiHat + heatFlux_vE0_psiHat;

heatFlux_withoutPhi1_psiHat = heatFlux_vm_psiHat + (5/3)*heatFlux_vE0_psiHat;

FSADensityPerturbation = FSADensityPerturbation / VPrimeHat;
FSABFlow = FSABFlow / VPrimeHat;
FSAPressurePerturbation = FSAPressurePerturbation / VPrimeHat;
FSABjHat = Zs(:)' * FSABFlow;
FSABjHatOverB0 = FSABjHat / B0OverBBar;
FSABjHatOverRootFSAB2 = FSABjHat / sqrt(FSABHat2);

if RHSMode==3
    % Monoenergetic transport matrix
    ispecies = 1;
    nHat = nHats(ispecies);
    THat = THats(ispecies);
    mHat = mHats(ispecies);
    sqrtTHat = sqrt(THat);
    sqrtmHat = sqrt(mHat);
    
    % The factors of THat, mHat, nHat, and Z are unnecessary below (since all are 1).
    switch whichRHS
        case 1
            transportMatrix(1,1) = 4/(Delta*Delta)*sqrtTHat/sqrtmHat*Zs(1)*Zs(1)*(GHat+iota*IHat) ...
                *particleFlux_vm_psiHat(1)*B0OverBBar/(THat*THat*GHat*GHat);
            transportMatrix(2,1) = 2*Zs(1)*FSABFlow(1)/(Delta*GHat*THat);
        case 2
            transportMatrix(1,2) = particleFlux_vm_psiHat(1)*2*FSABHat2/(nHat*gamma*Delta*GHat);
            transportMatrix(2,2) = FSABFlow(1)*sqrtTHat*sqrtmHat*FSABHat2/((GHat+iota*IHat)*gamma*Zs(1)*nHat*B0OverBBar);
    end
end

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
    format longg
    transportMatrix
end

end
