function sfincs_validateInput()

global Zs mHats THats nHats Nspecies
global inputRadialCoordinateForGradients
global dTHatdpsiHats dnHatdpsiHats
global dTHatdpsiNs dnHatdpsiNs
global dTHatdrHats dnHatdrHats
global dTHatdrNs dnHatdrNs
global collisionOperator Nx
global includeXDotTerm includeElectricFieldTermInXiDot useDKESExBDrift
global includePhi1 RHSMode

Nspecies = numel(Zs);
assert(Nspecies == numel(mHats))
assert(Nspecies == numel(THats))
assert(Nspecies == numel(nHats))
switch inputRadialCoordinateForGradients
    case 0
        assert(Nspecies == numel(dTHatdpsiHats))
        assert(Nspecies == numel(dnHatdpsiHats))
    case 1
        assert(Nspecies == numel(dTHatdpsiNs))
        assert(Nspecies == numel(dnHatdpsiNs))
    case 2
        assert(Nspecies == numel(dTHatdrHats))
        assert(Nspecies == numel(dnHatdrHats))
    case 3
        assert(Nspecies == numel(dTHatdrNs))
        assert(Nspecies == numel(dnHatdrNs))
    otherwise
        error('Invalid inputRadialCoordinateForGradients')
end

if RHSMode==3
    if Nx ~= 1
        fprintf('WARNING: Since RHSMode==3, Nx will be set to 1.\n')
        Nx = 1;
    end
    if collisionOperator ~= 1
        fprintf('WARNING: Since RHSMode==3, collisionOperator will be set to 1.\n')
        collisionOperator = 1;
    end
    if includePhi1
        fprintf('WARNING: Since RHSMode==3, includePhi1 will be set to false.\n')
        includePhi1 = false;
    end
    if ~useDKESExBDrift
        fprintf('WARNING: Since RHSMode==3, useDKESExBDrift will be set to true.\n')
        useDKESExBDrift = true;
    end
    if includeXDotTerm
        fprintf('WARNING: Since RHSMode==3, includeXDotTerm will be set to false.\n')
        includeXDotTerm = false;
    end
    if includeElectricFieldTermInXiDot
        fprintf('WARNING: Since RHSMode==3, includeElectricFieldTermInXiDot will be set to false.\n')
        includeElectricFieldTermInXiDot = false;
    end
end

global quasineutralityOption withAdiabatic
if includePhi1 && (quasineutralityOption==2) && (~withAdiabatic)
    error('If includePhi1=true and quasineutralityOption==2, you must include an adiabatic species')
end

end
