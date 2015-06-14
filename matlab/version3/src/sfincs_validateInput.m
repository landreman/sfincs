function sfincs_validateInput()

global Zs mHats THats nHats Nspecies
global inputRadialCoordinateForGradients
global dTHatdpsiHats dnHatdpsiHats
global dTHatdpsiNs dnHatdpsiNs
global dTHatdrHats dnHatdrHats
global dTHatdrNs dnHatdrNs

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

end