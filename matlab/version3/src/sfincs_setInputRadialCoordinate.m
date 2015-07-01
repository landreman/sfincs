function sfincs_setInputRadialCoordinate()

global inputRadialCoordinateForGradients
global psiHat_wish psiN_wish rHat_wish rN_wish
global psiHat psiN rHat rN
global dPhiHatdpsiHat dPhiHatdpsiN dPhiHatdrHat dPhiHatdrN
global dTHatdpsiHats dTHatdpsiNs dTHatdrHats dTHatdrNs
global dnHatdpsiHats dnHatdpsiNs dnHatdrHats dnHatdrNs
global psiAHat aHat
global ddpsiN2ddpsiHat ddrHat2ddpsiHat ddrN2ddpsiHat
global ddpsiHat2ddpsiN ddpsiHat2ddrHat ddpsiHat2ddrN

% By the time this function is called,
% rN should have been set by sfincs_computeBHat!  
if abs(rN+9999) < 1e-6
    error('It appears that rN was not set by sfincs_computeBHat()')
end
assert(rN >= 0)
assert(rN <= 1)

% Using rN, set the other radial variables:
psiHat = psiAHat * rN^2;
psiN = rN^2;
rHat = rN*aHat;

fprintf('Requested/actual flux surface for this calculation, in various radial coordinates:\n')
fprintf('  psiHat = %g / %g\n', psiHat_wish, psiHat)
fprintf('  psiN   = %g / %g\n', psiN_wish,   psiN)
fprintf('  rHat   = %g / %g\n', rHat_wish,   rHat)
fprintf('  rN     = %g / %g\n', rN_wish,     rN)

% Conversion factors
ddpsiN2ddpsiHat = 1/psiAHat;
ddrHat2ddpsiHat = aHat/(2*rN*psiAHat);
ddrN2ddpsiHat = 1/(2*rN*psiAHat);

ddpsiHat2ddpsiN = psiAHat;
ddpsiHat2ddrHat = (2*rN*psiAHat)/aHat;
ddpsiHat2ddrN = 2*rN*psiAHat;

% Set d/dpsiHat quantities:
switch inputRadialCoordinateForGradients
    case 0
        fprintf('Setting input gradients (of n, T, and Phi) from the specified dXXXdpsiHat values.\n')
        % Nothing more to do here.
    case 1
        fprintf('Setting input gradients (of n, T, and Phi) from the specified dXXXdpsiN values.\n')
        dPhiHatdpsiHat = dPhiHatdpsiN * ddpsiN2ddpsiHat;
        dnHatdpsiHats  = dnHatdpsiNs  * ddpsiN2ddpsiHat;
        dTHatdpsiHats  = dTHatdpsiNs  * ddpsiN2ddpsiHat;
    case 2
        fprintf('Setting input gradients (of n, T, and Phi) from the specified dXXXdrHat values.\n')
        dPhiHatdpsiHat = dPhiHatdrHat * ddrHat2ddpsiHat;
        dnHatdpsiHats  = dnHatdrHats  * ddrHat2ddpsiHat;
        dTHatdpsiHats  = dTHatdrHats  * ddrHat2ddpsiHat;
    case 3
        fprintf('Setting input gradients (of n, T, and Phi) from the specified dXXXdrN values.\n')
        dPhiHatdpsiHat = dPhiHatdrN * ddrN2ddpsiHat;
        dnHatdpsiHats  = dnHatdrNs  * ddrN2ddpsiHat;
        dTHatdpsiHats  = dTHatdrNs  * ddrN2ddpsiHat;
    otherwise
        error('Invalid inputRadialCoordinateForGradients')
end

% Finally, convert the input gradients from d/dpsiHat to all the other radial
% coordinates:

dPhiHatdpsiN = dPhiHatdpsiHat * ddpsiHat2ddpsiN;
dPhiHatdrHat = dPhiHatdpsiHat * ddpsiHat2ddrHat;
dPhiHatdrN   = dPhiHatdpsiHat * ddpsiHat2ddrN;

dnHatdpsiNs = dnHatdpsiHats * ddpsiHat2ddpsiN;
dnHatdrHats = dnHatdpsiHats * ddpsiHat2ddrHat;
dnHatdrNs   = dnHatdpsiHats * ddpsiHat2ddrN;

dTHatdpsiNs = dTHatdpsiHats * ddpsiHat2ddpsiN;
dTHatdrHats = dTHatdpsiHats * ddpsiHat2ddrHat;
dTHatdrNs   = dTHatdpsiHats * ddpsiHat2ddrN;

end