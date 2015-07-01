function sfincs_computeBIntegrals()

global thetaWeights zetaWeights geometryScheme
global VPrimeHat FSABHat2 BHat GHat IHat B0OverBBar DHat
global BHat_sub_theta BHat_sub_zeta

VPrimeHat = thetaWeights' * (1./DHat) * zetaWeights;
FSABHat2 = (1/VPrimeHat) * thetaWeights' * (BHat.*BHat./DHat) * zetaWeights;

if geometryScheme==5
    FSABHat3 = (1/VPrimeHat) * thetaWeights' * (BHat.*BHat.*BHat./DHat) * zetaWeights;
    B0OverBBar = FSABHat3/FSABHat2;
    GHat = thetaWeights' * (BHat_sub_zeta) * zetaWeights / (4*pi*pi);
    IHat = thetaWeights' * (BHat_sub_theta) * zetaWeights / (4*pi*pi);
end

end