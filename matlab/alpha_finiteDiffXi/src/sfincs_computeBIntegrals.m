function sfincs_computeBIntegrals()

global alphaWeights zetaWeights geometryScheme
global VPrimeHat FSABHat2 BHat GHat IHat B0OverBBar DHat
global BHat_sub_theta BHat_sub_zeta

VPrimeHat = alphaWeights' * (1./DHat) * zetaWeights;
FSABHat2 = (1/VPrimeHat) * alphaWeights' * (BHat.*BHat./DHat) * zetaWeights;

if geometryScheme==5
    FSABHat3 = (1/VPrimeHat) * alphaWeights' * (BHat.*BHat.*BHat./DHat) * zetaWeights;
    B0OverBBar = FSABHat3/FSABHat2;
    GHat = alphaWeights' * (BHat_sub_zeta) * zetaWeights / (4*pi*pi);
    IHat = alphaWeights' * (BHat_sub_theta) * zetaWeights / (4*pi*pi);
end

end
