function sfincs_initF0()

global f0 matrixSize x nHats mHats THats
global Nspecies Ntheta Nzeta Nx

f0 = zeros(matrixSize,1);
expx2 = exp(-x.*x);
allx = 1:Nx;

L = 0;
factors = nHats .* ((mHats./(pi*THats)) .^ (3/2));
for ispecies = 1:Nspecies
    factor = factors(ispecies)*expx2;
    for izeta = 1:Nzeta
        for itheta = 1:Ntheta
            indices = sfincs_indices(ispecies, allx, L+1, itheta, izeta, BLOCK_F);
            f0(indices) = factor;
        end
    end
end

end