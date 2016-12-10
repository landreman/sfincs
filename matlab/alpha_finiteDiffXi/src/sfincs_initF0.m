function sfincs_initF0()

global f0 matrixSize x
global Nspecies Nalpha Nzeta Nx Nxi
global BLOCK_F indexVars

f0 = zeros(matrixSize,1);
expx2 = exp(-x.*x);
allx = 1:Nx;

factor = 1/(pi*sqrt(pi));
for ispecies = 1:Nspecies
    factors = expx2 * factor;
    for izeta = 1:Nzeta
        for ialpha = 1:Nalpha
            for ixi = 1:Nxi
                indices = sfincs_indices(ispecies, allx, ixi, ialpha, izeta, BLOCK_F, indexVars);
                f0(indices) = factors;
            end
        end
    end
end

end
