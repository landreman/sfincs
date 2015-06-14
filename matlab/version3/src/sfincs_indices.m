function indices = sfincs_indices(ispecies, ix, ixi, itheta, izeta, whichBlock)

% Note: In contrast to the fortran version's indices.F90, which returns 0-based
% indices for PETSc, this function returns 1-based indices for matlab.

global Ntheta Nzeta Nxi Nx Nspecies matrixSize includePhi1 constraintScheme
global BLOCK_F BLOCK_QN BLOCK_PHI1_CONSTRAINT BLOCK_DENSITY_CONSTRAINT BLOCK_PRESSURE_CONSTRAINT BLOCK_F_CONSTRAINT

% Validation:
assert(all(ispecies>0))
assert(all(itheta>0))
assert(all(izeta>0))
assert(all(ixi>0))
assert(all(ix>0))

assert(all(ispecies <= Nspecies))
assert(all(itheta <= Ntheta))
assert(all(izeta <= Nzeta))
assert(all(ixi <= Nxi))
assert(all(ix <= Nx))

switch whichBlock
    case BLOCK_F
        indices = (ispecies-1)*Nx*Nxi*Ntheta*Nzeta ...
            +(ix-1)*Nxi*Ntheta*Nzeta ...
            +(ixi-1)*Ntheta*Nzeta ...
            +(itheta-1)*Nzeta ...
            +izeta;
        
    case BLOCK_QN
        if ~includePhi1
            error('whichBlock cannot be BLOCK_QN if includePhi1 is false.')
        end
        indices = Nspecies*Nx*Nxi*Ntheta*Nzeta ...
            +(itheta-1)*Nzeta ...
            +izeta;
        
    case BLOCK_PHI1_CONSTRAINT
        if ~includePhi1
            error('whichBlock cannot be BLOCK_QN if includePhi1 is false.')
        end
        indices = Nspecies*Nx*Nxi*Ntheta*Nzeta ...
            + Ntheta*Nzeta + 1;
        
    case BLOCK_DENSITY_CONSTRAINT
        if constraintScheme ~= 1
            error('whichBlock is BLOCK_DENSITY_CONSTRAINT but constraintScheme = %d',constraintScheme)
        end
        
        indices = Nspecies*Nx*Nxi*Ntheta*Nzeta ...
            + (ispecies-1)*2 + 1;
        if includePhi1
            indices = indices + Ntheta*Nzeta + 1;
        end
        
    case BLOCK_PRESSURE_CONSTRAINT
        if constraintScheme ~= 1
            error('whichBlock is BLOCK_PRESSURE_CONSTRAINT but constraintScheme = %d',constraintScheme)
        end
        
        indices = Nspecies*Nx*Nxi*Ntheta*Nzeta ...
            + (ispecies-1)*2 + 2;
        if includePhi1
            indices = indices + Ntheta*Nzeta + 1;
        end
        
    case BLOCK_F_CONSTRAINT
        if constraintScheme ~= 2
            error('whichBlock is BLOCK_F_CONSTRAINT but constraintScheme = %d',constraintScheme)
        end
        
        indices = Nspecies*Nx*Nxi*Ntheta*Nzeta ...
            + (ispecies-1)*Nx + ix;
        if includePhi1
            indices = indices + Ntheta*Nzeta + 1;
        end
        
    otherwise
        error('Invalid whichBlock')
end

% Final validation:
assert(all(indices>0))
assert(all(indices <= matrixSize))

end