function indices = sfincs_indices(ispecies, ix, ixi, itheta, izeta, whichBlock, v)

% Note: In contrast to the fortran version's indices.F90, which returns 0-based
% indices for PETSc, this function returns 1-based indices for matlab.

% The 'v' input structure is used instead of directly using the global
% variables because it seems extremely slow to access the global variables.

% You can un-comment the validation assertions below for debugging.
% They are commented out because they slow down this function a lot,
% slowing down sfincs substantially.

%{
% Validation:
global matrixSize 
%global constraintScheme
assert(all(ispecies>0))
assert(all(itheta>0))
assert(all(izeta>0))
assert(all(ixi>0))
assert(all(ix>0))

assert(all(ispecies <= v.Nspecies))
assert(all(itheta <= v.Ntheta))
assert(all(izeta <= v.Nzeta))
assert(all(ixi <= v.Nxi))
assert(all(ix <= v.Nx))
%}

switch whichBlock
    case v.BLOCK_F
        indices = (ispecies-1)*v.Nx*v.Nxi*v.Ntheta*v.Nzeta ...
            +(ix-1)*v.Nxi*v.Ntheta*v.Nzeta ...
            +(ixi-1)*v.Ntheta*v.Nzeta ...
            +(itheta-1)*v.Nzeta ...
            +izeta;
        
    case v.BLOCK_QN
        %{
        if ~v.includePhi1
            error('whichBlock cannot be v.BLOCK_QN if v.includePhi1 is false.')
        end
        %}
        indices = v.Nspecies*v.Nx*v.Nxi*v.Ntheta*v.Nzeta ...
            +(itheta-1)*v.Nzeta ...
            +izeta;
        
    case v.BLOCK_PHI1_CONSTRAINT
        %{
        if ~v.includePhi1
            error('whichBlock cannot be v.BLOCK_QN if v.includePhi1 is false.')
        end
        %}
        indices = v.Nspecies*v.Nx*v.Nxi*v.Ntheta*v.Nzeta ...
            + v.Ntheta*v.Nzeta + 1;
        
    case v.BLOCK_DENSITY_CONSTRAINT
        %{
        if constraintScheme ~= 1
            error('whichBlock is v.BLOCK_DENSITY_CONSTRAINT but constraintScheme = %d',constraintScheme)
        end
        %}
        indices = v.Nspecies*v.Nx*v.Nxi*v.Ntheta*v.Nzeta ...
            + (ispecies-1)*2 + 1;
        if v.includePhi1
            indices = indices + v.Ntheta*v.Nzeta + 1;
        end
        
    case v.BLOCK_PRESSURE_CONSTRAINT
        %{
        if constraintScheme ~= 1
            error('whichBlock is v.BLOCK_PRESSURE_CONSTRAINT but constraintScheme = %d',constraintScheme)
        end
        %}
        indices = v.Nspecies*v.Nx*v.Nxi*v.Ntheta*v.Nzeta ...
            + (ispecies-1)*2 + 2;
        if v.includePhi1
            indices = indices + v.Ntheta*v.Nzeta + 1;
        end
        
    case v.BLOCK_F_CONSTRAINT
        %{
        if constraintScheme ~= 2
            error('whichBlock is v.BLOCK_F_CONSTRAINT but constraintScheme = %d',constraintScheme)
        end
        %}
        indices = v.Nspecies*v.Nx*v.Nxi*v.Ntheta*v.Nzeta ...
            + (ispecies-1)*v.Nx + ix;
        if v.includePhi1
            indices = indices + v.Ntheta*v.Nzeta + 1;
        end
        
    otherwise
        error('Invalid whichBlock')
end

%{
% Final validation:
assert(all(indices>0))
assert(all(indices <= matrixSize))
%}

end