function sfincs_computeMatrixSize()

global Ntheta Nzeta Nxi Nx Nspecies includePhi1 constraintScheme matrixSize

matrixSize = Ntheta * Nzeta * Nxi * Nx * Nspecies;

if includePhi1
    matrixSize = matrixSize + Ntheta*Nzeta + 1;
end

switch constraintScheme
    case 0
    case 1
        matrixSize = matrixSize + 2 * Nspecies;
    case 2
        matrixSize = matrixSize + Nx * Nspecies;
    otherwise
        error('Invalid constraintScheme')
end

fprintf('The matrix is %d x %d elements.\n',matrixSize, matrixSize)

end