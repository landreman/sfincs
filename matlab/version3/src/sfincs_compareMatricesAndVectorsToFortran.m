function sfincs_compareMatricesAndVectorsToFortran(directory)

global Jacobian preconditionerMatrix initialResidual stateVector

filename = fullfile(directory,'sfincsBinary_iteration_000_residual');
fprintf('Attempting to read %s\n',filename)
residual_fortran = PetscBinaryRead(filename);

filename = fullfile(directory,'sfincsBinary_iteration_000_whichMatrix_0');
fprintf('Attempting to read %s\n',filename)
preconditionerMatrix_fortran = PetscBinaryRead(filename);

filename = fullfile(directory,'sfincsBinary_iteration_000_whichMatrix_1');
fprintf('Attempting to read %s\n',filename)
Jacobian_fortran = PetscBinaryRead(filename);

filename = fullfile(directory,'sfincsBinary_iteration_000_stateVector');
fprintf('Attempting to read %s\n',filename)
stateVector_fortran = PetscBinaryRead(filename);


figure(4)
clf
numRows = 2;
numCols = 2;

subplot(numRows,numCols,1)
plot(initialResidual,'.-','DisplayName','matlab')
hold on
plot(residual_fortran,'x:r','DisplayName','fortran')
legend show
title('Initial residual vector')

subplot(numRows,numCols,3)
plot(initialResidual - residual_fortran,'.-m')
title('Differences in initial residual vector')

subplot(numRows,numCols,2)
plot(stateVector,'.-','DisplayName','matlab')
hold on
plot(stateVector_fortran,'x:r','DisplayName','fortran')
legend show
title('stateVector')

subplot(numRows,numCols,4)
plot(stateVector - stateVector_fortran,'.-m')
title('Differences in stateVector')

figure(5)
clf

numRows = 2;
numCols = 4;
th1 = 1e-12;
th2 = 1e-6;

subplot(numRows,numCols,1)
spy(abs(Jacobian)>th1)
title('Matlab Jacobian')

subplot(numRows,numCols,2)
spy(abs(Jacobian_fortran)>th1)
title('Fortran Jacobian')

subplot(numRows,numCols,3)
spy(abs(Jacobian_fortran-Jacobian)>th1)
title(['Differences > ',num2str(th1)])

subplot(numRows,numCols,4)
spy(abs(Jacobian_fortran-Jacobian)>th2)
title(['Differences > ',num2str(th2)])

subplot(numRows,numCols,5)
spy(abs(preconditionerMatrix)>th1)
title('Matlab preconditioner')

subplot(numRows,numCols,6)
spy(abs(preconditionerMatrix_fortran)>th1)
title('Fortran preconditioner')

subplot(numRows,numCols,7)
spy(abs(preconditionerMatrix_fortran-preconditionerMatrix)>th1)
title(['Differences > ',num2str(th1)])

subplot(numRows,numCols,8)
spy(abs(preconditionerMatrix_fortran-preconditionerMatrix)>th2)
title(['Differences > ',num2str(th2)])

end