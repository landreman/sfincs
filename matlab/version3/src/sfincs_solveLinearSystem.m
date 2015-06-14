function soln = sfincs_solveLinearSystem(matrix, rhs, preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q)

global GMRES_restart GMRES_maxIterations solverTolerance useIterativeLinearSolver

if useIterativeLinearSolver
    tic
    fprintf('Beginning direct linear solve... ')
    soln = matrix \ rhs;
    fprintf('Done. Took %g seconds.\n',toc)
    
else
    tic
    fprintf('Applying GMRES... ')
    x0 = zeros(size(rhs));
    [soln,fl0,rr0,it0,rv0]=gmres(matrix,rhs,GMRES_restart,solverTolerance,GMRES_maxIterations/GMRES_restart,@preconditioner, [], x0);
    fprintf('Done. Took %g seconds.\n',toc)
    switch fl0
        case 0
            fprintf('Converged!\n')
            didItConverge = true;
        case 1
            fprintf('Did not converge :(\n')
            didItConverge = false;
        case 2
            fprintf('Preconditioner was ill-conditioned\n')
            didItConverge = false;
        case 3
            fprintf('Stagnated :(\n')
            didItConverge = false;
        otherwise
            error('Code should not reach this point')
    end
    residualVector = matrix*soln - rhs;
    L1_normalization = sum(abs(rhs));
    L2_normalization = sqrt(sum(rhs.*rhs));
    L_inf_normalization = max(abs(rhs));
    fprintf('Absolute GMRES residual in L1 norm: %g,   L2 norm: %g,   L_inf norm: %g\n',sum(abs(residualVector)), sum(residualVector.*residualVector), max(abs(residualVector)))
    fprintf('Relative GMRES residual in L1 norm: %g,   L2 norm: %g,   L_inf norm: %g\n',sum(abs(residualVector))/L1_normalization, sum(residualVector.*residualVector)/L2_normalization, max(abs(residualVector))/L_inf_normalization)
    
    figure(10)
    clf
    semilogy(rv0/rv0(1),'-o');
    xlabel('Iteration number');
    ylabel('Relative residual');
    title('Convergence of GMRES');
    drawnow
end

    function solnVector = preconditioner(rhsVector)
        solnVector = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * rhsVector)));
    end

end