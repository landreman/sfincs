function [soln, didItConverge, residual] ...
    = m20130728_01_solveLinearSystem(matrix, rhs, preconditionerHandle, ...
    tryIterativeSolvers, orderOfSolversToTry, requestedAbsoluteTolerance, maxIterations, restart, ...
    figureNum, tryDirectSolverIfIterativeSolversFail)
% Allowed values for orderOfSolversToTry:
% 1 = GMRES
% 2 = BiCGStab
% 3 = BiCGStab(l)
% 4 = TFQMR
% 5 = CGS
%
% rhs may have >1 column. soln will have the same size as rhs.
% If figureNum <= 0, no figure will be generated.

timeWhenSolveStarted = tic;

numSolversToTry = numel(orderOfSolversToTry);

if ~tryIterativeSolvers
    
    fprintf('Direct solution\n')
    soln = matrix \ rhs;
    didItConverge = true;
    residual = 0;
    
else
    % Use an iterative Krylov-space solver.
    
    soln = zeros(size(rhs));
    numCols = size(rhs,2);
    
    for col=1:numCols
        attempt=0;
        keepTrying = true;
        x0 = zeros(size(rhs));
        relativeTolerance = requestedAbsoluteTolerance;
        while keepTrying
            attempt = attempt+1;
            tic
            switch orderOfSolversToTry(attempt)
                case 1
                    solverName = 'GMRES';
                    fprintf('GMRES solution ***********\n')
                    [soln0,fl0,rr0,it0,rv0]=gmres(matrix,rhs(:,col),restart,relativeTolerance,maxIterations/restart,preconditionerHandle, [], x0);
                case 2
                    solverName = 'BiCGStab';
                    fprintf('BiCGStab solution ***********\n')
                    [soln0,fl0,rr0,it0,rv0]=bicgstab(matrix,rhs(:,col),relativeTolerance,maxIterations,preconditionerHandle, [], x0);
                case 3
                    solverName = 'BiCGStab(l)';
                    fprintf('BiCGStab(l) solution ***********\n')
                    [soln0,fl0,rr0,it0,rv0]=bicgstabl(matrix,rhs(:,col),relativeTolerance,maxIterations,preconditionerHandle, [], x0);
                case 4
                    solverName = 'TFQMR';
                    fprintf('TFQMR solution ***********\n')
                    [soln0,fl0,rr0,it0,rv0]=tfqmr(matrix,rhs(:,col),relativeTolerance,maxIterations,preconditionerHandle, [], x0);
                case 5
                    solverName = 'CGS';
                    fprintf('CGS solution ***********\n')
                    [soln0,fl0,rr0,it0,rv0]=cgs(matrix,rhs(:,col),relativeTolerance,maxIterations,preconditionerHandle, [], x0);
                otherwise
                    error('Invalid entry in orderOfSolversToTry')
                    
            end
            residualVector = matrix*soln0 - rhs;
            L1_normalization = sum(abs(rhs));
            L2_normalization = sqrt(sum(rhs.*rhs));
            L_inf_normalization = max(abs(rhs));
            fprintf('Absolute residual in L1 norm: %g,   L2 norm: %g,   L_inf norm: %g\n',sum(abs(residualVector)), sum(residualVector.*residualVector), max(abs(residualVector)))
            fprintf('Relative residual in L1 norm: %g,   L2 norm: %g,   L_inf norm: %g\n',sum(abs(residualVector))/L1_normalization, sum(residualVector.*residualVector)/L2_normalization, max(abs(residualVector))/L_inf_normalization)
            fprintf('Time to apply Krylov solver: %g seconds.\n',toc)
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
            if figureNum > 0
                figure(figureNum)
                if attempt==1
                    clf
                end
                subplot(1,numSolversToTry,attempt)
                semilogy(rv0/rv0(1),'-o');
                xlabel('Iteration number');
                ylabel('Relative residual');
                title(['Convergence of Krylov solver ',solverName]);
                drawnow
            end
            residual = min(rv0)/rv0(1);
            fprintf('minimum residual: %g.\n',residual)
            if fl0==0
                keepTrying=false;
                soln(:,col) = soln0;
            else
                if attempt >= numel(orderOfSolversToTry)
                    keepTrying=false;
                else
                    x0 = soln0;
                    fprintf('Iterative solver failed, so trying again with backup solver.\n')
                    % Adjust the tolerance based on the progress made so
                    % far:
                    relativeTolerance = relativeTolerance * rv0(1) / rv0(end);
                end
            end
        end
        
        % If last iterative solver failed, use direct solver.
        if fl0 ~= 0 
            % Last iterative solver failed
            if tryDirectSolverIfIterativeSolversFail
                fprintf('Switching to direct solution since iterative solver(s) failed.\n')
                tic
                soln = matrix \ rhs;
                residual = 0;
                didItConverge = true;
                fprintf('Time to solve system: %g seconds.\n',toc)
                break
            else
                % Accept the most recent solution, even though it doesn't
                % meet the desired tolerance:
                soln(:,col) = soln0;
            end
        end
    end
    
    
end

fprintf('Total time for solve: %g sec.\n',toc(timeWhenSolveStarted))

end