function sfincs_main()

startTime = tic;

fprintf('***************************************************************************\n')
fprintf('SFINCS: Stellarator Fokker-Plank Iterative Neoclassical Conservative Solver\n')
fprintf('Version 3 - Matlab serial edition\n')

sfincs_validateInput()

% Initialize NPeriods, psiAHat, and aHat.  We need to know NPeriods before
% we can initialize the zeta grid.
sfincs_initializeGeometry()

% Do various calculations that will not need to be repeated at each
% iteration, such as setting up the coordinate grids and evaluating
% the magnetic field and its derivatives on the spatial grid.
sfincs_createGrids()

global RHSMode
global nu_n nuPrime B0OverBBar GHat iota IHat
global dPhiHatdpsiHat gamma Delta EStar
if RHSMode == 3
    % Monoenergetic coefficient computation.
    % Overwrite nu_n and dPhiHatd* using nuPrime and EStar.
    
    nu_n = nuPrime * B0OverBBar / (GHat + iota * IHat);
    dPhiHatdpsiHat = 2 / (gamma * Delta) * EStar * iota * B0OverBBar / GHat;    
end


% For input quantities that depend on the radial coordinate, pick out the values for the selected
% radial coordinate, and use these values to over-write values for the other radial coordinates.
sfincs_setInputRadialCoordinate()

% Create HDF5 data structures, and save the quantities that will not change
% at each iteration of the solver (i.e. save all quantities except diagnostics.)
%call initializeOutputFile()

% Solve the main system, either linear or nonlinear.
% This step takes more time than everything else combined.
sfincs_solver()

%sfincs_finalizeHDF5()

fprintf('Total execution time: %g seconds.\n',toc(startTime))
fprintf('Good bye!\n')

end
