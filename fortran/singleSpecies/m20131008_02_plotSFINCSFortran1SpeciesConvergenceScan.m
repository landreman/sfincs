function m20131008_02_plotSFINCSFortran1SpeciesConvergenceScan()

% Name of .h5 HDF5 file from SFINCS:
filename='sfincsOutput.h5';
%filename='C:\Users\landreman\Documents\Ubuntu\sfincsOutput_20131008-04_convergenceScan.h5';

%excludeRunsThatDidntConverge = true;
excludeRunsThatDidntConverge = false;

figureOffset = 0;

colors = [1,0,0;
    0.8,0.6,0;
    0,0.7,0;
    0,0.8,0.9;
    0,0,1;
    1,0,1;
    0.6,0.6,0.6;
    0,0,0];


figureOffset = 50;

info = h5info(filename);
fprintf('Fields saved in this HDF5 file:\n')
for i=1:numel(info.Groups(1).Datasets)
    fprintf('  %s\n',info.Groups(1).Datasets(i).Name)
end

numRuns = numel(info.Groups);

programMode = h5read(filename,'/programMode');

switch programMode
    case 1
        error('This HDF5 file corresponds to a single run (programMode=1), not a convergence scan (programMode=2).')
    case {2}
    otherwise
        error('This HDF5 file does not correspond to a convergence scan (programMode=2).')
end

location  = getLocationString(1);
integerToRepresentTrue = h5read(filename,[location,'integerToRepresentTrue']);
    
RHSMode1s = 0;
RHSMode2s = 0;
didItConverges = false(numRuns,1);
for i=1:numRuns
    location  = getLocationString(i);
    didItConverge = h5read(filename,[location,'didItConverge']);
    didItConverges(i) = (didItConverge == integerToRepresentTrue);
    RHSMode = h5read(filename,[location,'RHSMode']);
    switch RHSMode
        case 1
            RHSMode1s = RHSMode1s + 1;
        case 2
            RHSMode2s = RHSMode2s + 1;            
        otherwise
            error('Unrecognized RHSMode')
    end
    if didItConverge ~= integerToRepresentTrue
        beep
        fprintf('Warning: run %d did not converge.\n',i)
    end
end

if RHSMode1s > 0 && RHSMode2s > 0
    error('Either all runs should have RHSMode1 or all runs should have RHSMode2')
end

if RHSMode1s > 0
    % Scan used RHSMode = 1
    yAxesLabels = {'Particle flux','q','<V|| B>','Did it converge','elapsed time'};
    numQuantities = numel(yAxesLabels);
    plotRows = 1:numQuantities;
    numRows=numQuantities;
else
    % Scan used RHSMode = 2
    yAxesLabels = {'L_{11}','L_{12}=L_{21}','L_{13}=L_{31}','L_{22}','L_{23}=L_{32}','L_{33}'};
    plotRows = [1, 2, 3, 2, 4, 5, 3, 5, 6];
    numQuantities = 9;
    numRows = 6;
end

runsToKeep = 1:numRuns;
if excludeRunsThatDidntConverge
    runsToKeep(~ didItConverges) = [];
    numRuns = numel(runsToKeep);
end

linespecs = {'.-r','.-g','.-b','.-m','.-c','.-r','.-r','.-b','.-m'};
%numQuantities = numel(quantitiesToRecord);

parametersToVary = {};
abscissae = {};
convergeds = {};
quantities = {};

Nthetas = zeros(numRuns,1);
Nzetas = zeros(numRuns,1);
Nxis = zeros(numRuns,1);
NLs = zeros(numRuns,1);
Nxs = zeros(numRuns,1);
NxPotentialsPerVths = zeros(numRuns,1);
xMaxs = zeros(numRuns,1);
log10tols = zeros(numRuns,1);

outputs = zeros(numRuns, numQuantities);
%elapsedTimes = zeros(numRuns,1);
%didItConverges = zeros(numRuns,1);

for runNum = 1:numRuns
    location  = getLocationString(runsToKeep(runNum));
    Nthetas(runNum) = h5read(filename,[location,'Ntheta']);
    Nzetas(runNum) = h5read(filename,[location,'Nzeta']);
    Nxis(runNum) = h5read(filename,[location,'Nxi']);
    NLs(runNum) = h5read(filename,[location,'NL']);
    Nxs(runNum) = h5read(filename,[location,'Nx']);
    NxPotentialsPerVths(runNum) = h5read(filename,[location,'NxPotentialsPerVth']);
    xMaxs(runNum) = h5read(filename,[location,'xMax']);
    log10tols(runNum) = -log10(h5read(filename,[location,'solverTolerance']));
 
    switch RHSMode
        case 1
            outputs(runNum,1) = h5read(filename,[location,'particleFlux']);
            outputs(runNum,2) = h5read(filename,[location,'heatFlux']);
            outputs(runNum,3) = h5read(filename,[location,'FSAFlow']);
            outputs(runNum,4) = h5read(filename,[location,'didItConverge']);
            outputs(runNum,5) = h5read(filename,[location,'elapsed time (s)']);
        case 2
            transportMatrix = h5read(filename,[location,'transportMatrix']);
            outputs(runNum,1) = transportMatrix(1,1);
            outputs(runNum,2) = transportMatrix(1,2);
            outputs(runNum,3) = transportMatrix(1,3);
            outputs(runNum,4) = transportMatrix(2,1);
            outputs(runNum,5) = transportMatrix(2,2);
            outputs(runNum,6) = transportMatrix(2,3);
            outputs(runNum,7) = transportMatrix(3,1);
            outputs(runNum,8) = transportMatrix(3,2);
            outputs(runNum,9) = transportMatrix(3,3);
    end
end

% Check whether Ntheta was scanned:
data=Nthetas;
[values, runIndices, scanIndices] = unique(data,'first');
if numel(values)>1
    parametersToVary{end+1} = 'N\theta';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error('More than one value of input parameter Ntheta was repeated.')
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = outputs(runIndices, :);
end

% Check whether Nzeta was scanned:
data=Nzetas;
[values, runIndices, scanIndices] = unique(data,'first');
if numel(values)>1
    parametersToVary{end+1} = 'N\zeta';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error('More than one value of input parameter Nzeta was repeated.')
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = outputs(runIndices, :);
end

% Check whether Nxi was scanned:
data = Nxis;
[values, runIndices, scanIndices] = unique(data,'first');
if numel(values)>1
    parametersToVary{end+1} = 'N\xi';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error('More than one value of input parameter Nxi was repeated.')
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = outputs(runIndices, :);
end

% Check whether NL was scanned:
data = NLs;
[values, runIndices, scanIndices] = unique(data,'first');
if numel(values)>1
    parametersToVary{end+1} = 'NL';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error('More than one value of input parameter NL was repeated.')
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = outputs(runIndices, :);
end

% Check whether Nx was scanned:
data = Nxs;
[values, runIndices, scanIndices] = unique(data,'first');
if numel(values)>1
    parametersToVary{end+1} = 'Nx';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error('More than one value of input parameter Nx was repeated.')
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = outputs(runIndices, :);
end

% Check whether NxPotentialsPerVth was scanned:
data = NxPotentialsPerVths;
[values, runIndices, scanIndices] = unique(data,'first');
if numel(values)>1
    parametersToVary{end+1} = 'NxPotentialsPerVth';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error('More than one value of input parameter NxPotentialsPerVths was repeated.')
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = outputs(runIndices, :);
end

% Check whether xMax was scanned:
data = xMaxs;
[values, runIndices, scanIndices] = unique(data,'first');
if numel(values)>1
    parametersToVary{end+1} = 'xMax';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error('More than one value of input parameter xMax was repeated.')
    end
    convergeds{end+1} = convergedValue;
    quantities{end+1} = outputs(runIndices, :);
end


% Check whether log10tol was scanned:
data = log10tols;
[values, runIndices, scanIndices] = unique(data,'first');
if numel(values)>1
    parametersToVary{end+1} = 'log_{10}tol';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error('More than one value of input parameter log10tol was repeated.')
    end
    convergeds{end+1} = convergedValue;
    quantities{end+1} = outputs(runIndices, :);
end


numParameters = numel(parametersToVary);


%maxs=ones(numQuantities,1)*(-1e10);
%mins=ones(numQuantities,1)*(1e10);
maxs=ones(numQuantities,1)*(-Inf);
mins=ones(numQuantities,1)*Inf;
for iParameter = 1:numParameters
    maxs = max([maxs, quantities{iParameter}'],[],2);
    mins = min([mins, quantities{iParameter}'],[],2);
end



figure(1+figureOffset)
numCols = numParameters;
clf
for iQuantity = 1:numQuantities
    if maxs(iQuantity) <= mins(iQuantity)
        maxs(iQuantity) = mins(iQuantity)+1;
    end
    for iParameter = 1:numParameters
        subplot(numRows, numCols, iParameter  + (plotRows(iQuantity) - 1)*numParameters)
        plot(1./abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', linespecs{iQuantity})
        hold on
        plot(1./[convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
        ylim([mins(iQuantity), maxs(iQuantity)])
        xlabel(['1/',parametersToVary{iParameter}])
        ylabel(yAxesLabels{plotRows(iQuantity)})
    end
end
stringForTop=sprintf('Convergence scan from fortran version of SFINCS (1 species)');

annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
    'Interpreter','none','VerticalAlignment','bottom',...
    'FontSize',12,'LineStyle','none','String',stringForTop);

figure(2+figureOffset)
clf
for iQuantity = 1:numQuantities
    if maxs(iQuantity) <= mins(iQuantity)
        maxs(iQuantity) = mins(iQuantity)+1;
    end
    for iParameter = 1:numParameters
        subplot(numRows, numCols, iParameter  + (plotRows(iQuantity) - 1)*numParameters)
        plot(abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', linespecs{iQuantity})
        hold on
        plot([convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
        ylim([mins(iQuantity), maxs(iQuantity)])
        xlabel(parametersToVary{iParameter})
        ylabel(yAxesLabels{plotRows(iQuantity)})
    end
end

annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
    'Interpreter','none','VerticalAlignment','bottom',...
    'FontSize',12,'LineStyle','none','String',stringForTop);

% --------------------------------------------------------

    function l = getLocationString(runNum)
        l = sprintf('/run%3d/',runNum);
    end

end