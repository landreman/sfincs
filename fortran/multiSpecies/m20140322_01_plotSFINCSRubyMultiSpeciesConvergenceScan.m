function m20140322_01_plotSFINCSRubyMultiSpeciesConvergenceScan()

% Name of .h5 HDF5 file from SFINCS:
h5filename='sfincsOutput.h5';

excludeRunsThatDidntConverge = true;
%excludeRunsThatDidntConverge = false;

figureOffset = 0;

colors = [1,0,0;
    0.8,0.6,0;
    0,0.7,0;
    0,0.8,0.9;
    0,0,1;
    1,0,1;
    0.6,0.6,0.6;
    0,0,0];

linespecs = {'.-r','.-g','.-b','.-m','.-c','.-r','.-r','.-b','.-m'};

numRuns = 0;
RHSMode1s = 1;
RHSMode2s = 0;
files = dir();
dumpedFieldsYet = false;
didItConverges = [];
Nthetas = [];
Nzetas = [];
Nxis = [];
NLs = [];
Nxs = [];
NxPotentialsPerVths = [];
xMaxs = [];
log10tols = [];
outputs = [];
Nspecies = -1;

for iFile = 1:size(files,1)
    if ~ files(iFile).isdir
        continue
    end

    % Skip the . and .. directories
    if strcmp(files(iFile).name,'.') | strcmp(files(iFile).name,'..')
        continue
    end

    % Try to open an HDF5 file
    filename = [files(iFile).name, '/', h5filename];
    try
        info = h5info(filename);
        fprintf('Successfully opened h5 file %s\n',filename)
    catch
        fprintf('Did not succeed in opening h5 file %s\n',filename)
        continue
    end

    if ~ dumpedFieldsYet
        fprintf('Fields saved in the HDF5 files:\n')
        for i=1:numel(info.Groups(1).Datasets)
            fprintf('  %s\n',info.Groups(1).Datasets(i).Name)
        end
        dumpedFieldsYet = true;
    end

    programMode = h5read(filename,'/programMode');
    if programMode ~= 1
        fprintf('Ignoring this run since programMode is not 1.\n')
        continue
    end

    location  = getLocationString(1);
    integerToRepresentTrue = h5read(filename,[location,'integerToRepresentTrue']);
    didItConverge = h5read(filename,[location,'didItConverge']);
    if excludeRunsThatDidntConverge && (didItConverge ~= integerToRepresentTrue)
        fprintf('Ignoring this run since it did not converge.\n')
        continue
    end


    % If we made it this far, then let's count the run.
    numRuns = numRuns + 1;
    didItConverges(numRuns) = (didItConverge == integerToRepresentTrue);
    if didItConverge ~= integerToRepresentTrue
        beep
        fprintf('Warning: run %d did not converge.\n',i)
    end

    %{
    RHSMode = h5read(filename,[location,'RHSMode']);
    switch RHSMode
        case 1
            RHSMode1s = RHSMode1s + 1;
        case 2
            RHSMode2s = RHSMode2s + 1;            
        otherwise
            error('Unrecognized RHSMode')
    end
    if RHSMode1s > 0 && RHSMode2s > 0
        error('Runs must either all have RHSMode=1 or RHSMode=2')
    end
    %}

    Nspecies_new = h5read(filename,[location,'Nspecies']);
    if Nspecies < 0
        Nspecies = Nspecies_new;
    else
        if Nspecies ~= Nspecies_new
            error('Number of species is not consistent among runs')
        end
    end

    runNum = numRuns;
    Nthetas(runNum) = h5read(filename,[location,'Ntheta']);
    Nzetas(runNum) = h5read(filename,[location,'Nzeta']);
    Nxis(runNum) = h5read(filename,[location,'Nxi']);
    NLs(runNum) = h5read(filename,[location,'NL']);
    Nxs(runNum) = h5read(filename,[location,'Nx']);
    NxPotentialsPerVths(runNum) = h5read(filename,[location,'NxPotentialsPerVth']);
    xMaxs(runNum) = h5read(filename,[location,'xMax']);
    log10tols(runNum) = -log10(h5read(filename,[location,'solverTolerance']));
    fprintf('%d %d %d %d %d %g %g %g\n',Nthetas(runNum), ...
            Nzetas(runNum),Nxis(runNum),NLs(runNum),Nxs(runNum), ...
            NxPotentialsPerVths(runNum), xMaxs(runNum), log10tols(runNum))

    outputs(runNum,((1:Nspecies)-1)*3+1) = h5read(filename,[location,'particleFlux']);
    outputs(runNum,((1:Nspecies)-1)*3+2) = h5read(filename,[location,'heatFlux']);
    outputs(runNum,((1:Nspecies)-1)*3+3) = h5read(filename,[location,'FSABFlow']);
    if Nspecies == 1
        outputs(runNum,4) = h5read(filename,[location,'didItConverge']);
        outputs(runNum,5) = h5read(filename,[location,'elapsed time (s)']);
    end
    
    %{
    if RHSMode1s > 0
        outputs(runNum,1) = h5read(filename,[location,'particleFlux']);
        outputs(runNum,2) = h5read(filename,[location,'heatFlux']);
        outputs(runNum,3) = h5read(filename,[location,'FSABFlow']);
        outputs(runNum,4) = h5read(filename,[location,'didItConverge']);
        outputs(runNum,5) = h5read(filename,[location,'elapsed time (s)']);
    else
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
    %}
    
end

if Nspecies == 1
    yAxesLabels = {'Particle flux','q','<V|| B>','Did it converge','elapsed time'};
    numQuantities = numel(yAxesLabels);
    plotRows = 1:numQuantities;
    numRows=numQuantities;
else
    yAxesLabels=cell(0);
    for i=1:Nspecies
        yAxesLabels{end+1} = ['Particle flux, species ', num2str(i)];
        yAxesLabels{end+1} = ['Heat flux, species ', num2str(i)];
        yAxesLabels{end+1} = ['<V|| B>, species ', num2str(i)];
    end
    numQuantities = numel(yAxesLabels);
    plotRows = 1:numQuantities;
    numRows=numQuantities;
end

%{
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
%}

parametersToVary = {};
abscissae = {};
convergeds = {};
quantities = {};

%elapsedTimes = zeros(numRuns,1);
%didItConverges = zeros(numRuns,1);

% Check whether Ntheta was scanned:
data=Nthetas;
[values, runIndices, scanIndices] = unique(data,'first');
if numel(values)>1
    parametersToVary{end+1} = 'Ntheta';
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
    parametersToVary{end+1} = 'Nzeta';
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
    parametersToVary{end+1} = 'Nxi';
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
clf
set(gcf,'Color','w')

numCols = numParameters;

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

temp=dbstack;
nameOfThisProgram=sprintf('%s',temp(1).file);
stringForTop = ['Convergence scan from fortran multi-species version of SFINCS, plotted using ',nameOfThisProgram];


annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
    'Interpreter','none','VerticalAlignment','bottom',...
    'FontSize',12,'LineStyle','none','String',stringForTop);

stringForBottom = ['Run in: ',pwd];

annotation('textbox',[0 0 1 .04],'HorizontalAlignment','center',...
           'Interpreter','none','VerticalAlignment','top',...
           'FontSize',12,'LineStyle','none','String', ...
           stringForBottom);

figure(2+figureOffset)
clf
set(gcf,'Color','w')

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

annotation('textbox',[0 0 1 .04],'HorizontalAlignment','center',...
           'Interpreter','none','VerticalAlignment','top',...
           'FontSize',12,'LineStyle','none','String', ...
           stringForBottom);


% --------------------------------------------------------

    function l = getLocationString(runNum)
        l = sprintf('/run%3d/',runNum);
    end

end