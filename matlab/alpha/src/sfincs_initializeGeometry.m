function sfincs_initializeGeometry()

% In this function, we initialize NPeriods, psiAHat, and aHat.  We need to know NPeriods before
% we can initialize the zeta grid.

global geometryScheme NPeriods helicity_n equilibriumFile aHat psiAHat rN_wish inputRadialCoordinate vmec

switch geometryScheme
    case 1
        % 3-helicity model.
        NPeriods = max([1, helicity_n]);
        % aHat and psiAHat are specified by the user.
        
    case 2
        % LHD standard configuration.
        % Values taken from Table 1 of
        % Beidler et al, Nuclear Fusion 51, 076001 (2011).
        NPeriods = 10;
        B0OverBBar = 1; % (Tesla)
        aHat = 0.5585; % (meters)
        psiAHat = B0OverBBar*aHat^2/2;
        rN_wish = 0.5;
        inputRadialCoordinate = 3;
        
    case 3
        % LHD inward-shifted configuration.
        % Values taken from Table 1 of
        % Beidler et al, Nuclear Fusion 51, 076001 (2011).
        NPeriods = 10;
        B0OverBBar = 1; % (Tesla)
        aHat = 0.5400; % (meters)
        psiAHat = B0OverBBar*aHat^2/2;
        rN_wish = 0.5;
        inputRadialCoordinate = 3;

    case 4
        % W7-X Standard configuration
        % Values taken from Table 1 of
        % Beidler et al, Nuclear Fusion 51, 076001 (2011).
        NPeriods = 5;
        aHat = 0.5109; % (meters)
        psiAHat = -0.384935;
        rN_wish = 0.5;
        inputRadialCoordinate = 3;

    case 5
        % Read VMEC netCDF file.
        sfincs_readVmec();
        
        NPeriods = vmec.nfp;
        psiAHat = vmec.phi(end)/(2*pi);
        aHat = vmec.Aminor_p;
    
    case 11
        fid = fopen(equilibriumFile);
        if fid<0
            error('Unable to open file %s\n',equilibriumFile)
        end
        try
            tmp_str=fgetl(fid);       %Skip comment line
            while strcmp(tmp_str(1:2),'CC');
                tmp_str=fgetl(fid);     %Skip comment line
            end
            header=fscanf(fid,'%d %d %d %d %f %f %f\n',7);
            NPeriods = header(4);
            psiAHat = header(5)/(2*pi);
            aHat = header(6);
            fclose(fid);
        catch me
            error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec .bc output file.\n',...
                me.message, equilibriumFile)
        end
        
    case 12
        fid = fopen(equilibriumFile);
        if fid<0
            error('Unable to open file %s\n',equilibriumFile)
        end
        try
            tmp_str=fgetl(fid);       %Skip comment line
            while strcmp(tmp_str(1:2),'CC');
                tmp_str=fgetl(fid);     %Skip comment line
            end
            header=fscanf(fid,'%d %d %d %d %f %f %f\n',7);
            NPeriods = header(4);
            psiAHat = header(5)/(2*pi);
            aHat = header(6);
            fclose(fid);
        catch me
            error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec .bc output file.\n',...
                me.message, equilibriumFile)
        end
    otherwise
        error('Invalid setting for geometryScheme')
end

end