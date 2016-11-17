function sfincs_computeBHat()

global rN geometryScheme plotB
global BDotCurlB DHat BHat_sub_theta BHat_sub_zeta BHat_sub_psi
global dBHat_sub_psi_dzeta dBHat_sub_zeta_dpsiHat dBHat_sub_theta_dpsiHat dBHat_sub_psi_dtheta dBHat_sub_zeta_dtheta dBHat_sub_theta_dzeta
global force0RadialCurrentInEquilibrium


% Using the selected radial coordinate, set input quantities for the other radial coordinates:
sfincs_setInputRadialCoordinateWish()
% Note that this call only sets the "wish" radial coordinates, not the final radial coordinates
% or the input gradients. These quantities will be set later as we load the magnetic
% geometry, in case the final radial coordinate is different from the "wish" values.

% Set the radius to a silly value here to make sure the proper value is set eventually:
rN = -9999;

switch geometryScheme
    case {1,2,3,4,11,12}
        computeBHat_Boozer()
    case 5
        computeBHat_VMEC()
    otherwise
        error('Invalid geometryScheme')
end

BDotCurlB = DHat .* (...
    BHat_sub_theta .* (dBHat_sub_psi_dzeta - dBHat_sub_zeta_dpsiHat) ...
    + BHat_sub_zeta .* (dBHat_sub_theta_dpsiHat - dBHat_sub_psi_dtheta));

if ~force0RadialCurrentInEquilibrium
    BDotCurlB = BDotCurlB + DHat .* BHat_sub_psi .* (dBHat_sub_zeta_dtheta - dBHat_sub_theta_dzeta);
end

if plotB
    sfincs_plotB()
end

end

% ***********************************************************
% ***********************************************************

function computeBHat_Boozer()

global B0OverBBar GHat IHat iota epsilon_t epsilon_h epsilon_antisymm
global helicity_l helicity_n helicity_antisymm_l helicity_antisymm_n
global psiAHat Nalpha Nzeta NPeriods rN rN_wish
global equilibriumFile min_Bmn_to_load
global BHat dBHatdtheta dBHatdzeta dBHatdpsiHat
global DHat BHat_sub_psi BHat_sub_theta BHat_sub_zeta BHat_sup_theta BHat_sup_zeta
global dBHat_sub_psi_dtheta dBHat_sub_psi_dzeta
global dBHat_sub_theta_dpsiHat dBHat_sub_theta_dzeta
global dBHat_sub_zeta_dpsiHat dBHat_sub_zeta_dtheta
global dBHat_sup_theta_dpsiHat dBHat_sup_theta_dzeta
global dBHat_sup_zeta_dpsiHat dBHat_sup_zeta_dtheta
global alpha2D theta2D zeta2D geometryScheme

switch geometryScheme
    case 1
        % 2-helicity model:
        BHarmonics_l = [1, helicity_l, helicity_antisymm_l];
        if helicity_n==0
            BHarmonics_n = [0, 0, helicity_antisymm_n];
        else
            BHarmonics_n = [0, 1, helicity_antisymm_n/helicity_n];
        end
        BHarmonics_amplitudes = [epsilon_t, epsilon_h, epsilon_antisymm];
        BHarmonics_parity = [true, true, false];
        
        if helicity_n == 0
            if helicity_antisymm_n ~= 0
                fprintf('WARNING: Typically, helicity_antisymm_n should be an integer multiple of helicity_n (possibly 0).\n')
            end
        else
            if mod(helicity_antisymm_n, helicity_n) ~= 0
                fprintf('WARNING: Typically, helicity_antisymm_n should be an integer multiple of helicity_n (possibly 0).\n')
            end
        end
        rN = rN_wish;
        
    case 2
        % LHD standard configuration.
        % Values taken from Table 1 of
        % Beidler et al, Nuclear Fusion 51, 076001 (2011).
        iota = 0.4542;
        BHarmonics_l = [1, 2, 1];
        BHarmonics_n = [0, 1, 1];
        BHarmonics_amplitudes = [-0.07053, 0.05067, -0.01476];
        BHarmonics_parity = [1, 1, 1];
        
        B0OverBBar = 1; % (Tesla)
        R0 = 3.7481; % (meters)
        GHat = B0OverBBar * R0;
        IHat = 0;
        dGdpHat=NaN;
        rN = 0.5;
        
    case 3
        % LHD inward-shifted configuration.
        % Values taken from Table 1 of
        % Beidler et al, Nuclear Fusion 51, 076001 (2011).
        iota = 0.4692;
        BHarmonics_l = [1, 2, 1, 0];
        BHarmonics_n = [0, 1, 1, 1];
        BHarmonics_amplitudes = [-0.05927, 0.05267, -0.04956, 0.01045];
        BHarmonics_parity = [1, 1, 1, 1];
        
        B0OverBBar = 1; % (Tesla)
        R0 = 3.6024; % (meters)
        GHat = B0OverBBar * R0;
        IHat = 0;
        dGdpHat=NaN;
        rN = 0.5;
        
    case 4
        % W7-X Standard configuration
        % Values taken from Table 1 of
        % Beidler et al, Nuclear Fusion 51, 076001 (2011).
        iota=0.8700;
        BHarmonics_l = [0, 1, 1];
        BHarmonics_n = [1, 1, 0];
        BHarmonics_amplitudes = [0.04645, -0.04351, -0.01902];
        BHarmonics_parity = [1, 1, 1];
        
        B0OverBBar = 3.089; % (Tesla)
        GHat = -17.885;
        IHat = 0;
        dGdpHat=NaN;
        rN = 0.5;
        
    case 10
        fid = fopen(equilibriumFile);
        % File description:
        % 1st line: 2 integers:     nfp,ns
        % 2nd line: 4 real numbers: aspect,rmax,rmin,betaxis
        % 3rd line: 3 integers:     mboz, nboz, mnboz
        % 4th line: 7 real numbers: iota,pres,beta,phip,phi,bvco,buco
        %
        % Then, you have 'mnboz' lines.
        % If 'mn' is a dummy integer variable that goes from 1 to mnboz,
        % for each value of mn you read
        %
        % m(mn),n(mn),bmn(mn),rmnc(mn),zmns(mn)pmns(m,n),gmn(mn)
        try
            header=fscanf(fid,'%d %d\n %f %f %f %f\n %d %d %d %f %f %f %f %f %f %f',16);
            mnboz=header(9);
            modes =fscanf(fid,'%d %d %g %g %g %g %g',[7,mnboz]);
            fclose(fid);
            
            % scalar values
            %Nper = header(1); %number of field periods
            iota = header(10);
            IHat = header(16);  % Covariant theta comp. of B, known as I in sfincs (in meter * Tesla)
            GHat = header(15);  % Covariant phi comp. of B, known as G in sfincs (in meter * Tesla)
            % Note that the flux at the separatrix is not stored in the
            % file, so we set PsiAHat in the Physics parameters
            % section in the beginning of the program
            
            % mode amplitudes
            if modes(1,1)==0 && modes(2,1)==0
                B0OverBBar=modes(3,1); %The B00 component in Tesla
            else
                error('The first equilibriumFile entry is not the B00 component')
            end
            BHarmonics_l = modes(1,2:end);
            BHarmonics_n = modes(2,2:end) / NPeriods;
            % Make sure all toroidal mode numbers are integers:
            assert(all(BHarmonics_n == round(BHarmonics_n)))
            BHarmonics_amplitudes = modes(3,2:end)/B0OverBBar; % Store the values normalised to the B00 component.
            BHarmonics_parity = ones(1,length(BHarmonics_amplitudes));
            dGdpHat=NaN; %Not implemented yet
        catch me
            error('%s\n\nFile\n\t%s\ndoes not seem to be a valid vmec fort.996 output file.\n',...
                me.message, equilibriumFile)
        end
    case 11
        
        fid = fopen(equilibriumFile);
        if fid<0
            error('Unable to open file %s\n',equilibriumFile)
        end
        fprintf('Reading magnetic geometry from file %s\n',equilibriumFile)
        
        try
            tmp_str=fgetl(fid);
            while strcmp(tmp_str(1:2),'CC');
                tmp_str=fgetl(fid); %Skip comment line
            end
            line = fgetl(fid);
            header=sscanf(line,'%d %d %d %d %f %f %f\n',7);
            fgetl(fid);
            
            max_no_of_modes=500;
            modesm_new=NaN*zeros(1,max_no_of_modes);
            modesn_new=NaN*zeros(1,max_no_of_modes);
            modesb_new=NaN*zeros(1,max_no_of_modes);
            rN_new=-inf;
            no_of_modes_new=NaN;
            iota_new=NaN;
            G_new=NaN;
            I_new=NaN;
            pPrimeHat_new=NaN;
            end_of_file=0;
            
            while (rN_new<rN_wish) && not(end_of_file)
                rN_old=rN_new;
                no_of_modes_old=no_of_modes_new;
                modesm_old=modesm_new;
                modesn_old=modesn_new;
                modesb_old=modesb_new;
                iota_old=iota_new;
                G_old=G_new;
                I_old=I_new;
                pPrimeHat_old=pPrimeHat_new;
                
                fgetl(fid);
                surfheader=fscanf(fid,'%f %f %f %f %f %f\n',6);
                
                rN_new=sqrt(surfheader(1));
                iota_new=surfheader(2);
                % Note that G and I has a minus sign in the following two lines
                % because Ampere's law comes with a minus sign in the left-handed
                % (r,pol,tor) system.
                G_new=-surfheader(3)*NPeriods/2/pi*(4*pi*1e-7); %Tesla*meter
                I_new=-surfheader(4)/2/pi*(4*pi*1e-7);          %Tesla*meter
                pPrimeHat_new=surfheader(5)*(4*pi*1e-7);       % p=pHat \bar{B}^2 / \mu_0
                
                fgetl(fid); %Skip units line
                proceed=1;
                modeind=0;
                while proceed
                    tmp_str=fgetl(fid);
                    if length(tmp_str)==1
                        if tmp_str==-1 %End of file has been reached
                            proceed=0;
                            end_of_file=1;
                        end
                    elseif not(isempty(find(tmp_str=='s'))) %Next flux surface has been reached
                        proceed=0;
                    else
                        tmp=sscanf(tmp_str,'%d %d %f %f %f %f',6);
                        if abs(tmp(6))>min_Bmn_to_load
                            modeind=modeind+1;
                            %if modeind > max_no_of_modes %Unnecessary to check this in matlab
                            %  error(' modeind > max_no_of_modes !')
                            %end
                            modesm_new(modeind)=tmp(1);
                            modesn_new(modeind)=tmp(2);
                            modesb_new(modeind)=tmp(6);
                        end
                    end
                end
                no_of_modes_new=modeind;
                modesm_new(no_of_modes_new+1:end)=NaN;
                modesn_new(no_of_modes_new+1:end)=NaN;
                modesb_new(no_of_modes_new+1:end)=NaN;
            end
            fclose(fid);
        catch me
            error('%s\n\nFile\n\t%s\ndoes not seem to be a valid .bc geometry file.\n',...
                me.message, equilibriumFile)
        end
        
        [~,minind]=min([(rN_old-rN_wish)^2,...
            (rN_new-rN_wish)^2]);
        if minind==1
            BHarmonics_l = modesm_old(1:no_of_modes_old);
            BHarmonics_n = modesn_old(1:no_of_modes_old);
            BHarmonics_amplitudes = modesb_old(1:no_of_modes_old);
            iota=iota_old;
            GHat=G_old;
            IHat=I_old;
            pPrimeHat=pPrimeHat_old;
            rN=rN_old;
        else %minind=2
            BHarmonics_l = modesm_new(1:no_of_modes_new);
            BHarmonics_n = modesn_new(1:no_of_modes_new);
            BHarmonics_amplitudes = modesb_new(1:no_of_modes_new);
            iota=iota_new;
            GHat=G_new;
            IHat=I_new;
            pPrimeHat=pPrimeHat_new;
            rN=rN_new;
        end
        dGdpHat=(G_new-G_old)/(rN_new^2-rN_old^2)/pPrimeHat;
        
        m0inds=find(BHarmonics_l==0);
        n0m0inds=find(BHarmonics_n(m0inds)==0);
        if isempty(n0m0inds)
            error(' B00 component is missing!')
        end
        nm00ind=m0inds(n0m0inds);
        B0OverBBar=BHarmonics_amplitudes(nm00ind); % Assumes \bar{B} = 1 Tesla
        BHarmonics_amplitudes=[BHarmonics_amplitudes(1:nm00ind-1),...
            BHarmonics_amplitudes(nm00ind+1:end)]...
            /B0OverBBar;
        BHarmonics_l = [BHarmonics_l(1:nm00ind-1),...
            BHarmonics_l(nm00ind+1:end)];
        BHarmonics_n = [BHarmonics_n(1:nm00ind-1),...
            BHarmonics_n(nm00ind+1:end)];
        BHarmonics_parity = ones(1,length(BHarmonics_amplitudes));
        
        % Sign correction for files from Joachim Geiger
        if GHat*psiAHat<0
            disp(['This is a stellarator symmetric file from Joachim Geiger.'...
                ' It will now be turned 180 degrees around a ' ...
                'horizontal axis <=> flip the sign of G and I, so that it matches the sign ' ...
                'of its total toroidal flux.'])
            GHat = -GHat;
            IHat = -IHat;
            dGdpHat=-dGdpHat;
        end
        
        %Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
        psiAHat=psiAHat*(-1);           %toroidal direction switch sign
        GHat = GHat*(-1);               %toroidal direction switch sign
        iota = iota*(-1);               %toroidal direction switch sign
        BHarmonics_n=BHarmonics_n*(-1); %toroidal direction switch sign
        
    case 12
        %Non-stellarator symmetric case
        fid = fopen(equilibriumFile);
        if fid<0
            error('Unable to open file %s\n',equilibriumFile)
        end
        
        try
            tmp_str=fgetl(fid);
            while strcmp(tmp_str(1:2),'CC');
                tmp_str=fgetl(fid); %Skip comment line
            end
            header=fscanf(fid,'%d %d %d %d %f %f %f\n',7);
            fgetl(fid);  %Skip variable name line
            
            NPeriods = header(4);
            psiAHat  = header(5)/2/pi; %Convert the flux from Tm^2 to Tm^2/rad
            a        = header(6);      %minor radius %m
            
            max_no_of_modes=1000;
            modesm_new=NaN*zeros(1,max_no_of_modes);
            modesn_new=NaN*zeros(1,max_no_of_modes);
            modesb_new=NaN*zeros(1,max_no_of_modes);
            rN_new=-inf;
            no_of_modes_new=NaN;
            iota_new=NaN;
            G_new=NaN;
            I_new=NaN;
            pPrimeHat_new=NaN;
            end_of_file=0;
            
            while (rN_new<rN_wish) && not(end_of_file)
                rN_old=rN_new;
                no_of_modes_old=no_of_modes_new;
                modesm_old=modesm_new;
                modesn_old=modesn_new;
                modesb_old=modesb_new;
                iota_old=iota_new;
                G_old=G_new;
                I_old=I_new;
                pPrimeHat_old=pPrimeHat_new;
                
                fgetl(fid);
                surfheader=fscanf(fid,'%f %f %f %f %f %f\n',6);
                
                rN_new=sqrt(surfheader(1)); %r/a=sqrt(psi/psi_a)
                iota_new=surfheader(2);
                % Note that G and I has a minus sign in the following two lines
                % because Ampere's law comes with a minus sign in the left-handed
                % (r,pol,tor) system.
                G_new=-surfheader(3)*NPeriods/2/pi*(4*pi*1e-7); %Tesla*meter
                I_new=-surfheader(4)/2/pi*(4*pi*1e-7);          %Tesla*meter
                pPrimeHat_new=surfheader(5)*(4*pi*1e-7);       % p=pHat \bar{B}^2 / \mu_0
                
                fgetl(fid); %Skip units line
                proceed=1;
                modeind=0;
                while proceed
                    tmp_str=fgetl(fid);
                    if length(tmp_str)==1
                        if tmp_str==-1 %End of file has been reached
                            proceed=0;
                            end_of_file=1;
                        end
                    elseif not(isempty(find(tmp_str=='s'))) %Next flux surface has been reached
                        proceed=0;
                    else
                        tmp=sscanf(tmp_str,'%d %d %f %f %f %f %f %f %f %f',10);
                        if (abs(tmp(9))>min_Bmn_to_load) || (abs(tmp(10))>min_Bmn_to_load)
                            modeind=modeind+1;
                            modesm_new(modeind)=tmp(1);
                            modesn_new(modeind)=tmp(2);
                            modesb_new(modeind)=tmp(9); %Cosinus component
                            
                            modeind=modeind+1;
                            modesm_new(modeind)=tmp(1);
                            modesn_new(modeind)=tmp(2);
                            modesb_new(modeind)=tmp(10); %Sinus component
                        end
                    end
                end
                no_of_modes_new=modeind;
                modesm_new(no_of_modes_new+1:end)=NaN;
                modesn_new(no_of_modes_new+1:end)=NaN;
                modesb_new(no_of_modes_new+1:end)=NaN;
            end
            fclose(fid);
        catch me
            error('%s\n\nFile\n\t%s\ndoes not seem to be a valid .bc geometry file.\n',...
                me.message, equilibriumFile)
        end
        
        [~,minind]=min([(rN_old-rN_wish)^2,...
            (rN_new-rN_wish)^2]);
        if minind==1
            BHarmonics_l = modesm_old(1:no_of_modes_old);
            BHarmonics_n = modesn_old(1:no_of_modes_old);
            BHarmonics_amplitudes = modesb_old(1:no_of_modes_old);
            iota=iota_old;
            GHat=G_old;
            IHat=I_old;
            pPrimeHat=pPrimeHat_old;
            rN=rN_old;
        else %minind=2
            BHarmonics_l = modesm_new(1:no_of_modes_new);
            BHarmonics_n = modesn_new(1:no_of_modes_new);
            BHarmonics_amplitudes = modesb_new(1:no_of_modes_new);
            iota=iota_new;
            GHat=G_new;
            IHat=I_new;
            pPrimeHat=pPrimeHat_new;
            rN=rN_new;
        end
        dGdpHat=(G_new-G_old)/(rN_new^2-rN_old^2)/pPrimeHat;
        
        disp(['The calculation is performed for radius ' ...
            ,num2str(rN*a),' m , r/a=',num2str(rN)])
        
        m0inds=find(BHarmonics_l==0);
        n0m0inds=find(BHarmonics_n(m0inds)==0);
        if isempty(n0m0inds)
            error(' B00 component is missing!')
        end
        nm00ind=m0inds(n0m0inds(1));
        B0OverBBar=BHarmonics_amplitudes(nm00ind); %Assumes \bar{B}=1T
        BHarmonics_amplitudes=[BHarmonics_amplitudes(1:nm00ind-1),...
            BHarmonics_amplitudes(nm00ind+2:end)]...
            /B0OverBBar;
        BHarmonics_l = [BHarmonics_l(1:nm00ind-1),...
            BHarmonics_l(nm00ind+2:end)];
        BHarmonics_n = [BHarmonics_n(1:nm00ind-1),...
            BHarmonics_n(nm00ind+2:end)];
        BHarmonics_parity=((-1).^(0:length(BHarmonics_n)-1)+1)/2; %[1,0,1,0,1,0,1,0,...], i.e. cos,sin.cos,sin,...
        
        %Switch from a left-handed to right-handed (radial,poloidal,toroidal) system
        psiAHat=psiAHat*(-1);           %toroidal direction switch sign
        GHat = GHat*(-1);               %toroidal direction switch sign
        iota = iota*(-1);               %toroidal direction switch sign
        BHarmonics_n=BHarmonics_n*(-1); %toroidal direction switch sign
        
    otherwise
        error('Invalid setting for geometryScheme')
end

theta2D = alpha2D + iota*zeta2D;

NHarmonics = numel(BHarmonics_amplitudes);
BHat = B0OverBBar * ones(Nalpha,Nzeta);
dBHatdtheta = zeros(Nalpha,Nzeta);
dBHatdzeta = zeros(Nalpha,Nzeta);
for i=1:NHarmonics
    if BHarmonics_parity(i) %The cosine components of BHat
        BHat = BHat + B0OverBBar * BHarmonics_amplitudes(i) *...
            cos(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
        dBHatdtheta = dBHatdtheta - B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_l(i) *...
            sin(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
        dBHatdzeta = dBHatdzeta + B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_n(i) * NPeriods *...
            sin(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods ...
            * zeta2D);
    else  %The sine components of BHat
        BHat = BHat + B0OverBBar * BHarmonics_amplitudes(i) *...
            sin(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
        dBHatdtheta = dBHatdtheta + B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_l(i) *...
            cos(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods * zeta2D);
        dBHatdzeta = dBHatdzeta - B0OverBBar * BHarmonics_amplitudes(i) * BHarmonics_n(i) * NPeriods *...
            cos(BHarmonics_l(i) * theta2D - BHarmonics_n(i) * NPeriods ...
            * zeta2D);
    end
end
% ---------------------------------------------------------------------------------------
% Calculate parallel current u from harmonics of 1/B^2. Used in NTV calculation.
% \nabla_\parallel u = (2/B^4) \nabla B \times \vector{B} \cdot \iota \nabla \psi
% ---------------------------------------------------------------------------------------
uHat = zeros(Nalpha,Nzeta);
duHatdtheta = zeros(Nalpha,Nzeta);
duHatdzeta = zeros(Nalpha,Nzeta);
hHat=1./(BHat.^2);
if any(BHarmonics_parity==0) %sine components exist
    for m=0:floor(Nalpha/2)-1 %Nyquist max freq.
        if m==0
            nrange=1:floor(Nzeta/2)-1;
        else
            nrange=-floor(Nzeta/2):(floor(Nzeta/2)-1);
        end
        for n=nrange
            %cos
            hHatHarmonics_amplitude = 2/(Nalpha*Nzeta) *...
                sum(sum(cos(m * theta2D  - n * NPeriods * zeta2D).*hHat));
            uHatHarmonics_amplitude = ...
                iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude;
            uHat = uHat + uHatHarmonics_amplitude * cos(m * theta2D - n * NPeriods * zeta2D);
            duHatdtheta = duHatdtheta ...
                - uHatHarmonics_amplitude * m * sin(m * theta2D - n * NPeriods * zeta2D);
            duHatdzeta = duHatdzeta ...
                + uHatHarmonics_amplitude * n * NPeriods * sin(m * theta2D - n * NPeriods * zeta2D);
            
            %sin
            hHatHarmonics_amplitude = 2/(Nalpha*Nzeta) *...
                sum(sum(sin(m * theta2D  - n * NPeriods * zeta2D).*hHat));
            uHatHarmonics_amplitude = ...
                iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude;
            uHat = uHat + uHatHarmonics_amplitude * sin(m * theta2D - n * NPeriods * zeta2D);
            duHatdtheta = duHatdtheta ...
                + uHatHarmonics_amplitude * m * cos(m * theta2D - n * NPeriods * zeta2D);
            duHatdzeta = duHatdzeta ...
                - uHatHarmonics_amplitude * n * NPeriods * cos(m * theta2D - n * NPeriods * zeta2D);
        end
    end
else %only cosinus components
    for m=0:floor(Nalpha/2)-1 %Nyquist max freq.
        if m==0
            nrange=1:floor(Nzeta/2)-1;
        else
            nrange=-floor(Nzeta/2):(floor(Nzeta/2)-1);
        end
        for n=nrange
            hHatHarmonics_amplitude = 2/(Nalpha*Nzeta) *...
                sum(sum(cos(m * theta2D  - n * NPeriods * zeta2D).*hHat));
            uHatHarmonics_amplitude = ...
                iota*(GHat*m + IHat*n * NPeriods)/(n * NPeriods - iota*m) * hHatHarmonics_amplitude;
            uHat = uHat + uHatHarmonics_amplitude * cos(m * theta2D - n * NPeriods * zeta2D);
            duHatdtheta = duHatdtheta ...
                - uHatHarmonics_amplitude * m * sin(m * theta2D - n * NPeriods * zeta2D);
            duHatdzeta = duHatdzeta ...
                + uHatHarmonics_amplitude * n * NPeriods * sin(m * theta2D - n * NPeriods * zeta2D);
        end
    end
end

% This next line needs fixing before the NTV calculation works:
dGdpHat = 0;

NTVkernel = 2/5 * ( ...
    dGdpHat ./ BHat .* (iota * dBHatdtheta + dBHatdzeta) + ...
    1/2 * (iota * (duHatdtheta + uHat * 2./BHat .* dBHatdtheta) ...
    + duHatdzeta + uHat * 2./BHat .* dBHatdzeta) );

% Expressions relating general coordinates and Boozer coordinates:

DHat = BHat.*BHat / (GHat + iota*IHat);
BHat_sup_theta = iota*DHat;
BHat_sup_zeta = DHat;
BHat_sub_theta = IHat*ones(size(BHat));
BHat_sub_zeta = GHat*ones(size(BHat));

% Eventually these next lines could be replaced with a proper calculation:

zero2D = zeros(size(BHat));

dBHatdpsiHat = zero2D;
BHat_sub_psi = zero2D;

dBHat_sub_psi_dtheta = zero2D;
dBHat_sub_psi_dzeta = zero2D;

dBHat_sub_theta_dpsiHat = zero2D;
dBHat_sub_theta_dzeta = zero2D;

dBHat_sub_zeta_dpsiHat = zero2D;
dBHat_sub_zeta_dtheta = zero2D;

dBHat_sup_theta_dpsiHat = zero2D;
dBHat_sup_theta_dzeta = zero2D;

dBHat_sup_zeta_dpsiHat = zero2D;
dBHat_sup_zeta_dtheta = zero2D;


end

% ***********************************************************
% ***********************************************************

function computeBHat_VMEC()

global psiAHat iota psiN_wish rN
global equilibriumFile min_Bmn_to_load VMECRadialOption vmec
global BHat dBHatdtheta dBHatdzeta dBHatdpsiHat
global DHat BHat_sub_psi BHat_sub_theta BHat_sub_zeta BHat_sup_theta BHat_sup_zeta
global dBHat_sub_psi_dtheta dBHat_sub_psi_dzeta
global dBHat_sub_theta_dpsiHat dBHat_sub_theta_dzeta
global dBHat_sub_zeta_dpsiHat dBHat_sub_zeta_dtheta
global dBHat_sup_theta_dpsiHat dBHat_sup_theta_dzeta
global dBHat_sup_zeta_dpsiHat dBHat_sup_zeta_dtheta
global alpha2D theta2D zeta2D

psiN_full = linspace(0,1,vmec.ns);
psiN_half = psiN_full(2:end) - 0.5*(psiN_full(2)-psiN_full(1));

% Options for VMECRadialOption:
% 0: use exact radius requested
% 1: use nearest point in vmec's half mesh
% 2: use nearest point in vmec's full mesh

% First set the actual radius that will be used.
switch VMECRadialOption
    case 0
        psiN = psiN_wish;
    case 1
        errors = abs(psiN_wish - psiN_half);
        [~,index] = min(errors);
        psiN = psiN_half(index);
    case 2
        errors = abs(psiN_wish - psiN_full);
        [~,index] = min(errors);
        psiN = psiN_full(index);
    otherwise
        error('Invalid VMECRadialOption')
end

% Compute vmecRadialIndex_full and vmecRadialWeight_full:
if psiN<0
    error('psiN<0')
elseif psiN>1
    error('psiN>1')
elseif psiN==1
    vmecRadialIndex_full = [vmec.ns-1,vmec.ns];
    vmecRadialWeight_full = [0,1];
else
    % This is the usual case: psiN is >= 0 and < 1.
    index = find(psiN >= psiN_full, 1, 'last');
    if numel(index) ~= 1 || index<1 || index >= vmec.ns
        index
        error('Should not get here')
    end
    vmecRadialIndex_full = [index, index+1];
    weight = (psiN_full(index+1)-psiN) / (psiN_full(index+1)-psiN_full(index));
    vmecRadialWeight_full = [weight, 1-weight];
end

% Compute vmecRadialIndex_half and vmecRadialWeight_half:
if psiN < psiN_half(1)
    fprintf('Warning: extrapolating off the end of the vmec half grid (toward the axis.)\n')
    vmecRadialIndex_half = [2, 3];
elseif psiN > psiN_half(end)
    fprintf('Warning: extrapolating off the end of the vmec half grid (toward the edge.)\n')
    vmecRadialIndex_half = [vmec.ns-1, vmec.ns];
elseif psiN == psiN_half(end)
    vmecRadialIndex_half = [vmec.ns-1, vmec.ns];
else
    %This is the usual case: psiN is in the interior of the half grid, or
    %possibly the first point on the half grid
    index = find(psiN >= psiN_half, 1, 'last');
    if numel(index) ~= 1 || index<1 || index >= vmec.ns-1
        index
        error('Should not get here')
    end
    % In the next line, we add an extra 1 to the index since psiN_half only
    % has vmec.ns-1 elements, and the vmec arrays for half-grid quantities
    % have vmec.ns elements (with the 1st element ignored.)
    vmecRadialIndex_half = [index+1, index+2];
end
weight1 = (psiN_half(vmecRadialIndex_half(2)-1)-psiN) / (psiN_half(vmecRadialIndex_half(2)-1)-psiN_half(vmecRadialIndex_half(1)-1));
vmecRadialWeight_half = [weight1, 1-weight1];

fprintf('vmecRadialIndex_half: %d, %d\n',vmecRadialIndex_half(1), vmecRadialIndex_half(2))
fprintf('vmecRadialWeight_half: %g, %g\n',vmecRadialWeight_half(1), vmecRadialWeight_half(2))
fprintf('vmecRadialIndex_full: %d, %d\n',vmecRadialIndex_full(1), vmecRadialIndex_full(2))
fprintf('vmecRadialWeight_full: %g, %g\n',vmecRadialWeight_full(1), vmecRadialWeight_full(2))

assert(vmec.xn_nyq(1)==0)
assert(vmec.xm_nyq(1)==0)
iota = 0;
B00 = 0;
for i=1:2
    iota = iota + vmecRadialWeight_half(i) * vmec.iotas(vmecRadialIndex_half(i));
    B00 = B00 + vmecRadialWeight_half(i) * vmec.bmnc(1, vmecRadialIndex_half(i));
end

theta2D = alpha2D + iota*zeta2D

BHat = zeros(size(alpha2D));
DHatInverse = zeros(size(alpha2D));
dBHatdpsiHat = zeros(size(alpha2D));
dBHatdtheta = zeros(size(alpha2D));
dBHatdzeta = zeros(size(alpha2D));

BHat_sub_psi = zeros(size(alpha2D));
BHat_sub_theta = zeros(size(alpha2D));
BHat_sub_zeta = zeros(size(alpha2D));
BHat_sup_theta = zeros(size(alpha2D));
BHat_sup_zeta = zeros(size(alpha2D));

dBHat_sub_psi_dtheta = zeros(size(alpha2D));
dBHat_sub_psi_dzeta = zeros(size(alpha2D));

dBHat_sub_theta_dpsiHat = zeros(size(alpha2D));
dBHat_sub_theta_dzeta = zeros(size(alpha2D));

dBHat_sub_zeta_dtheta = zeros(size(alpha2D));
dBHat_sub_zeta_dpsiHat = zeros(size(alpha2D));

dBHat_sup_theta_dpsiHat = zeros(size(alpha2D));
dBHat_sup_theta_dzeta = zeros(size(alpha2D));

dBHat_sup_zeta_dtheta = zeros(size(alpha2D));
dBHat_sup_zeta_dpsiHat = zeros(size(alpha2D));

numSymmetricModesIncluded = 0;
numAntisymmetricModesIncluded = 0;

if vmec.lasym
    fprintf('VMEC file is not stellarator-symmetric.\n')
    error('Not yet implemented')
else
    fprintf('VMEC file is stellarator-symmetric.\n')
end

dpsi = vmec.phi(2)/(2*pi);

% Set to true for exact agreement with fortran version.
% Set to false to include every harmonic in the vmec file.
excludeAliasedFrequencies = true;

for imode = 1:numel(vmec.xn_nyq)
    m = vmec.xm_nyq(imode);
    n = vmec.xn_nyq(imode);
    if excludeAliasedFrequencies && (abs(n)>vmec.ntor*vmec.nfp || m >= vmec.mpol)
        %fprintf('Excluding mode m=%d, n=%d\n',m,n)
        continue
    end
    
    b = vmecRadialWeight_half(1) * vmec.bmnc(imode, vmecRadialIndex_half(1)) ...
        + vmecRadialWeight_half(2) * vmec.bmnc(imode, vmecRadialIndex_half(2));
    
    if abs(b/B00) >= min_Bmn_to_load
        numSymmetricModesIncluded = numSymmetricModesIncluded + 1;

        % Evaluate the radial derivatives we will need.
        
        vmec_dBHatdpsiHat = (vmec.bmnc(imode,3:end) - vmec.bmnc(imode,2:(end-1))) / dpsi;        
        vmec_dBHat_sub_theta_dpsiHat = (vmec.bsubumnc(imode,3:end) - vmec.bsubumnc(imode,2:(end-1))) / dpsi;
        vmec_dBHat_sub_zeta_dpsiHat = (vmec.bsubvmnc(imode,3:end) - vmec.bsubvmnc(imode,2:(end-1))) / dpsi;
        vmec_dBHat_sup_theta_dpsiHat = (vmec.bsupumnc(imode,3:end) - vmec.bsupumnc(imode,2:(end-1))) / dpsi;
        vmec_dBHat_sup_zeta_dpsiHat = (vmec.bsupvmnc(imode,3:end) - vmec.bsupvmnc(imode,2:(end-1))) / dpsi;
        
        % Simplistic extrapolation at the endpoints of the radial grid:
        vmec_dBHatdpsiHat = [vmec_dBHatdpsiHat(1), vmec_dBHatdpsiHat, vmec_dBHatdpsiHat(end)];
        vmec_dBHat_sub_theta_dpsiHat = [vmec_dBHat_sub_theta_dpsiHat(1), vmec_dBHat_sub_theta_dpsiHat, vmec_dBHat_sub_theta_dpsiHat(end)];
        vmec_dBHat_sub_zeta_dpsiHat = [vmec_dBHat_sub_zeta_dpsiHat(1), vmec_dBHat_sub_zeta_dpsiHat, vmec_dBHat_sub_zeta_dpsiHat(end)];
        vmec_dBHat_sup_theta_dpsiHat = [vmec_dBHat_sup_theta_dpsiHat(1), vmec_dBHat_sup_theta_dpsiHat, vmec_dBHat_sup_theta_dpsiHat(end)];
        vmec_dBHat_sup_zeta_dpsiHat = [vmec_dBHat_sup_zeta_dpsiHat(1), vmec_dBHat_sup_zeta_dpsiHat, vmec_dBHat_sup_zeta_dpsiHat(end)];

        % End of radial derivatives.
        
        angle = m * theta2D - n * zeta2D;
        cc = cos(angle);
        ss = sin(angle);
        
        BHat = BHat + b*cc;
        dBHatdtheta = dBHatdtheta - ss*m*b;
        dBHatdzeta = dBHatdzeta - ss*(-n)*b;
        
        for isurf = 1:2
             % Handle Jacobian:                                                                                                                               
             % SFINCS's DHat is the INVERSE Jacobian of the (psiHat, theta, zeta) coordinates.
             % VMEC's gmnc and gmns are the Jacobian of the (psiN, theta, zeta) coordinates.
             % Because one uses psiHat and the other uses psiN, we need a factor of psiAHat for conversion.
             % We will also set DHat = 1 / DHat at the end of this loop.
             DHatInverse = DHatInverse + cc * vmec.gmnc(imode,vmecRadialIndex_half(isurf))*vmecRadialWeight_half(isurf) / psiAHat;
             
             % Handle B sup theta:
             % Note that VMEC's bsupumnc and bsupumns are exactly the same as SFINCS's BHat_sup_theta, with no conversion factors of 2pi needed.
             BHat_sup_theta = BHat_sup_theta + cc * vmec.bsupumnc(imode,vmecRadialIndex_half(isurf))*vmecRadialWeight_half(isurf);
             dBHat_sup_theta_dzeta = dBHat_sup_theta_dzeta - ss*(-n) * vmec.bsupumnc(imode,vmecRadialIndex_half(isurf))*vmecRadialWeight_half(isurf);
             
             % Handle B sup zeta:
             % Note that VMEC's bsupvmnc and bsupvmns are exactly the same as SFINCS's BHat_sup_zeta, with no conversion factors of 2pi or Nperiods needed.
             BHat_sup_zeta = BHat_sup_zeta + cc * vmec.bsupvmnc(imode,vmecRadialIndex_half(isurf))*vmecRadialWeight_half(isurf);
             dBHat_sup_zeta_dtheta = dBHat_sup_zeta_dtheta - ss*m * vmec.bsupvmnc(imode,vmecRadialIndex_half(isurf))*vmecRadialWeight_half(isurf);
             
             % Handle B sub theta:
             % Note that VMEC's bsubumnc and bsubumns are exactly the same as SFINCS's BHat_sub_theta, with no conversion factors of 2pi needed.
             BHat_sub_theta = BHat_sub_theta + cc * vmec.bsubumnc(imode,vmecRadialIndex_half(isurf))*vmecRadialWeight_half(isurf);
             dBHat_sub_theta_dzeta = dBHat_sub_theta_dzeta - ss*(-n) * vmec.bsubumnc(imode,vmecRadialIndex_half(isurf))*vmecRadialWeight_half(isurf);
             
             % Handle B sub zeta:
             % Note that VMEC's bsubvmnc and bsubvmns are exactly the same as SFINCS's BHat_sub_zeta, with no conversion factors of 2pi needed.
             BHat_sub_zeta = BHat_sub_zeta + cc * vmec.bsubvmnc(imode,vmecRadialIndex_half(isurf))*vmecRadialWeight_half(isurf);
             dBHat_sub_zeta_dtheta = dBHat_sub_zeta_dtheta - ss*m * vmec.bsubvmnc(imode,vmecRadialIndex_half(isurf))*vmecRadialWeight_half(isurf);
             
             % Handle B sub psi:
             temp = vmec.bsubsmns(imode, vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf) / psiAHat;
             BHat_sub_psi = BHat_sub_psi + ss * temp;
             dBHat_sub_psi_dtheta = dBHat_sub_psi_dtheta + cc*m * temp;
             dBHat_sub_psi_dzeta = dBHat_sub_psi_dzeta + cc*(-n) * temp;
             
             % Radial derivatives are on the full mesh:
             dBHatdpsiHat = dBHatdpsiHat + cc * vmec_dBHatdpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf);
             dBHat_sub_theta_dpsiHat = dBHat_sub_theta_dpsiHat + cc * vmec_dBHat_sub_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf);
             dBHat_sub_zeta_dpsiHat = dBHat_sub_zeta_dpsiHat + cc * vmec_dBHat_sub_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf);
             dBHat_sup_theta_dpsiHat = dBHat_sup_theta_dpsiHat + cc * vmec_dBHat_sup_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf);
             dBHat_sup_zeta_dpsiHat = dBHat_sup_zeta_dpsiHat + cc * vmec_dBHat_sup_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf);
        end
    end
end

fprintf('%d of %d stellarator-symmetric modes included.\n',numSymmetricModesIncluded, numel(vmec.xn_nyq))
if vmec.lasym
    fprintf('%d of %d stellarator-antisymmetric modes included.\n',numAntisymmetricModesIncluded, numel(vmec.xn_nyq))
end

% Convert Jacobian to inverse Jacobian:
DHat = 1./DHatInverse;

rN = sqrt(psiN);

fprintf('Successfully read VMEC geometry file %s\n',equilibriumFile)

end

