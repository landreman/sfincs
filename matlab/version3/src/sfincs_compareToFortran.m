function sfincs_compareToFortran(filename)

global theta zeta x
global geometryScheme GHat IHat VPrimeHat FSABHat2 B0OverBBar iota BDotCurlB
global BHat dBHatdtheta dBHatdzeta dBHatdpsiHat
global DHat BHat_sub_psi BHat_sub_theta BHat_sub_zeta BHat_sup_theta BHat_sup_zeta
global dBHat_sub_psi_dtheta dBHat_sub_psi_dzeta
global dBHat_sub_theta_dpsiHat dBHat_sub_theta_dzeta
global dBHat_sub_zeta_dpsiHat dBHat_sub_zeta_dtheta
global dBHat_sup_theta_dpsiHat dBHat_sup_theta_dzeta
global dBHat_sup_zeta_dpsiHat dBHat_sup_zeta_dtheta
global psiHat psiN rHat rN

%filename = 'sfincsOutput.h5';

fprintf('Comparing matlab results to fortran results in %s\n',filename)

small = 1e-12;

compare('theta')
compare('zeta')
compare('x')
compare('psiHat')
compare('psiN')
compare('rHat')
compare('rN')
compare('geometryScheme')
compare('GHat')
compare('IHat')
compare('B0OverBBar')
compare('VPrimeHat')
compare('FSABHat2')
compare('iota')
compare('BHat')
compare('DHat')
compare('BDotCurlB')
compare('dBHatdpsiHat')
compare('dBHatdtheta')
compare('dBHatdzeta')
compare('BHat_sub_psi')
compare('BHat_sub_theta')
compare('BHat_sub_zeta')
compare('BHat_sup_theta')
compare('BHat_sup_zeta')

compare('dBHat_sub_psi_dtheta')
compare('dBHat_sub_psi_dzeta')

compare('dBHat_sub_theta_dpsiHat')
compare('dBHat_sub_theta_dzeta')

compare('dBHat_sub_zeta_dtheta')
compare('dBHat_sub_zeta_dpsiHat')

compare('dBHat_sup_theta_dpsiHat')
compare('dBHat_sup_theta_dzeta')

compare('dBHat_sup_zeta_dtheta')
compare('dBHat_sup_zeta_dpsiHat')

    function compare(varName)
        matlabVar = eval(varName);
        try
            fortranVar = h5read(filename,['/',varName]);
        catch
            fprintf('** WARNING: Unable to test variable %s since it does not exist in %s\n',varName,fortranVar)
            return
        end
        
        try
            difference = max(max(abs(matlabVar-fortranVar)));
        catch
            fprintf('** ERROR! Matlab and fortran variables %s are different sizes.\n',varName)
            fprintf('Size of matlab version:\n')
            size(matlabVar)            
            fprintf('Size of fortran version:\n')
            size(fortranVar)
            return
        end
        if difference<small
            fprintf('  %s agrees.\n',varName)
        else
            fprintf('** ERROR!  %s disagrees! max(abs(difference)) = %g\n',varName, difference)
            if numel(matlabVar)==1
                fprintf('  Matlab version: %g,  Fortran version: %g\n',matlabVar, fortranVar)
            elseif numel(size(matlabVar)) <= 2 && numel(matlabVar)<100
                fprintf('Here comes matlab version:\n')
                matlabVar
                fprintf('Here comes fortran version:\n')
                fortranVar
            end
            assignin('base',[varName,'_matlab'],matlabVar)
            assignin('base',[varName,'_fortran'],fortranVar)
        end
    end

end