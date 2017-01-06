function sfincs_setInputRadialCoordinateWish()

global inputRadialCoordinate
global psiHat_wish psiN_wish rHat_wish rN_wish
global psiAHat aHat

% First, use the requeted *_wish coordinate to set psiHat_wish.
switch inputRadialCoordinate
    case 0
        % Nothing to do
        fprintf('Setting radial coordinate based on psiHat_wish.\n')
    case 1
        fprintf('Setting radial coordinate based on psiN_wish.\n')
        assert(psiN_wish >= 0)
        assert(psiN_wish <= 1)
        psiHat_wish = psiN_wish * psiAHat;
    case 2
        fprintf('Setting radial coordinate based on rHat_wish.\n')
        psiHat_wish = (rHat_wish/aHat)^2 * psiAHat;
    case 3
        fprintf('Setting radial coordinate based on rN_wish.\n')
        assert(rN_wish >= 0)
        assert(rN_wish <= 1)
        psiHat_wish = rN_wish^2 * psiAHat;
    otherwise
        error('Invalid inputRadialCoordinate')
end

% Now that psiHat_wish is set, set the other *_wish variables.
psiN_wish = psiHat_wish / psiAHat;
rHat_wish = aHat * sqrt(psiHat_wish / psiAHat);
rN_wish = sqrt(psiHat_wish / psiAHat);

end