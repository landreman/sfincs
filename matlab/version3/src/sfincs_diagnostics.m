function sfincs_diagnostics()

global stateVector f0
global Ntheta Nzeta Nspecies Nx Nxi
global Phi1Hat dPhi1Hatdtheta dPhi1Hatdzeta ddtheta ddzeta includePhi1

if includePhi1
    for itheta = 1:Ntheta
        indices = sfincs_indices(1,1,1,itheta,1:Nzeta,BLOCK_QN);
        Phi1Hat(itheta,:) = stateVector(indices);
    end
    
    dPhi1Hatdtheta = ddtheta * Phi1Hat;
    dPhi1Hatdzeta = (ddzeta * (Phi1Hat'))';
end

fullf = stateVector + f0;

end