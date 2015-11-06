function sfincs_readVmec()

global vmec equilibriumFile

vmec = struct();
vmec.nfp = double(ncread(equilibriumFile, 'nfp'));
vmec.ns = double(ncread(equilibriumFile, 'ns'));
vmec.ntor = double(ncread(equilibriumFile, 'ntor'));
vmec.mpol = double(ncread(equilibriumFile, 'mpol'));
vmec.phi = ncread(equilibriumFile, 'phi');
vmec.Aminor_p = ncread(equilibriumFile,'Aminor_p');
vmec.iotas = ncread(equilibriumFile,'iotas'); % iota on the half mesh
vmec.bmnc = ncread(equilibriumFile,'bmnc');
vmec.gmnc = ncread(equilibriumFile,'gmnc');
vmec.bsupumnc = ncread(equilibriumFile,'bsupumnc');
vmec.bsupvmnc = ncread(equilibriumFile,'bsupvmnc');
vmec.bsubumnc = ncread(equilibriumFile,'bsubumnc');
vmec.bsubvmnc = ncread(equilibriumFile,'bsubvmnc');
vmec.bsubsmns = ncread(equilibriumFile,'bsubsmns');
vmec.xm_nyq = ncread(equilibriumFile,'xm_nyq');
vmec.xn_nyq = ncread(equilibriumFile,'xn_nyq');
vmec.xn = ncread(equilibriumFile,'xn');
fprintf('VMEC num modes: %d,  num nyq modes: %d\n',numel(vmec.xn), numel(vmec.xn_nyq))

vmec.lasym = (ncread(equilibriumFile, 'lasym__logical__') ~= 0);
if vmec.lasym
    vmec.bmns = ncread(equilibriumFile,'bmns');
    vmec.gmns = ncread(equilibriumFile,'gmns');
    vmec.bsupumns = ncread(equilibriumFile,'bsupumns');
    vmec.bsupvmns = ncread(equilibriumFile,'bsupvmns');
    vmec.bsubumns = ncread(equilibriumFile,'bsubumns');
    vmec.bsubvmns = ncread(equilibriumFile,'bsubvmns');
    vmec.bsubsmnc = ncread(equilibriumFile,'bsubsmnc');
end

end