function sfincs_plotB()

global iota BDotCurlB
global BHat dBHatdtheta dBHatdzeta dBHatdpsiHat
global DHat BHat_sub_psi BHat_sub_theta BHat_sub_zeta BHat_sup_theta BHat_sup_zeta
global dBHat_sub_psi_dtheta dBHat_sub_psi_dzeta
global dBHat_sub_theta_dpsiHat dBHat_sub_theta_dzeta
global dBHat_sub_zeta_dpsiHat dBHat_sub_zeta_dtheta
global dBHat_sup_theta_dpsiHat dBHat_sup_theta_dzeta
global dBHat_sup_zeta_dpsiHat dBHat_sup_zeta_dtheta
global theta2D zeta2D Nzeta zeta

figureOffset = 0;
numContours = 20;

if Nzeta==1
    % This is a tokamak
    
else
    % Nzeta > 1, so this is a stellarator.
    
    % ********************************************************************
    
    figure(1+figureOffset)
    clf
    numRows = 2;
    numCols = 4;
    plotNum = 1;
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, BHat, numContours, 'EdgeColor','none')
    colorbar
    hold on
    if iota>0
        plot([0, max(zeta)], [0, max(zeta)*iota], 'k')
    else
        plot([0, max(zeta)], [-max(zeta)*iota, 0], 'k')
    end
    xlabel('\theta')
    ylabel('\zeta')
    title('BHat with a field line','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, BHat, numContours, 'EdgeColor','none')
    colorbar
    hold on
    plot(zeta2D, theta2D, '.k')
    xlabel('\theta')
    ylabel('\zeta')
    title('BHat. Dots show (theta,zeta) grid.','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, BHat_sup_theta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('BHat_sup_theta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, BHat_sup_zeta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('BHat_sup_zeta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, DHat, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('DHat (inverse Jacobian)','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, BHat_sub_psi, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('BHat_sub_psi','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, BHat_sub_theta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('BHat_sub_theta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, BHat_sub_zeta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('BHat_sub_zeta','Interpreter','none')
        
    % ********************************************************************
    
    figure(2+figureOffset)
    clf
    numRows = 2;
    numCols = 3;
    plotNum = 1;
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, dBHat_sub_psi_dtheta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dBHat_sub_psi_dtheta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, dBHat_sub_psi_dzeta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dBHat_sub_psi_dzeta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, dBHat_sub_theta_dpsiHat, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dBHat_sub_theta_dpsiHat','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, dBHat_sub_theta_dzeta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dBHat_sub_theta_dzeta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, dBHat_sub_zeta_dpsiHat, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dBHat_sub_zeta_dpsiHat','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, dBHat_sub_zeta_dtheta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dBHat_sub_zeta_dtheta','Interpreter','none')
    
    % ********************************************************************
    
    figure(3+figureOffset)
    clf
    numRows = 2;
    numCols = 4;
    plotNum = 1;
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, dBHatdpsiHat, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dBHatdpsiHat','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, dBHatdtheta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dBHatdtheta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, dBHatdzeta, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dBHatdzeta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    mu0 = 1.25663706143592e-06;
    data = (BHat_sup_theta.*(dBHat_sub_psi_dtheta-dBHat_sub_theta_dpsiHat)+BHat_sup_zeta.*(dBHat_sub_psi_dzeta-dBHat_sub_zeta_dpsiHat))/mu0;
    contourf(zeta2D, theta2D, data, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('dp/dpsi computed from curl(B) cross B','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    data = DHat.*(dBHat_sub_zeta_dtheta-dBHat_sub_theta_dzeta);
    contourf(zeta2D, theta2D, data, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('curl(B) dot grad psiHat','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    data = DHat.*(dBHat_sub_psi_dzeta-dBHat_sub_zeta_dpsiHat);
    contourf(zeta2D, theta2D, data, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('curl(B) dot grad theta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    data = DHat.*(dBHat_sub_theta_dpsiHat-dBHat_sub_psi_dtheta);
    contourf(zeta2D, theta2D, data, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('curl(B) dot grad zeta','Interpreter','none')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    contourf(zeta2D, theta2D, BDotCurlB, numContours, 'EdgeColor','none')
    colorbar
    xlabel('\theta')
    ylabel('\zeta')
    title('BDotCurlB','Interpreter','none')
    

end

end