function sfincs_RosenbluthPotentialResponse(polynomials_a, polynomials_b, polynomials_c)

global Nx x xWeights xGrid_k NL Nspecies THats mHats xGridScheme
global Rosenbluth_H Rosenbluth_dHdxb Rosenbluth_d2Gdxb2

a = polynomials_a;
b = polynomials_b;
c = polynomials_c;
integrationPower = 0;

tic
fprintf('Computing Rosenbluth potential response matrices...  ')

collocation2modal = zeros(Nx);
for j = 1:Nx
    for i = 1:Nx
        collocation2modal(j,i) = xWeights(i) * (x(i) .^ xGrid_k) * evaluatePolynomial(x(i))/c(j);
    end
end

tempMatrix_H = zeros(Nx);
tempMatrix_dHdxb = zeros(Nx);
tempMatrix_d2Gdxb2 = zeros(Nx);

Rosenbluth_H = zeros(Nspecies, Nspecies, NL, Nx, Nx);
Rosenbluth_dHdxb = zeros(Nspecies, Nspecies, NL, Nx, Nx);
Rosenbluth_d2Gdxb2 = zeros(Nspecies, Nspecies, NL, Nx, Nx);

if xGridScheme==6
    ixMin = 2;
else
    ixMin = 3;
end

for L = 0:(NL-1)
    alpha = -(2*L-1)/(2*L+3);
    
    for iSpeciesA = 1:Nspecies
        for iSpeciesB = 1:Nspecies
            speciesFactor = sqrt(THats(iSpeciesA)*mHats(iSpeciesB) / (THats(iSpeciesB) * mHats(iSpeciesA)));
            for ix = ixMin:Nx
                xb = x(ix) * speciesFactor;
                for j = 1:Nx
                    integrationPower = L+2;
                    I_2pL = integral(@integrandWithPower, 0, xb);
                    
                    integrationPower = L+4;
                    I_4pL = integral(@integrandWithPower, 0, xb);
                    
                    partition = max([10,2*xb]);
                    
                    integrationPower = 1-L;
                    I_1mL = integral(@integrandWithPower, xb, partition) + integral(@integrandWithPower, partition, Inf);
                    
                    integrationPower = 3-L;
                    I_3mL = integral(@integrandWithPower, xb, partition) + integral(@integrandWithPower, partition, Inf);
                    
                    % Thes next equations can be found in my notes 20150330-03:
                    
                    tempMatrix_H(ix,j) = 4*pi/(2*L+1) * (I_2pL/(xb^(L+1)) + (xb^L)*I_1mL);
                    
                    tempMatrix_dHdxb(ix,j) = 4*pi/(2*L+1) ...
                        * (-(L+1)*I_2pL/(xb^(L+2)) + L*(xb^(L-1))*I_1mL);
                    
                    tempMatrix_d2Gdxb2(ix,j) = -4*pi/(4*L*L-1) * ( ...
                        L*(L-1)*(xb^(L-2))*I_3mL ...
                        + alpha*(L+1)*(L+2)*(xb^L)*I_1mL ...
                        + alpha*(L+1)*(L+2)/(xb^(L+3))*I_4pL ...
                        + L*(L-1)/(xb^(L+1))*I_2pL);

                end
            end
            Rosenbluth_H(iSpeciesA, iSpeciesB, L+1, :, :)       = tempMatrix_H * collocation2modal;
            Rosenbluth_dHdxb(iSpeciesA, iSpeciesB, L+1, :, :)   = tempMatrix_dHdxb * collocation2modal;
            Rosenbluth_d2Gdxb2(iSpeciesA, iSpeciesB, L+1, :, :) = tempMatrix_d2Gdxb2 * collocation2modal;
            
            %{
            fprintf('For iSpeceisA=%d, iSpeciesB=%d, L=%d\n',iSpeciesA,iSpeciesB,L)
            fprintf('Rosenbluth_H=\n')
            squeeze(Rosenbluth_H(iSpeciesA,iSpeciesB,L+1,:,:))
            fprintf('Rosenbluth_dHdxb=\n')
            squeeze(Rosenbluth_dHdxb(iSpeciesA,iSpeciesB,L+1,:,:))
            fprintf('Rosenbluth_d2Gdxb2=\n')
            squeeze(Rosenbluth_d2Gdxb2(iSpeciesA,iSpeciesB,L+1,:,:))
            %}
        end
    end
end

fprintf('Done. Took %g sec.\n',toc)
%collocation2modal

    function y = evaluatePolynomial(xx)
        % Note: this function uses 1-based indices for j
        if j == 1
            y = ones(size(xx));
        else
            pjMinus1 = zeros(size(xx));
            pj = ones(size(xx));
            for ii = 1:j-1
                y = (xx-a(ii)).*pj - b(ii)*pjMinus1;
                pjMinus1 = pj;
                pj = y;
            end
        end
    end

    function y = integrandWithPower(x)
        % Note that x.^xGrid_k should not be included in the next line!
        y = (x .^ integrationPower) .* evaluatePolynomial(x) .* exp(-x.*x);
    end
end