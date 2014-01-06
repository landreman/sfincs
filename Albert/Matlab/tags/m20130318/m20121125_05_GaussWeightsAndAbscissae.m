function [x,w]=m20121125_05_GaussWeightsAndAbscissae(N, weight, xMin, xMax)
% This script computes the nodes and weights for Gaussian quadrature for any
% weight function on any interval.  The algorithm used is the Stieltjes
% procedure described in section 4.6.3 of Numerical Recipes, 3rd edition,
% together with the Jacobi matrix procedure of section 4.6.2.  See also
% www.nr.com/webnotes?3
%
% Inputs:
% N = number of nodes
% weight = function handle to the weight function
% [xmin, xmax] = interval for the integration. -Inf and Inf are allowed.
% 
% Outputs:
% x = nodes (i.e. abscissae) in a column vector
% w = weights in a column vector
%
% To compute the integral of any function f(), just form w' * f(x)
%
% Written by Matt Landreman - Dec 17, 2011
%

a = zeros(1,N);
b = zeros(1,N);
c = zeros(1,N);
d = zeros(1,N);

intervalDivider=10;
oldc=1;
for i=1:N
    j=i-1;
    if isinf(xMax)
        c(i) = quadgk(@integrandWithoutX, xMin, 5)+quadgk(@integrandWithoutX, 5, intervalDivider)+quadgk(@integrandWithoutX, intervalDivider, Inf);
        d(i) = quadgk(@integrandWithX, xMin, 5)+quadgk(@integrandWithX, 5, intervalDivider)+quadgk(@integrandWithX, intervalDivider, Inf); 
    else
        c(i) = quadgk(@integrandWithoutX, xMin, xMax);
        d(i) = quadgk(@integrandWithX, xMin, xMax);
    end
    b(i) = c(i)/oldc;
    a(i) = d(i)/c(i);
    oldc=c(i);
    if isnan(c(i))
        fprintf('Nan in c for i=%d.\n',i)
    end
    if isnan(d(i))
        fprintf('Nan in d for i=%d.\n',i)
        figure(33)
        clf
        xs=linspace(0,14,400);
        plot(xs, integrandWithX(xs));
        hold on
        plot(xs, integrandWithoutX(xs),'r');
        
    end
end



    function y=evaluatePolynomial(x)
        if j==0
            y=ones(size(x));
        else
            pjMinus1 = zeros(size(x));
            pj = ones(size(x));
            for ii=1:j
                y = (x-a(ii)).*pj - b(ii)*pjMinus1;
                pjMinus1 = pj;
                pj = y;
            end
        end
    end

    function z=integrandWithoutX(x)
        p = evaluatePolynomial(x);
        z=p.*weight(x).*p;
        ii=isnan(z);
        z(ii)=0;
    end

    function z=integrandWithX(x)
        p = evaluatePolynomial(x);
        z=x.*p.*weight(x).*p;
        ii=isnan(z);
        z(ii)=0;
    end

sqrtb=sqrt(b);
% Discard the first element:
sqrtb(1)=[];
% Form the tri-diagonal matrix:
JacobiMatrix = diag(a) + diag(sqrtb,1) + diag(sqrtb,-1);

[eigenvectors, eigenvalues] = eig(JacobiMatrix);
x=diag(eigenvalues);
w=eigenvectors(1,:)';
w=c(1)*w.*w;
end