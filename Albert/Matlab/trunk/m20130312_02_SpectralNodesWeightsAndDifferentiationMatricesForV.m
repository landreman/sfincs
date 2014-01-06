function [x, ddx, d2dx2, weights] = m20130312_02_SpectralNodesWeightsAndDifferentiationMatricesForV(N, k, scale, pointAtZero)
% Computes the abscissae (i.e. grid points), Gaussian integration weights,
% and spectral colocation differentiation matrices for functions on the
% domain [0, Inf] which behave roughly as exp(-x^2) for large x.
%
% Written by Matt Landreman - July 6, 2012
%
% Inputs:
% N = desired number of abscissae
% k = extra power of x in the orthogonality relation and weight. 
%        (Usually, k=0 works well.)
% scale = factor by which to scale all the abscissae, (typically in the range 0.7 to 1)
%               so as to not waste points at very large x.
% pointAtZero = true or false, depending on whether or not you require a
%               grid point at x=0.
%
% Outputs:
% x = column vector of abscissae
% ddx = first derivative matrix
% d2dx2 = second derivative matrix
% weights = Gaussian integration weights
%
% This function depends on
% m20121125_05_GaussWeightsAndAbscissae.m
% as well as
% poldif.m
% The latter is part of the DMSuite package by S.C. Reddy and J.A.C. Weideman, available at
% http://www.mathworks.com/matlabcentral/fileexchange/29
% or here:
% http://dip.sun.ac.za/~weideman/research/differ.html

if pointAtZero
    [x, weights]=m20121007_01_generateGaussRadauPointsAndWeights(N, @(x) (x.^k) .* exp(-x.*x), 0, Inf);
else
    [x, weights]=m20121125_05_GaussWeightsAndAbscissae(N, @(x) (x.^k) .* exp(-x.*x), 0, Inf);
end
weights=weights./((x.^k) .* exp(-x.*x));
x = x*scale;
weights = weights*scale;
weightFactors=[k./x-2*x, k*(k-1)./(x.*x)-2*(2*k+1)+4*x.*x]';
differentiationMatrices=poldif(x,(x.^k) .* exp(-x.*x),weightFactors);
ddx=differentiationMatrices(:,:,1);
d2dx2=differentiationMatrices(:,:,2);

end
