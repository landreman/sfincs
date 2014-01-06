function interpMatrix = m20120703_03_polynomialInterpolationMatrix(xk, x, alpxk, alpx)
% Returns the matrix for polynomial interpolation from the xk grid to the x
% grid.
%
% This script is a modification of the polint.m routine from DMSuite,
% by J.A.C. Weideman and S.C. Reddy, available here:
% http://www.mathworks.com/matlabcentral/fileexchange/29
% or here:
% http://dip.sun.ac.za/~weideman/research/differ.html
%
% The original polint.m routine performed the interpolation explicitly.
% This modified function returns the interpolation matrix, for use in
% implicit codes.
% 
% Modification by Matt Landreman, MIT Plasma Science and Fusion Center,
% July 3, 2012.
%
%  Inputs
%  xk:    Vector of x-coordinates of data (assumed distinct).
%  x:     Vector of x-values where polynomial interpolant is to be evaluated.
%  alpxk: Vector of weight values sampled at the points xk.
%  alpx:  Vector of weight values sampled at the points x.
%
%  Output:
%  interpolation matrix
%
%  The code implements the barycentric formula; see page 252 in
%  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
%  (Note that if some fk > 1/eps, with eps the machine epsilon,
%  the value of eps in the code may have to be reduced.)
%  Note added May 2003:  Except for certain nice node distributions
%  polynomial interpolation of high-degree is an ill-conditioned
%  problem.  This code does not test for conditioning so use with
%  care. 

%  J.A.C. Weideman, S.C. Reddy 1998.

N = length(xk);
     M = length(x);

if nargin == 2
    alpx = 1; alpxk = 1;
    deWeightify = eye(N);
elseif nargin > 2
    deWeightify = diag(1./alpxk);
else
    error('Unexpected number of inputs.')
end

     x = x(:);                          % Make sure the data are column vectors
    xk = xk(:);  
 alpxk = alpxk(:); alpx = alpx(:);

     L = logical(eye(N));

     D = xk(:,ones(1,N))-xk(:,ones(1,N))';  % Compute the weights w(k)
  D(L) = ones(N,1);
     w = 1./prod(D)';                       
 
     
     D = x(:,ones(1,N)) - xk(:,ones(1,M))'; % Compute quantities x-x(k) 
     D = 1./(D+eps*(D==0));                % and their reciprocals. 
  
     % Evaluate interpolant as matrix-vector products:
     interpMatrix = diag(alpx./(D*w))*D*diag(w)*deWeightify;
