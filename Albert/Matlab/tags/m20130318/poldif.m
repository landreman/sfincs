function DM = poldif(x, malpha, B)

%  The function DM =  poldif(x, maplha, B) computes the
%  differentiation matrices D1, D2, ..., DM on arbitrary nodes.
%
%  The function is called with either two or three input arguments.
%  If two input arguments are supplied, the weight function is assumed 
%  to be constant.   If three arguments are supplied, the weights should 
%  be defined as the second and third arguments.
%
%  Input (constant weight):
%
%  x:        Vector of N distinct nodes.
%  malpha:   M, the number of derivatives required (integer).
%  B:        Omitted.
%
%  Note:     0 < M < N-1.
%
%  Input (non-constant weight):
%
%  x:        Vector of N distinct nodes.
%  malpha:   Vector of weight values alpha(x), evaluated at x = x(k).
%  B:        Matrix of size M x N,  where M is the highest 
%            derivative required.  It should contain the quantities 
%            B(ell,j) = beta(ell,j) = (ell-th derivative
%            of alpha(x))/alpha(x),   evaluated at x = x(j).
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.

%  J.A.C. Weideman, S.C. Reddy 1998

       N = length(x);                      
       x = x(:);                     % Make sure x is a column vector 

if nargin == 2                       % Check if constant weight function
       M = malpha;                   % is to be assumed.
   alpha = ones(N,1);              
       B = zeros(M,N);
elseif nargin == 3
   alpha = malpha(:);                % Make sure alpha is a column vector
       M = length(B(:,1));           % First dimension of B is the number 
end                                  % of derivative matrices to be computed
 
        I = eye(N);                  % Identity matrix.
        L = logical(I);              % Logical identity matrix.

       XX = x(:,ones(1,N));
       DX = XX-XX';                  % DX contains entries x(k)-x(j).

    DX(L) = ones(N,1);               % Put 1's one the main diagonal.

        c = alpha.*prod(DX,2);       % Quantities c(j).

        C = c(:,ones(1,N)); 
        C = C./C';                   % Matrix with entries c(k)/c(j).
   
        Z = 1./DX;                   % Z contains entries 1/(x(k)-x(j))
     Z(L) = zeros(N,1);              % with zeros on the diagonal.

        X = Z';                      % X is same as Z', but with 
     X(L) = [];                      % diagonal entries removed.
        X = reshape(X,N-1,N);

        Y = ones(N-1,N);             % Initialize Y and D matrices.
        D = eye(N);                  % Y is matrix of cumulative sums,
                                     % D differentiation matrices.
for ell = 1:M
        Y   = cumsum([B(ell,:); ell*Y(1:N-1,:).*X]); % Diagonals
        D   = ell*Z.*(C.*repmat(diag(D),1,N) - D);   % Off-diagonals
     D(L)   = Y(N,:);                                % Correct the diagonal
DM(:,:,ell) = D;                                     % Store the current D
end
