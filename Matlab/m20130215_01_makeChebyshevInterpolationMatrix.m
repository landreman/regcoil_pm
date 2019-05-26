function matrix = m20130215_01_makeChebyshevInterpolationMatrix(N, xMin, xMax, x)
%function p = chebint2(fk, xMin, xMax, x)

%  The function p = chebint(fk, , xMin, xMax, x) computes the polynomial interpolant
%  of the data (xk, fk), where xk are the Chebyshev nodes on [xMin, xMax].
%  Two or more data points are assumed.
%
%  Input:
%  fk:  Vector of y-coordinates of data, at Chebyshev points 
%       x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
%  x:   Vector of x-values where polynomial interpolant is to be evaluated.
%
%  Output:
%  p:    Vector of interpolated values.
%
%  The code implements the barycentric formula; see page 252 in
%  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
%  (Note that if some fk > 1/eps, with eps the machine epsilon,
%  the value of eps in the code may have to be reduced.)

%  J.A.C. Weideman, S.C. Reddy 1998

% Beginning of changes by MJL:
if any(x>xMax)
    error('All entries in x must be <= xMax.')
end
if any(x<xMin)
    error('All entries in x must be >= xMin.')
end

if N==1
    matrix = ones(numel(x),1);
    return;
end

xMid = 0.5*(xMin+xMax);
x = -2*(x-xMid)/(xMax-xMin);
% End of changes by MJL.

  %fk = fk(:); 
  x = x(:);                    % Make sure data are column vectors.

   %N = length(fk); 
   M = length(x);
     
  xk = sin(pi*[N-1:-2:1-N]'/(2*(N-1)));    % Compute Chebyshev points.

   w = ones(N,1).*(-1).^[0:N-1]';          % w = weights for Chebyshev formula
w(1) = w(1)/2; w(N) = w(N)/2;
 
x(:,ones(1,N))
xk(:,ones(1,M))'
   D = x(:,ones(1,N)) - xk(:,ones(1,M))';  % Compute quantities x-x(k)
   D = 1./(D+eps*(D==0));                  % and their reciprocals.
  
   %p = D*(w.*fk)./(D*w);                   % Evaluate interpolant as
                                           % matrix-vector products.
                                           fprintf('D:\n')
                                           D
                                           fprintf('w:\n')
                                           w
                                           fprintf('Denominator:\n')
                                           denominator=D*w
matrix = diag(1./(D*w)) * D * diag(w);