% Vincente Pericoli
% UC Davis
% 25 Jan 2015
%
% Finite difference approximant for non-uniformly spaced 1D grid.
% Uses the Fornberg method for all points to calculate almost-arbitrarily 
% high-order derivatives. It uses a stencil size of 5 points in the 
% interior, and 6 for each of the two boundary points.
%
%
% Fornberg, Bengt. "Calculation of Weights in Finite Difference Formulas." 
%       SIAM Review 40, no. 3 (January 1998): 685-91.
%       http://dx.doi.org/10.1137/S0036144596322507.
%
% Fornberg, Bengt. "Generation of finite difference formulas on arbitrarily
%       spaced grids." Mathematics of Computation 51, no. 184 (1988): 
%       699-99. http://dx.doi.org/10.1090/S0025-5718-1988-0935077-0
%
%
% This code is a  modified version of:
%       http://www.mathworks.com/matlabcentral/fileexchange/50415
%       (c) 2015 Jess (BSD License)
%

function [ du ] = mderiv_fornberg( m, X, u )
%NDERIV_FORNBERG
% This routine computes the m^th derivative of a vector (u) defined on a 
% set of points (X) using Fornberg's algorithm.
%
% Inputs:
%   m = the order of the derivative you wish to compute
%   X = vector set of points on which u is evaluated
%   u = vector of values u(x) which are evaluated on X (u(i) := u(X(i)))
% 
% Requires/Assumes:
%   * X should be sorted in ascending order
%   * m < length(X)
% RESTRICTIONS: 
%   * k must be less than Nx, but greater than 1
%   * length(X) must be equal to or greater than 6

% check inputs
Nx = length(X);
if m >= Nx, error('*** length(x) must be larger than k'); end
if Nx ~= length(u), error('*** length(x) must be equal length(u)'); end

% preallocate
du = zeros(size(u));

% Low boundary 2 points
du(1) = dot( fornweights(m,X(1),X(1:6)), u(1:6) );
du(2) = dot( fornweights(m,X(2),X(1:6)), u(1:6) );

% interior points
for ix = 3:Nx-2
   du(ix) = dot( fornweights(m,X(ix),X(ix-2:ix+2)), u(ix-2:ix+2) );
end

% Upper boundary 2 points
du(Nx-1) = dot( fornweights(m,X(Nx-1),X(Nx-5:Nx)), u(Nx-5:Nx) );
du(Nx  ) = dot( fornweights(m,X(Nx  ),X(Nx-5:Nx)), u(Nx-5:Nx) );

end


function c = fornweights(m,xbar,x)
% Compute coefficients for finite difference approximation for the
% derivative of order k at xbar based on grid values at points in x.
%
% This function returns a row vector c of dimension 1 by n, where
% n = length(x), containing coefficients to approximate u^{(k)}(xbar), 
% the k'th derivative of u evaluated at xbar, based on n values
% of u at x(1), x(2), ... x(n).  
%
% If U is a column vector containing u(x) at these n points, then 
% c*U will give the approximation to u^{(k)}(xbar).
%
% Note for k=0 this can be used to evaluate the interpolating polynomial 
% itself.
%
% Requires length(x) > k.  
% Usually the elements x(i) are monotonically increasing
% and x(1) <= xbar <= x(n), but neither condition is required.
% The x values need not be equally spaced but must be distinct.  
%
% This program should give the same results as fdcoeffV.m, but for large
% values of n is much more stable numerically.
%
% Based on the program "weights" in 
%   B. Fornberg, "Calculation of weights in finite difference formulas",
%   SIAM Review 40 (1998), pp. 685-691.
%
% Note: Forberg's algorithm can be used to simultaneously compute the
% coefficients for derivatives of order 0, 1, ..., m where m <= n-1.
% This gives a coefficient matrix C(1:n,1:m) whose k'th column gives
% the coefficients for the k'th derivative.
%
% In this version we set m=k and only compute the coefficients for
% derivatives of order up to order k, and then return only the k'th column
% of the resulting C matrix (converted to a row vector).  
% This routine is then compatible with fdcoeffV.   
% It can be easily modified to return the whole array if desired.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)


n = length(x);
if m >= n, error('*** length(x) must be larger than m'); end

c1 = 1;
c4 = x(1) - xbar;
C = zeros(n-1,m+1);
C(1,1) = 1;
for i=1:n-1
    i1 = i+1;
    mn = min(i,m);
    c2 = 1;
    c5 = c4;
    c4 = x(i1) - xbar;
    for j=0:i-1
        j1 = j+1;
        c3 = x(i1) - x(j1);
        c2 = c2*c3;
        if j==i-1
            for s=mn:-1:1
                s1 = s+1;
                C(i1,s1) = c1*(s*C(i1-1,s1-1) - c5*C(i1-1,s1))/c2;
            end
            C(i1,1) = -c1*c5*C(i1-1,1)/c2;
        end
        for s=mn:-1:1
            s1 = s+1;
            C(j1,s1) = (c4*C(j1,s1) - s*C(j1,s1-1))/c3;
        end
        C(j1,1) = c4*C(j1,1)/c3;
    end
    c1 = c2;
end

c = C(:,end)'; % last column of c gives desired row vector

return;
end