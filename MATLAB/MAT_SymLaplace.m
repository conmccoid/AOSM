function A = MAT_SymLaplace(N)
% SYMLAPLACE generates the symmetric Laplacian matrix of size NxN for 
% domain [-1,1]^2
% A=SymLaplace(N) constructs the Laplacian for the domain [-1,1]^2 with
%   Dirichlet boundary conditions, such that A is of size NxN. N must be a
%   square.
%
%   The matrix A can be used for any domain with h=2/(sqrt(N)-1).
%   
%   Laplace also produces optimized transmission conditions T for 
%   non-overlapping domains.

n=sqrt(N); %NN=(n-1)^2;
h=2/(n+1);
d=ones(n,1)/h^2;
d=[d,-2*d,d];
d=spdiags(d,[-1,0,1],n,n);

I = speye(n);
A=kron(d,I)+kron(I,d);

end