function [A,ind,T1,T2] = MAT_Laplace(N)
% LAPLACE generates the Laplacian matrix of size NxN for domain [-1,1]^2
% A=Laplace(N) constructs the Laplacian for the domain [-1,1]^2 with
%   Dirichlet boundary conditions, such that A is of size NxN. N must be a
%   square.
%
%   The matrix A can be used for any domain with h=2/(sqrt(N)-1).
%   
%   Laplace also produces optimized transmission conditions T for 
%   non-overlapping domains.

n=sqrt(N); %NN=(n-1)^2;
h=2/(n-1);
d=ones(n,1)/h^2;
d=[d,-2*d,d];
d=spdiags(d,[-1,0,1],n,n);
d(1,:)=0; d(end,:)=0;

I=speye(n); Iint=I; Iint(1)=0; Iint(end)=0; Iext=I-Iint;
Boundaries=kron(Iint,Iext)+kron(Iext,Iint)+kron(Iext,Iext);
A=kron(d,Iint)+kron(Iint,d)+Boundaries;
ind=1:N; ind=ind(diag(Boundaries)==1);

if nargout>2
    p=pi/sqrt(h)^3;
    T1=ones(n,1)/h^2;
    T1=[0.5*T1,-2*T1,0.5*T1];
    T1=spdiags(T1,[-1,0,1],n,n);
    T2=-T1 - p*speye(n);
    T1=-T1 - p*speye(n);
    T1(1,:)=0; T1(end,:)=0;
    T2(1,:)=0; T2(end,:)=0;
end