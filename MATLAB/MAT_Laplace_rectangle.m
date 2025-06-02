function [A,ind] = MAT_Laplace_rectangle(nx,ny,h)
% LAPLACE_rectangle generates the Laplacian matrix for a rectangular
% domain, nx points long by ny points tall, with mesh parameter h

dx=ones(nx,1)/h^2;
dy=ones(ny,1)/h^2;
dx=[dx,-2*dx,dx];
dy=[dy,-2*dy,dy];
dx=spdiags(dx,[-1,0,1],nx,nx);
dy=spdiags(dy,[-1,0,1],ny,ny);
dx(1,:)=0; dx(end,:)=0;
dy(1,:)=0; dy(end,:)=0;

Iintx=speye(nx); Iintx(1)=0; Iintx(end)=0;
Iinty=speye(ny); Iinty(1)=0; Iinty(end)=0;
Iextx=speye(nx) - Iintx;
Iexty=speye(ny) - Iinty;

BCs=kron(Iintx,Iexty) + kron(Iextx,Iinty) + kron(Iextx,Iexty);
A = kron(dx,Iinty) + kron(Iintx,dy) + BCs;
ind=1:(nx*ny); ind=ind(diag(BCs)==1);