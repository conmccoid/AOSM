%=== EX_Laplace ===%
% example using the Laplacian in 2D

clear

N = 101^2;
[A,ind] = MAT_Laplace(N);
f = -ones(N,1); f(ind)=0;
u_exact = A \ f;

x = linspace(-1,1,101);
xx= kron(x,ones(1,101));
yy= kron(ones(1,101),x);
ind=1:N;
ind1 = ind(xx<0);
ind2 = ind(xx>0 & yy>0);
ind3 = ind(xx>0 & yy<0);
indtr= ind(xx==0 | (xx>0 & yy==0));
spy(A([ind1,ind2,ind3,indtr],[ind1,ind2,ind3,indtr]))
 
% A11 = A(ind1,ind1);
% A22 = A(ind2,ind2);
% A33 = A(ind3,ind3);
% A1t = A(ind1,indtr);
% A2t = A(ind2,indtr);
% A3t = A(ind3,indtr);
% At1 = A(indtr,ind1);
% At2 = A(indtr,ind2);
% At3 = A(indtr,ind3);

u = ALGO_trAOSM(A,f,{ind1,ind2,ind3,indtr},rand(length(indtr),1));

figure(1)
subplot(1,2,1)
surf(reshape(u_exact,101,101))
subplot(1,2,2)
surf(reshape(u,101,101))