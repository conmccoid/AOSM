%=== EX_Laplace ===%
% example using the Laplacian in 2D

close all
clear all

N = 101^2;
[A,ind] = MAT_Laplace(N);
f = -ones(N,1); f(ind)=0;
u_exact = A \ f;

x = linspace(-1,1,101);
xx= kron(x,ones(1,101));
yy= kron(ones(1,101),x);
ind=1:N;

%% 2 subdomains
% method can't handle a trace with an interior
% ind1 = ind(xx<-0.1);
% ind2 = ind(xx>0.1);
% indtr= ind(xx>=-0.1 & xx<=0.1);
ind1 = ind(xx<0);
ind2 = ind(xx>0);
indtr= ind(xx==0);

% S1 = -A(indtr,[ind2,ind3]) * ( A([ind2,ind3],[ind2,ind3]) \ A([ind2,ind3],indtr) );
% T2 = -A(indtr,ind2) * ( A(ind2,ind2) \ A(ind2,indtr) );
% T3 = -A(indtr,ind3) * ( A(ind3,ind3) \ A(ind3,indtr) );
% surf(log10(abs(S1 - T2 - T3)./abs(S1)))
 
% A11 = A(ind1,ind1);
% A22 = A(ind2,ind2);
% A33 = A(ind3,ind3);
% A1t = A(ind1,indtr);
% A2t = A(ind2,indtr);
% A3t = A(ind3,indtr);
% At1 = A(indtr,ind1);
% At2 = A(indtr,ind2);
% At3 = A(indtr,ind3);

[u,err] = ALGO_trAOSM(A,f,{ind1,ind2,indtr},rand(length(indtr),1));

figure(1)
subplot(1,2,1)
surf(reshape(u_exact,101,101))
subplot(1,2,2)
surf(reshape(u,101,101))

figure(2)
semilogy(err,'r.-')
pause
%% 3 subdomains
close all

ind1 = ind(xx<0);
ind2 = ind(xx>0 & yy>0);
ind3 = ind(xx>0 & yy<0);
indtr= ind(xx==0 | (xx>0 & yy==0));

[u,err] = ALGO_trAOSM(A,f,{ind1,ind2,ind3,indtr},rand(length(indtr),1));

figure(1)
% subplot(1,2,1)
% surf(reshape(u_exact,101,101))
% subplot(1,2,2)
% surf(reshape(u,101,101))
surf(reshape(u_exact - u, 101,101))

figure(2)
semilogy(err,'r.-')
pause
%% 4 subdomains
% this doesn't converge with one crosspoint, only when the subdomains slip
close all

ind1 = ind(xx<0 & yy<-0.1);
ind2 = ind(xx>0 & yy>0.1);
ind3 = ind(xx>0 & yy<0.1);
ind4 = ind(xx<0 & yy>-0.1);
indtr= ind(xx==0 | (xx<0 & yy==-0.1) | (xx>0 & yy==0.1));

[u,err] = ALGO_trAOSM(A,f,{ind1,ind2,ind3,ind4,indtr},rand(length(indtr),1));

figure(1)
% subplot(1,2,1)
% surf(reshape(u_exact,101,101))
% subplot(1,2,2)
% surf(reshape(u,101,101))
surf(reshape(u_exact - u, 101,101))

figure(2)
semilogy(err,'r.-')
pause
%% 4 subdomains, with 5th for crosspoint
% there is some instability here: with bad initial guesses (criterion
% unclear) the convergence rate slows significantly
close all

ind1 = ind(xx<0 & yy<0);
ind2 = ind(xx>0 & yy>0);
ind3 = ind(xx>0 & yy<0);
ind4 = ind(xx<0 & yy>0);
ind5 = ind(xx==0 & yy==0);
indtr= ind(xor(xx==0,yy==0));

u0 = rand(length(indtr),1);
% u0 = f(indtr); % needs to be modified
[u,err] = ALGO_trAOSM(A,f,{ind1,ind2,ind3,ind4,ind5,indtr},u0);

figure(1)
% subplot(1,2,1)
% surf(reshape(u_exact,101,101))
% subplot(1,2,2)
% surf(reshape(u,101,101))
surf(reshape(u_exact - u, 101,101))

figure(2)
semilogy(err,'r.-')
pause
%% 9 subdomains, with 10th for crosspoints
% unstable, apparently heavily dependent on initial guess
% re-stabilized by putting the crosspoints in larger subdomains
close all
ind1 = ind((xx<-0.3 & yy<-0.3) | (xx==-0.3 & yy==-0.3));
ind2 = ind(xx>-0.3 & xx<0.3 & yy<-0.3);
ind3 = ind((xx>0.3 & yy<-0.3) | (xx==0.3 & yy==-0.3));
ind4 = ind(xx<-0.3 & yy>-0.3 & yy<0.3);
ind5 = ind(xx>-0.3 & xx<0.3 & yy>-0.3 & yy<0.3);
ind6 = ind(xx>0.3 & yy>-0.3 & yy<0.3);
ind7 = ind((xx<-0.3 & yy>0.3) | (xx==-0.3 & yy==0.3));
ind8 = ind(xx>-0.3 & xx<0.3 & yy>0.3);
ind9 = ind((xx>0.3 & yy>0.3) | (xx==0.3 & yy==0.3));
% ind10= ind((xx==0.3 & yy==0.3) | (xx==0.3 & yy==-0.3) | (xx==-0.3 & yy==0.3) | (xx==-0.3 & yy==-0.3));
indtr= ind(xor( xx==-0.3 | xx==0.3,yy==-0.3 | yy==0.3 ));

u0 = rand(length(indtr),1);
[u,err] = ALGO_trAOSM(A,f,{ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8,ind9,indtr},u0);

figure(1)
% subplot(1,2,1)
% surf(reshape(u_exact,101,101))
% subplot(1,2,2)
% surf(reshape(u,101,101))
surf(reshape(u_exact - u, 101,101))

figure(2)
semilogy(err,'r.-')
pause
%% 3 strips
close all
ind1 = ind(xx<-0.3);
ind2 = ind(xx>-0.3 & xx<0.3);
ind3 = ind(xx>0.3);
indtr= ind(xx==-0.3 | xx==0.3);

[u,err] = ALGO_trAOSM(A,f,{ind1,ind2,ind3,indtr},rand(length(indtr),1));

figure(1)
% subplot(1,2,1)
% surf(reshape(u_exact,101,101))
% subplot(1,2,2)
% surf(reshape(u,101,101))
surf(reshape(u_exact - u, 101,101))

figure(2)
semilogy(err,'r.-')
pause
%% 4 strips
close all
ind1 = ind(xx<-0.5);
ind2 = ind(xx>-0.5 & xx<0);
ind3 = ind(xx>0 & xx<0.5);
ind4 = ind(xx>0.5);
indtr= ind(xx==-0.5 | xx==0.5 | xx==0);

[u,err] = ALGO_trAOSM(A,f,{ind1,ind2,ind3,ind4,indtr},rand(length(indtr),1));

figure(1)
% subplot(1,2,1)
% surf(reshape(u_exact,101,101))
% subplot(1,2,2)
% surf(reshape(u,101,101))
surf(reshape(u_exact - u, 101,101))

figure(2)
semilogy(err,'r.-')
pause