%=== EX_Laplace ===%
% example using the Laplacian in 2D

close all
clear all

N = 99^2;
A = MAT_SymLaplace(N);
f = -ones(N,1);
u_exact = A \ f;

x = linspace(-1,1,101); x=x(2:end-1);
xx= kron(x,ones(1,99));
yy= kron(ones(1,99),x);
ind=1:N;

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
% [u_FOM,r_global_FOM,r_local_FOM] = ALGO_trAOSM(A,f,{ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8,ind9,indtr},u0,101,101);
[u_CG,r_global_CG,r_local_CG] = ALGO_trAOSM_CG(A,f,{ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8,ind9,indtr},u0,30,1e-16);

%%
figure(1)
% subplot(1,2,1)
% surf(reshape(u_exact,101,101))
% subplot(1,2,2)
% surf(reshape(u,101,101))
subplot(1,2,1)
surf(reshape(u_exact - u_FOM, 99,99))
subplot(1,2,2)
surf(reshape(u_exact - u_CG, 99,99))

figure(2)
semilogy(1:length(r_global_FOM),r_global_FOM,'r.--',...
    1:length(r_global_CG),r_global_CG,'bo--')
legend('FOM','CG')

figure(3)
subplot(1,2,1)
r=r_local_FOM{1}; r=r(r>0);
semilogy(1:length(r),r,'--')
hold on
for i=2:length(r_local_FOM)
    r=r_local_FOM{i}; r=r(r>0);
    semilogy(1:length(r),r,'--')
end
hold off
subplot(1,2,2)
r=r_local_CG{1}; r=r(r>0);
semilogy(1:length(r),r,'--')
hold on
for i=2:length(r_local_CG)
    r=r_local_CG{i}; r=r(r>0);
    semilogy(1:length(r),r,'--')
end
hold off