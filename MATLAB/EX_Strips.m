%=== EX_Strips ===%
% example using the Laplacian in 2D

close all
clear all

n=41;
N = n^2;
[A,ind] = MAT_Laplace(N);
f = -ones(N,1); f(ind)=0;
u_exact = A \ f;

x = linspace(-1,1,n);
xx= kron(x,ones(1,n));
ind=1:N;

n_strips=2;
ind_strips=cell(n_strips+1,1);
x_step = floor(n/n_strips)+1;
x_strips=x(x_step:x_step:end); x_strips=[x_strips,1.01];
ind_strips{1} = ind(xx<x_strips(1));
indtr=ind(xx==x_strips(1));
for strip=2:n_strips
    ind_strips{strip} = ind(xx>x_strips(strip-1) & xx<x_strips(strip));
    indtr=[indtr,ind(xx==x_strips(strip))];
end
ind_strips{n_strips+1} = indtr;

[u,err] = ALGO_trAOSM(A,f,ind_strips,rand(n*(n_strips-1),1));

figure(1)
subplot(1,2,1)
surf(reshape(u_exact,n,n))
subplot(1,2,2)
surf(reshape(u,n,n))

figure(2)
semilogy(err,'r.-')
pause

for n=[81,161,321,641]
    N=n^2;
    [A,ind] = MAT_Laplace(N);
    f = -ones(N,1); f(ind)=0;
    u_exact = A \ f;

    x = linspace(-1,1,n);
    xx= kron(x,ones(1,n));
    ind=1:N;

    n_strips=8*n_strips;
    ind_strips=cell(n_strips+1,1);
    x_step = floor(n/n_strips)+1;
    x_strips=x(x_step:x_step:end); x_strips=[x_strips,1.01];
    ind_strips{1} = ind(xx<x_strips(1));
    indtr=ind(xx==x_strips(1));
    for strip=2:n_strips
        ind_strips{strip} = ind(xx>x_strips(strip-1) & xx<x_strips(strip));
        indtr=[indtr,ind(xx==x_strips(strip))];
    end
    ind_strips{n_strips+1} = indtr;
    
    [u,err] = ALGO_trAOSM(A,f,ind_strips,rand(n*(n_strips-1),1));
    
    figure(1)
    subplot(1,2,1)
    surf(reshape(u_exact,n,n))
    subplot(1,2,2)
    surf(reshape(u,n,n))
    
    figure(2)
    semilogy(err,'r.-')
    pause
end
