%=== EX_Strips ===%
% example using the Laplacian in 2D

%% Weak scaling

close all
clear all

ny=11;
h=0.1;

list_strips=2:5; list_strips=2.^list_strips;
conv_rate = zeros(length(list_strips),1); k=1;
err_master = cell(length(list_strips),1);
for n_strips=list_strips
    nx=n_strips*(ny-1) + 1;
    [A,ind] = MAT_Laplace_rectangle(nx,ny,h);
    f = -ones(nx*ny,1); f(ind)=0;
    u_exact = A \ f;

%     figure(1)
%     surf(reshape(u_exact,ny,nx))

    % ind_strips subroutine
    ind_strips=cell(n_strips+1,1);
    ind_strip=1:ny*(ny-2);
    ind_strips{1}=1:ny*(ny-1);
    ind_i=ny*(ny-1);
    indtr=ind_i+(1:ny);
    ind_i=ind_i+ny;
    for i=2:n_strips-1
        ind_strips{i}=ind_i+ind_strip;
        ind_i=ind_i+ny*(ny-2);
        indtr=[indtr,ind_i+(1:ny)];
        ind_i=ind_i+ny;
    end
    ind_strips{n_strips}=ind_i+(1:ny*(ny-1));
    ind_strips{end}=indtr;

    [u,err] = ALGO_trAOSM(A,f,ind_strips,rand(length(indtr),1),nx,ny);
    P = polyfit(1:length(err),log(err),1);
    conv_rate(k)=exp(P(1));
    err_master{k}=err;
    k=k+1;

%     figure(2)
%     surf(reshape(u,ny,nx))
%     figure(3)
%     semilogy(err,'r.--')
% 
%     pause
end

figure
semilogx(list_strips,abs(conv_rate),'bo--')
xlabel('Number of strips')
ylabel('Convergence rate')

figure
legend_list=cell(length(list_strips),1);
for i=1:length(list_strips)
    semilogy(err_master{i},'.--')
    xlabel('Iteration')
    ylabel('Error')
    legend_list{i}=num2str(list_strips(i));
    hold on
end
legend(legend_list)
hold off

%% Strong scaling

close all
clear all

ny=11;
h=0.1;

list_strips=2:5; list_strips=2.^list_strips;
conv_rate = zeros(length(list_strips),1); k=1;
err_master = cell(length(list_strips),1);

nx=max(list_strips)*(ny-1) + 1;
[A,ind] = MAT_Laplace_rectangle(nx,ny,h);
f = -ones(nx*ny,1); f(ind)=0;
u_exact = A \ f;

for n_strips=list_strips

%     figure(1)
%     surf(reshape(u_exact,ny,nx))

    % ind_strips subroutine - nb: redo for strong scaling
    ind_strips=cell(n_strips+1,1);
    ind_strip=1:ny*(ny-2);
    ind_strips{1}=1:ny*(ny-1);
    ind_i=ny*(ny-1);
    indtr=ind_i+(1:ny);
    ind_i=ind_i+ny;
    for i=2:n_strips-1
        ind_strips{i}=ind_i+ind_strip;
        ind_i=ind_i+ny*(ny-2);
        indtr=[indtr,ind_i+(1:ny)];
        ind_i=ind_i+ny;
    end
    ind_strips{n_strips}=ind_i+(1:ny*(ny-1));
    ind_strips{end}=indtr;

    [u,err] = ALGO_trAOSM(A,f,ind_strips,rand(length(indtr),1),nx,ny);
    P = polyfit(1:length(err),log(err),1);
    conv_rate(k)=exp(P(1));
    err_master{k}=err;
    k=k+1;

%     figure(2)
%     surf(reshape(u,ny,nx))
%     figure(3)
%     semilogy(err,'r.--')
% 
%     pause
end

figure
semilogx(list_strips,abs(conv_rate),'bo--')
xlabel('Number of strips')
ylabel('Convergence rate')

figure
legend_list=cell(length(list_strips),1);
for i=1:length(list_strips)
    semilogy(err_master{i},'.--')
    xlabel('Iteration')
    ylabel('Error')
    legend_list{i}=num2str(list_strips(i));
    hold on
end
legend(legend_list)
hold off