%=== EX_Blocks ===%
% example using the Laplacian in 2D

%% Weak scaling

close all
clear all

list_blocks=2:5; block_size=10;
conv_rate = zeros(length(list_blocks),1); k=1;
err_master = cell(length(list_blocks),1);
for n_blocks=list_blocks
    ny=n_blocks*block_size + 1;
    nx=ny;
    h=1/(ny-1);

    [A,ind] = MAT_Laplace_rectangle(nx,ny,h);
    f = -ones(nx*ny,1); f(ind)=0;
    u_exact = A \ f;

%     figure(1)
%     surf(reshape(u_exact,ny,nx))
%     pause

    % ind_blocks subroutine
    x = linspace(0,1,nx);
    xx= kron(x,ones(1,nx));
    yy= kron(ones(1,nx),x);
    ind_master=1:nx*ny;

    ind_blocks=cell(n_blocks+1,1);
    ind_tr=[];
    count=1;
    x(1)=-1; x(end)=2;
    for i=1:n_blocks
        ind_x = xx < x(i*block_size+1) & xx > x((i-1)*block_size+1);
        for j=1:n_blocks
            ind_y = yy < x(j*block_size+1) & yy > x((j-1)*block_size+1);
            block_int=ind_master(ind_x & ind_y);
            crosspoint=ind_master(xx==x(i*block_size+1) & yy==x(j*block_size+1));
            ind_blocks{count}=[block_int, crosspoint]; count=count+1;
            ind_tr=[ind_tr,ind_master(xx==x(i*block_size+1) & ind_y)];
            ind_tr=[ind_tr,ind_master(ind_x & yy==x(j*block_size+1))];
        end
    end
    ind_blocks{count}=ind_tr;

    [u,err] = ALGO_trAOSM(A,f,ind_blocks,rand(length(ind_tr),1),nx,ny);
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
plot(list_blocks.^2,abs(conv_rate),'bo--')
xlabel('Number of blocks')
ylabel('Convergence rate')

figure
legend_list=cell(length(list_blocks),1);
for i=1:length(list_blocks)
    semilogy(err_master{i},'.--')
    xlabel('Iteration')
    ylabel('Error')
    legend_list{i}=num2str(list_blocks(i));
    hold on
end
legend(legend_list)
hold off

%% Strong scaling

close all
clear all

list_blocks=2:5; block_size=3*4;
conv_rate = zeros(length(list_blocks),1); k=1;
err_master = cell(length(list_blocks),1);

ny=max(list_blocks)*block_size + 1;
nx=ny;
h=1/(ny-1);

[A,ind] = MAT_Laplace_rectangle(nx,ny,h);
f = -ones(nx*ny,1); f(ind)=0;
u_exact = A \ f;

for n_blocks=list_blocks

    block_size=(ny-1)/n_blocks;

%     figure(1)
%     surf(reshape(u_exact,ny,nx))

    % ind_blocks subroutine
    x = linspace(0,1,nx);
    xx= kron(x,ones(1,nx));
    yy= kron(ones(1,nx),x);
    ind_master=1:nx*ny;

    ind_blocks=cell(n_blocks+1,1);
    ind_tr=[];
    count=1;
    x(1)=-1; x(end)=2;
    for i=1:n_blocks
        ind_x = xx < x(i*block_size+1) & xx > x((i-1)*block_size+1);
        for j=1:n_blocks
            ind_y = yy < x(j*block_size+1) & yy > x((j-1)*block_size+1);
            block_int=ind_master(ind_x & ind_y);
            crosspoint=ind_master(xx==x(i*block_size+1) & yy==x(j*block_size+1));
            ind_blocks{count}=[block_int, crosspoint]; count=count+1;
            ind_tr=[ind_tr,ind_master(xx==x(i*block_size+1) & ind_y)];
            ind_tr=[ind_tr,ind_master(ind_x & yy==x(j*block_size+1))];
        end
    end
    ind_blocks{count}=ind_tr;

    [u,err] = ALGO_trAOSM(A,f,ind_blocks,rand(length(ind_tr),1),nx,ny);
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
plot(list_blocks.^2,abs(conv_rate),'bo--')
xlabel('Number of blocks')
ylabel('Convergence rate')

figure
legend_list=cell(length(list_blocks),1);
for i=1:length(list_blocks)
    semilogy(err_master{i},'.--')
    xlabel('Iteration')
    ylabel('Error')
    legend_list{i}=num2str(list_blocks(i));
    hold on
end
legend(legend_list)
hold off