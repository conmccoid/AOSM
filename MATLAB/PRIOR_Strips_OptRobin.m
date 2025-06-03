%===Strips: optimizing Robin parameter===%
% needs to be re-done with restricted trace procedures

ny=11;
h=0.1;
n_strips=3;
nx=n_strips*(ny-1) + 1;
[A,ind] = MAT_Laplace_rectangle(nx,ny,h);
f = -ones(nx*ny,1); f(ind)=0;
u_exact = A \ f;

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

list_p=-pi/sqrt(h)^3 + (-10:2:10);
% list_p=(n_strips-1)*list_p;
% list_p=-[200,100,50,25,5,1,0];
legend_list=cell(length(list_p),1);
conv_rate=zeros(length(list_p),1);
k=1;
for p=list_p
    u0 = rand(length(indtr),1);
    u0(ismember(indtr,ind))=0;
    [u,err] = ALGO_trOSM(A,f,ind_strips,u0,p,nx,ny);
    P=polyfit(1:length(err),log(err),1);
    conv_rate(k)=exp(P(1));

    figure(2)
    semilogy(err,'.--')
    xlabel('Iteration')
    ylabel('Error')
    legend_list{k}=num2str(p);
    k=k+1;
    hold on
end
legend(legend_list)
hold off

%%
figure(3)
semilogy(list_p,conv_rate,'r.--')
xlabel('p')
ylabel('Convergence rate')