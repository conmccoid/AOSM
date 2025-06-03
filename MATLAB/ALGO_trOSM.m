function [u_master, err_out] = ALGO_trOSM(A,f,sub_indices, u0, p, nx, ny)
% trOSM runs a multi-subdomain optimized Schwarz method that shares the trace with all
% subdomains
%   - A: system matrix
%   - f: right hand side
%   - sub_indices: array (?) of subdomain indices, the last of which is the
%   trace
%   - u0: initial guess on the trace
%   - p: Robin parameter

n = length(sub_indices);        % number of subdomains (+1 for the trace)
N = length(A);                  % size of global problem
ind_tr = sub_indices{n};        % indices for the trace
M = length(ind_tr);             % size of trace
A_tr = A(ind_tr,ind_tr);        % matrix block for trace
f_tr = f(ind_tr);               % right hand side for trace
F_tr = zeros(M,n-1);            % mods to rhs for trace
for i=1:n-1
    ind_i = sub_indices{i};
    F_tr(:,i) = A(ind_tr,ind_i) * (A(ind_i,ind_i) \ f(ind_i));
end
T_master = cell(n,1);           % master list of adapted transmission conditions
t_master = cell(n,1);           % master list for soln on trace
S_master = zeros(n*M,M);            % master list of sums of ATCs
u_master = zeros(N,1);          % master soln of global problem
itermax = M;                    % maximum number of iterations (nb: make input)

for i=1:n-1                     % for each subdomain
    T_master{i} = -0.5*A_tr + p*speye(M); % store adapted transmission conditions
    t_master{i} = u0;
    for j=setdiff(1:n-1,i)
        S_master((1:M)+M*(j-1),1:M) = S_master((1:M)+M*(j-1),1:M) + T_master{i};
    end
end

%===Global iterations===%
r = f_tr - A(ind_tr,:)*u_master;% residual on the trace
u_new = zeros(N,1);               % storage of new soln
t_new = cell(n,1);
iter = 1;                       % init interation count
err_out(iter)=norm(r);
while norm(r)>1e-10 && iter<itermax
    for i=1:n-1
        ind_i = sub_indices{i};     % indices for this subdomain
        N_i = length(ind_i);        % size of subdomain
        A_ii = A(ind_i,ind_i);      % main block for this subdomain
        A_it = A(ind_i,ind_tr);     % top right block
        A_ti = A(ind_tr,ind_i);     % bottom left block
        f_i = f(ind_i);             % right hand side for this subdomain
%         S_i = S_master((1:M)+M*(i-1),1:M);
        S_i = T_master{i};
        y_i = zeros(M,1);
        for j=setdiff(1:n-1,i)
            ind_j = sub_indices{j};
            T_i = T_master{j};
            t_i = t_master{j};
            y_i = y_i - A(ind_tr,ind_j)*u_master(ind_j) + T_i*t_i;
        end
        u_i = [A_ii, A_it; A_ti, A_tr+S_i] \ [f_i; f_tr + y_i];
        t_new{i} = u_i(N_i+1:end);
        u_new(ind_i) = u_i(1:N_i);
        u_new(ind_tr)= u_new(ind_tr) + t_new{i}/(n-1);
    end
    t_master = t_new;
    u_master = u_new;
    u_new = zeros(N,1);
    r = f_tr - A(ind_tr,:)*u_master;% residual on the trace
    iter = iter+1;

    figure(1)
    title(['Residual: ',num2str(norm(r))])
    surf(reshape(u_master,ny,nx))
    pause(0.1)
    err_out(iter)=norm(r);
end