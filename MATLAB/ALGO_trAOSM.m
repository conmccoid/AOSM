function u_master = ALGO_trAOSM(A,f,sub_indices, u0)
% trAOSM runs a multi-subdomain AOSM that shares the trace with all
% subdomains
%   - A: system matrix
%   - f: right hand side
%   - sub_indices: array (?) of subdomain indices, the last of which is the
%   trace
%   - u0: initial guess on the trace

n = length(sub_indices);        % number of subdomains (+1 for the trace)
N = length(A);                  % size of global problem
itermax = N;                    % maximum number of iterations (nb: make input)
tol = 1e-12;                     % tolerance for KSP and Schwarz (nb: make input)
ind_tr = sub_indices{n};        % indices for the trace
M = length(ind_tr);             % size of trace
A_tr = A(ind_tr,ind_tr);        % matrix block for trace
f_tr = f(ind_tr);               % right hand side for trace
T_master = cell(n,1);           % master list of adapted transmission conditions
t_master = cell(n,1);           % master list for soln on trace
S_master = zeros(n*M,M);            % master list of sums of ATCs
u_master = zeros(N,1);          % master soln of global problem

for i=1:n-1                     % for each subdomain
    ind_i = sub_indices{i};     % indices for this subdomain
    N_i = length(ind_i);        % size of subdomain
    A_ii = A(ind_i,ind_i);      % main block for this subdomain
    A_it = A(ind_i,ind_tr);     % top right block
    A_ti = A(ind_tr,ind_i);     % bottom left block
    f_i = f(ind_i);             % right hand side for this subdomain
    T_i = zeros(M);             % init for adapted transmission conditions

    W_it = zeros(M);            % init W in trace
    W_ii = zeros(N_i,M);        % init W in subdomain

    u_it = u0;                  % init guess
    u_ii = -A_ii \ A_it * u_it; % init soln in subdomain
    u_i = [A_ii, A_it; A_ti, A_tr + T_i] \ [f_i; f_tr - A_ti*u_ii + T_i*u_it];
    d_ii = u_i(1:N_i) - u_ii;   % init d in subdomain
    d_it = u_i(N_i+1:end) - u_it; % init d in trace

    a = norm(d_it);             % normalization
    W_it(:,1) = d_it/a;         % first w vector on trace
    W_ii(:,1) = d_ii/a;         % first w vector in subdomain
    v_i = -A_ti*W_ii(:,1) + T_i*W_it(:,1); % first v vector
    T_i = T_i - v_i*W_it(:,1)'; % update to transmission conditions

    iter = 1;                   % init iteration count
    while a>tol && iter<itermax
        d_i = [A_ii, A_it; A_ti, A_tr + T_i] \ [sparse(N_i,1);a*v_i];
        u_i = u_i + d_i;        % solve and soln update
        %---MGS---%
        W_ii(:,iter) = d_i(1:N_i);
        W_it(:,iter) = d_i(N_i+1:end);
        for k=1:iter
            h_i = W_it(:,k)'*W_it(:,iter);
            W_it(:,iter) = W_it(:,iter) - h_i*W_it(:,k);
            W_ii(:,iter) = W_ii(:,iter) - h_i*W_ii(:,k);
        end
        a = norm(W_it(:,iter));
        W_it(:,iter) = W_it(:,iter)/a;
        W_ii(:,iter) = W_ii(:,iter)/a;
        v_i = -A_ti*W_ii(:,iter) + T_i*W_it(:,iter);
        T_i = T_i - v_i*W_it(:,iter)';
        iter = iter+1;
    end
    T_master{i} = T_i; % store adapted transmission conditions (nb: turn this list into cells for easier access)
    t_master{i} = u_i(N_i+1:end);
    for j=setdiff(1:n-1,i)
        S_master((1:M)+M*(j-1),1:M) = S_master((1:M)+M*(j-1),1:M) + T_i;
    end
    u_master(ind_i) = u_i(1:N_i);
    u_master(ind_tr)= u_master(ind_tr) + t_master{i}/2;

    figure(1)
    surf(reshape(u_master,sqrt(N),sqrt(N)))

    figure(2)
    contourf(T_i)
    pause
end

%===Global iterations===%
r = f_tr - A(ind_tr,:)*u_master;% residual on the trace
u_new = zeros(N,1);               % storage of new soln
t_new = cell(n,1);
iter = 1;                       % init interation count
while norm(r)>tol && iter<itermax
    for i=1:n-1
        ind_i = sub_indices{i};     % indices for this subdomain
        N_i = length(ind_i);        % size of subdomain
        A_ii = A(ind_i,ind_i);      % main block for this subdomain
        A_it = A(ind_i,ind_tr);     % top right block
        A_ti = A(ind_tr,ind_i);     % bottom left block
        f_i = f(ind_i);             % right hand side for this subdomain
%         S_i = S_master((1:M)+M*(i-1),1:M);
        S_i = zeros(M);
        y_i = zeros(M,1);
        for j=setdiff(1:n-1,i)
            ind_j = sub_indices{j};
            T_i = T_master{j};
            t_i = t_master{j};
            S_i = S_i + T_i;
            y_i = y_i - A(ind_tr,ind_j)*u_master(ind_j) + T_i*t_i;
        end
        u_i = [A_ii, A_it; A_ti, A_tr+S_i] \ [f_i; f_tr + y_i];
        t_new{i} = u_i(N_i+1:end);
        u_new(ind_i) = u_i(1:N_i);
        u_new(ind_tr)= u_new(ind_tr) + t_new{i}/2;
    end
    t_master = t_new;
    u_master = u_new;
    u_new = zeros(N,1);
    r = f_tr - A(ind_tr,:)*u_master;% residual on the trace
    iter = iter+1;

    figure(1)
    surf(reshape(u_master,sqrt(N),sqrt(N)))
    pause
end