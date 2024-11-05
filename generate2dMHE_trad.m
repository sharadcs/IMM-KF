function [Solvers,n_solvers] = generate2dMHE_trad(window_len,N,M,n_switches,Ts)
    
% set time steps for optimization
Tvariable t_bar [1];

% enumerate possible switching signals
max_switches = n_switches;
n_solvers = window_len; % TODO: set to correct number of solvers needed
%num_switch_signals = factorial(n_modes)/factorial(n_modes-1-max_switches)
%switch_bank = {};
%for i = 1:max_switches+1
%    switch_bank{i} = nchoosek(1:n_modes,i);
%end

%Tvariable active_modes window_len;
Tvariable x [N,window_len];
Tvariable xm [N,1];
x_aug = [x xm];

Tvariable A [N,N,n_switches+1];
Tvariable B [N,n_switches+1];

Tvariable y [M,window_len];
Tvariable beacon_pos [2,M];
%Tvariable H [M,N];
Tvariable Rinv [];

% measurement cost function (max likelihood)
% pred_meas = H*x;
% meas_diff = pred_meas - y;
% cost_meas = norm2(pred_meas-y);% + 100*norm2(pred_meas(:,1)-y(:,1));%sum(tprod(meas_diff,[1 -1],meas_diff,[1 -1]),1);
pp=reshape(x([1,3],:),[2,1,window_len]);
pp=pp(:,ones(M,1),:);
bb_pos=reshape(beacon_pos,[2,M,1]);
bb_pos=bb_pos(:,:,ones(window_len,1));
rel_pos=pp-bb_pos; %[3,nBeacons,nT]
distance=sqrt(tprod(rel_pos,[-1,1,2],rel_pos,[-1,1,2])); %[nBeacons,nT]
cost_meas = Rinv*norm2(y-distance);

% State and Input transition matrices
F1 = [1 Ts 0 0; 0 1 0 0; 0 0 1 Ts; 0 0 0 1];
F2 = [1 Ts 0 0; 0 1 0 0; 0 0 1 Ts; 0 0 0 1];
G1 = (Ts*eye(N) + Ts^2*reshape(A(:,:,1),[N,N])/2)*reshape(B(:,1),[N,1]);
G2 = (Ts*eye(N) + Ts^2*reshape(A(:,:,2),[N,N])/2)*reshape(B(:,2),[N,1]);

Fs1 = [1 t_bar 0 0; 0 1 0 0; 0 0 1 t_bar; 0 0 0 1];
Fs2 = [1 Ts-t_bar 0 0; 0 1 0 0; 0 0 1 Ts-t_bar; 0 0 0 1];
Gs1 = (reshape(t_bar,[])*eye(N) + reshape(t_bar,[])*reshape(t_bar,[])*reshape(A(:,:,1),[N,N])/2)*reshape(B(:,1),[N,1]);
Gs2 = (reshape((Ts-t_bar),[])*eye(N) + reshape((Ts-t_bar),[])*reshape((Ts-t_bar),[])*reshape(A(:,:,2),[N,N])/2)*reshape(B(:,2),[N,1]);

% dynamics contstraints
Tvariable u [1,window_len-1];
Tvariable d [N,window_len-1];
%Tvariable dm [N,max_switches];
%Tvariable active_const [3*window_len];

cost_proc = 100*(norm2(d));%+norm2(dm));

constraints = {};
n_solvers
window_len
for solv_num = 1:n_solvers
curr_constraints = {};

% constraints on dynamics
for i = 1:window_len-1
    %constraints{4*i-3} = reshape(active_const(3*i-2),[])*(x(:,i+1) == F1*x(:,i) + reshape(G1*u(:,i),[N,1]) + d(:,i));
    %constraints{4*i-2} = reshape(active_const(3*i-1),[])*(x(:,i+1) == F2*x(:,i) + reshape(G2*u(:,i),[N,1]) + d(:,i));
    %constraints{4*i-1} = reshape(active_const(3*i),[])*(xm == Fs1*x(:,i) + reshape(Gs1*u(:,i),[N,1]) + d(:,window_len));
    %constraints{4*i}   = reshape(active_const(3*i),[])*(x(:,i+1) == Fs2*xm + reshape(Gs2*u(:,i),[N,1]) + d(:,window_len+1));
    
    if i == solv_num
        curr_constraints{end+1} = (xm == Fs1*x(:,i) + reshape(Gs1*u(:,i),[N,1]))% + d(:,i));
        curr_constraints{end+1}   = (x(:,i+1) == Fs2*xm + reshape(Gs2*u(:,i),[N,1]))% + dm);%d(:,window_len));
    elseif i < solv_num
        curr_constraints{end+1} = (x(:,i+1) == F1*x(:,i) + reshape(G1*u(:,i),[N,1]))% + d(:,i));
    else
        curr_constraints{end+1} = (x(:,i+1) == F2*x(:,i) + reshape(G2*u(:,i),[N,1]))% + d(:,i));
    end
end

% constraints on switch times
tmax_constr = t_bar <= Ts;
tmin_constr = t_bar >= 0;

if solv_num<window_len
    curr_constraints{end+1} = tmax_constr;
    curr_constraints{end+1} = tmin_constr;
else
    %curr_constraints{end+1} = tmax_constr;
end

constraints{end+1} = curr_constraints;
end



% constraints on state
Tvariable x_max [N,window_len+max_switches];
Tvariable x_min [N,window_len+max_switches];
x_max_constr = x_aug <= x_max;
x_min_constr = x_aug >= x_min;
         

%constraints{end+1} = x_max_constr;
%constraints{end+1} = x_min_constr;
         
%Tvariable Hess_ [97,97];%[121,121]%(2*N*(window_len+max_switches) + 1 + 2*N*(window_len+max_switches) +2)*[1,1];
%Tvariable dHess_ size(Hess_,1);  % smaller

cost = cost_meas+cost_proc;

Solvers = {};
for solv_num = 1:n_solvers
this_constraints = constraints{solv_num};
name = sprintf('tmpTradFilter_%d',solv_num);
if solv_num<n_solvers
    classname=class2optimizeCS('classname',name,...
                          'objective',cost,...
                          'optimizationVariables',{...
                              x,xm,t_bar,d},...
                          'parameters',{...
                              y,Rinv,A,B,u,beacon_pos},...
                          'constraints',this_constraints,...
                          'outputExpressions',{...
                              cost,x,xm,t_bar},...
                          'addEye2Hessian',false,...
                          'solverVerboseLevel',1);
else
    classname=class2optimizeCS('classname',name,...
                          'objective',cost,...
                          'optimizationVariables',{...
                              x,d},...
                          'parameters',{...
                              y,Rinv,A,B,u,xm,t_bar,beacon_pos},...
                          'constraints',this_constraints,...
                          'outputExpressions',{...
                              cost,x,d},...
                          'addEye2Hessian',false,...
                          'solverVerboseLevel',1);
end
this_solver = feval(classname);
Solvers{solv_num} = this_solver;
end

end

