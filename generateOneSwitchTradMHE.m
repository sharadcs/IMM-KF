function [obj,label_list] = generateOneSwitchTradMHE(window_len,N,M,interval,A,B,Ts)

% set time steps for optimization
Tvariable t_bar [1];

% enumerate possible switching signals
max_switches = 1;
%num_switch_signals = factorial(n_modes)/factorial(n_modes-1-max_switches)
%switch_bank = {};
%for i = 1:max_switches+1
%    switch_bank{i} = nchoosek(1:n_modes,i);
%end

%Tvariable active_modes window_len;
Tvariable x [N,window_len];
Tvariable xm [N,1];
x_aug = [x(:,1:interval) xm x(:,interval+1:end)];

Tvariable y [M,window_len];
Tvariable H [M,N];
Tvariable R [M,M];

% measurement cost function (max likelihood)
pred_meas = H*x;
meas_diff = pred_meas - y;
cost_meas = 100*norm2(pred_meas-y)%sum(tprod(meas_diff,[1 -1],meas_diff,[1 -1]),1);

% State and Input transition matrices
F1 = expm(A{1}*Ts);
F2 = expm(A{2}*Ts);
G1 = (Ts*eye(N) + Ts^2*A{1}/2)*B{1};
G2 = (Ts*eye(N) + Ts^2*A{2}/2)*B{2};

Fs1 = [1 t_bar; 0 1];
Gs1 = [t_bar*t_bar/2; t_bar];
Fs2 = [1 Ts-t_bar; 0 1];
Gs2 = [(Ts-t_bar)*(Ts-t_bar); 2*Ts-t_bar];

% dynamics contstraints
Tvariable u [1,window_len+max_switches];
Tvariable d [N,window_len+max_switches];

cost_proc = 100*norm2(d);

constraints = {}
for i = 1:interval-1
    constraints{3(i-1)+1} = x(:,i+1) == F1*x(:,i) + reshape(G1*u(:,i),[2,1]) + d(:,i);
    constraints{3(i-1)+2} = xm == Fs1*x(:,interval)+reshape(Gs1*u(:, interval),[2,1]) + d(:,end-1);
    constraints{3(i-1)+2} = x(:,interval+1) == Fs2*xm + reshape(Gs2*u(:,interval),[2,1]) + d(:,interval+1);
end

constraints{interval+1} = x(:,interval+1) == Fs2*xm + reshape(Gs2*u(:,interval),[2,1]) + d(:,interval+1);
for i = interval+1:window_len-1
    constraints{i+1} = x(:,i+1) == Fs2*x(:,i) + reshape(Gs2*u(:,i),[2,1]) + d(:,i+1);
end

% constraints on switch times
tmax_constr = t_bar <= Ts;
tmin_constr = t_bar >= 0;

% constraints on state
Tvariable x_max [2,window_len+max_switches];
Tvariable x_min [2,window_len+max_switches];
x_max_constr = x_aug <= x_max;
x_min_constr = x_aug >= x_min;
         
constraints{end+1} = tmax_constr;
constraints{end+1} = tmin_constr;
%constraints{end+1} = x_max_constr;
%constraints{end+1} = x_min_constr;
         
%Tvariable Hess_ [97,97];%[121,121]%(2*N*(window_len+max_switches) + 1 + 2*N*(window_len+max_switches) +2)*[1,1];
%Tvariable dHess_ size(Hess_,1);  % smaller

cost = cost_meas+cost_proc;
         
classname=class2optimizeCS('classname','tmpTradFilt',...
                          'objective',cost,...
                          'optimizationVariables',{...
                              x,xm,t_bar,d},...
                          'parameters',{...
                              y,H,R,u,x_max,x_min},...
                          'constraints',constraints,...
                          'outputExpressions',{...
                              cost,x,xm,t_bar},...
                          'addEye2Hessian',false,...
                          'solverVerboseLevel',3);
obj = feval(classname);        
         
end