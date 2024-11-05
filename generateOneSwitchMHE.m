function [obj] = generateOneSwitchMHE(window_len,N,M,interval)

% set time steps for optimization
Tvariable Ts [1];
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
Tvariable xm [N];

Tvariable y [M,window_len];
Tvariable H [M,N];
Tvariable R [M,M];

% measurement cost function (max likelihood)
pred_meas = H*x;
meas_diff = pred_meas - y;
cost = sum(tprod(meas_diff,[1 -1],meas_diff,[1 -1]),1);

% State and Input transition matrices
Tvariable A1 [N,N];
Tvariable A2 [N,N];
Tvariable B1 [N,1];
Tvariable B2 [N,1];

% dynamics contstraints
Tvariable u [1,window_len+max_switches];
Tvariable d [N,window_len+max_switches];
cost = cost + 100*norm2(d);

x_aug = [x(:,1:interval) xm x(:,interval+1:end)];
dx = [x_aug(:,2)-x_aug(:,1) x_aug(:,3:end)-x_aug(:,1:end-2) x_aug(:,end)-x_aug(:,end-1)];
dt = reshape([Ts; repmat(Ts,[interval-2]); Ts+t_bar; Ts; 2*Ts-t_bar; ...
              repmat(Ts,[window_len-interval-2]); Ts],[1,window_len+max_switches]);

dxdt = [A1*x(:,1:interval) A2*xm A2*x(:,interval+1:end)] + ...
       [B1*u(:,1:interval) B2*u(:,interval+1:end)];

dyn_constraint = tprod(repmat(dt,[N 1]),[1,2],dxdt,[1,2]) == dx + d;

% constraints on switch times
tmax_constr = t_bar <= Ts;
tmin_constr = t_bar >= 0;

% constraints on state
Tvariable x_max [2,window_len+max_switches];
Tvariable x_min [2,window_len+max_switches];
x_max_constr = x_aug <= x_max;
x_min_constr = x_aug >= x_min;
         
constraints = {
    dyn_constraint,
    tmax_constr,
    tmin_constr,
    x_max_constr,
    x_min_constr };
         
%Tvariable Hess_ (10*nT+8   +   4*nT)*[1,1];  % smaller
%Tvariable dHess_ size(Hess_,1);  % smaller
         
classname=class2optimizeCS('classname','tmpFilt',...
                           'objective',cost,...
                          'optimizationVariables',{...
                              x,xm,t_bar,d},...
                          'parameters',{...
                              Ts,A1,A2,B1,B2,y,H,R,u,x_max,x_min},...
                          'constraints',constraints,...
                          'outputExpressions',{...
                              cost,x,xm,t_bar},...
                          'addEye2Hessian',false,...
                          'solverVerboseLevel',3);
         
obj = feval(classname);        
         
end