function [obj,switch_list] = generateMHE(window_len,N,M,n_modes,Ts,dwell_time)

% set time steps for optimization
max_switches = floor((window_len*Ts)/dwell_time);
Tvariable tm window_len;
Tvariable ts max_switches;
Tvariable switch_intervals max_switches; 

% enumerate possible switching signals
num_switch_signals = factorial(n_modes)/factorial(n_modes-1-max_switches)
switch_bank = {};
for i = 1:max_switches+1
    switch_bank{i} = nchoosek(1:n_modes,i);
end

obj_list = {};
obj_ind = {};

Tvariable active_modes window_len;
Tvariable x [N,window_len+max_switches];

Tvariable y [M,window_len];
Tvariable H [N,M];
Tvariable R [M,M];

% measurement cost function (max likelihood)
pred_meas = tprod(H,[-1,1,2],x,[-1,1,2]);
cost = tprod(pred_meas-y,);

% State and Input transition matrices


% dynamics contstraints
Tvariable A [N,N,window_len+max_switches];
Tvariable B [N,window_len+max_switches];
Tvariable u window_len+max_switches;


dyn_constr = 

% constraints on switch times
Tvariable tmax max_switches;
Tvariable tmin max_switches;

tmax_constr = tm <= tmax;
tmin_constr = tm >= tmin;

end