function [obj] = generateOneSwitchTradMHE(window_len,N,M,interval,A,B,Ts)

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


end