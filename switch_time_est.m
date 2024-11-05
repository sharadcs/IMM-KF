close all;
delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');

u = 1000;
R1 = 2;
L1 = 500e-6;
C0 = 470e-6;
R0 = 50;

%continuous time state dynamics
A1 = [-R1/L1 0; 0 -1/(R0*C0)];
A2 = [-R1/L1 -1/L1; 1/C0 -1/(R0*C0)];

B = [1/L1; 0];

A = {A1,A2};

sigma = [1,2,1,2,1];
switch_periods = 0.0021:0.0021:0.05;
switch_times = [switch_periods];% (switch_periods-0.001)];%[0.4444,2.1111,3.8888,5.6666];
T=0.05;
ts=0.0005;
K=T/ts;

init_std = 10;
P_init = init_std^2*eye(2);

%discretization of system
c2d_switch = @(t,A_beg,A_end) expm(A_end*(ts-t))*expm(A_beg*t);
input_switch = @(t,A_beg,A_end,B_b) inv(A_end)*(expm(A_end*(ts-t))-eye(2))*B_b ...
    + expm(A_end*(ts-t))*inv(A_beg)*(expm(A_beg*(t))-eye(2 ))*B_b;

%discrete time measurements
H = [1 0];%eye(2);
R = 50000;%*eye(2);

A_switch = c2d_switch(A1,A2,0.0001);

x_0 = [0,100]';
x_1 = A_switch*x_0;

y_0 = H*x_0 + sqrt(R)*randn(1);
y_1 = H*x_1 + sqrt(R)*randn(1);

% Solver setup

Tvariable x_beg [2,1];
%Tvariable x_end [2,1];
Tvariable y_meas_beg [1];
Tvariable y_meas_end [1];

Tvariable t_switch [];

y_beg = H*x_beg;
y_end = H*expm(A2*(ts-t_switch))*expm(A1*t_switch)*x_end;

J = reshape(R*norm2(y_beg-y_meas_beg) + R*norm2(y_end-y_meas_end),1);

classname=class2optimizeCS('classname',name,...
                         'objective',J,...
                         'optimizationVariables',{t_switch},...
                         'outputExpressions',{J,t_switch},...
                         'parameters',{x_beg,y_meas_beg,y_meas_end},...
                         'solverVerboseLevel',1);
% Create object
solver=feval(classname);

% Solve
solver.setP_x_beg(x_0);
solver.setP_y_meas_beg(y_0);
solver.setP_y_meas_end(y_1);

solver.setV_t_switch(0.0001);

mu0=1; % for IPM
maxIter=200; % external stopping condition
saveIter=-1;
[stat,iter,time]=solve(solver,mu0,int32(maxIter),int32(saveIter)); % runs solver

[J_hat,t_hat] = getOutputs(solver); % returns
                                                          % specified