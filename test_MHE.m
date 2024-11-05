clear all
close all

regenerate_solvers = true; % false to use previously generated solver
if regenerate_solvers
    delete('toremove.m','tmp*');rc=rmdir('@tmp*','s');
end

A1 = [0 1; 0 0];
A2 = [0 1; 0 0];
B1 = [0; 1];
B2 = [0; 2];

sys1 = ss(A1,B1,eye(2),zeros(2,1));
sys2 = ss(A2,B2,eye(2),zeros(2,1));

Tsim = 0.01;
t_end = 1;
t = 0:Tsim:t_end;

t_switch = 0.45;
ind_switch = t_switch/Tsim;

t_1 = 0:Tsim:t_switch;
t_2 = 0:Tsim:t_end-t_switch;
x_1 = lsim(sys1,ones(size(t_1)),t_1);
x_2 = lsim(sys2,ones(size(t_2)),t_2,x_1(end,:));
x = [x_1; x_2(2:end,:)];

downsample = 10;
tk = 0:10*Tsim:t_end;
xk = x(1:10:end,:);
yk = xk(:,1) + 0.05*randn(size(xk(:,1)));

mhe = generateOneSwitchMHE(11,2,1,5);

setP_Ts(mhe,downsample*Tsim);
setP_A1(mhe,A1);
setP_A2(mhe,A2);
setP_B1(mhe,B1);
setP_B2(mhe,B2);
setP_y(mhe,yk');
setP_H(mhe,[1 0]);
setP_R(mhe,1);
setP_u(mhe,ones(1,size(yk,1)+1))
setP_x_max(mhe,10*ones(size(xk,2),size(xk,1)+1));
setP_x_min(mhe,-10*ones(size(xk,2),size(xk,1)+1));

setV_x(mhe,xk');
setV_xm(mhe,x(ind_switch,:)');
setV_t_bar(mhe,0.05);
setV_d(mhe,zeros(2,12));

mu0=1; % for IPM
maxIter=200; % external stopping condition
saveIter=-1;
[status,iter,runtime]=solve(mhe,mu0,int32(maxIter),int32(saveIter)); % runs solver

[J_hat,x_hat,xm_hat,t_bar_hat] = getOutputs(mhe); % returns
                                                  % specified
                                                  % outputs

trad_mhe = generateOneSwitchTradMHE(11,2,1,5,{A1,A2},{B1,B2},downsample*Tsim)

setP_y(trad_mhe,yk');
setP_H(trad_mhe,[1 0]);
setP_R(trad_mhe,1);
setP_u(trad_mhe,ones(1,size(yk,1)+1))
setP_x_max(trad_mhe,10*ones(size(xk,2),size(xk,1)+1));
setP_x_min(trad_mhe,-10*ones(size(xk,2),size(xk,1)+1));

setV_x(trad_mhe,xk');
setV_xm(trad_mhe,x(ind_switch,:)');
setV_t_bar(trad_mhe,0.05);
setP_d(trad_mhe,zeros(2,12));

mu0=1; % for IPM
maxIter=200; % external stopping condition
saveIter=-1;
[status,iter,runtime]=solve(trad_mhe,mu0,int32(maxIter),int32(saveIter)); % runs solver

[J_hat_trad,x_hat_trad,xm_hat_trad,t_bar_hat_trad] = getOutputs(trad_mhe); % returns
                                                  % specified
                                                  % outputs

figure
plot(tk,xk(:,1))
hold on
plot(tk,x_hat(1,:))
plot(tk,x_hat_trad(1,:))