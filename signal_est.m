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

%simulate state development
x=zeros(2,K);
F=zeros(2,2,K);
G=zeros(2,1,K);
x_0 = [0;0];
x_e = [10;200];
c_e = [-0.1182 0.5726; -1.0506 0.1111]*x_e;
x_k = x_0;
i = 1;
for k=1:K
    tk = ts*k;
    %default system
     F(:,:,k) = expm(A{sigma(i)}*ts);
     G(:,:,k) = inv(A{sigma(i)})*(F(:,:,k)-eye(2))*B;
     %c_e'*(F(:,:,k)*x_k+G(:,:,k)*u-x_e)
    
    %check if switch needed and if so when
% $$$     if sigma(i)==1
% $$$         if c_e'*(F(:,:,k)*x_k+G(:,:,k)*u-x_e)>0
% $$$             syms 't_switch';
% $$$             A_switch = expm(A{sigma(i)}*t_switch);
% $$$             B_switch = inv(A{sigma(i)})*(A_switch-eye(2))*B;
% $$$             x_t_switch = A_switch*x_k + B_switch*u;
% $$$             eqn=c_e'*(x_t_switch-x_e)==0;
% $$$             sol = solve(eqn,t_switch,'Real',true);
% $$$             t_sols = eval(sol);
% $$$             t_sol = min(t_sols(t_sols>=0));
% $$$             switch_times(end+1) = tk-ts+t_sol;
% $$$             sigma(i+1) = 2;
% $$$         end
% $$$     else
% $$$         if c_e'*(F(:,:,k)*x_k+G(:,:,k)*u-x_e)<=0
% $$$             syms 't_switch';
% $$$             A_switch = expm(A{sigma(i)}*t_switch);
% $$$             B_switch = inv(A{sigma(i)})*(A_switch-eye(2))*B;
% $$$             x_t_switch = A_switch*x_k + B_switch*u;
% $$$             eqn=c_e'*(x_t_switch-x_e)==0;
% $$$             sol = solve(eqn,t_switch,'Real',true);
% $$$             t_sols = eval(sol);
% $$$             t_sol = min(t_sols(t_sols>=0));
% $$$             switch_times(end+1) = tk-ts+t_sol;
% $$$             sigma(i+1) = 1;
% $$$         end
% $$$     end
    
    curr_switches=switch_times(switch_times>tk-ts & switch_times<=tk);
    if ~isempty(curr_switches);
        sigma(i+1)=3-sigma(i);
        F(:,:,k) = c2d_switch(curr_switches-(tk-ts),A{sigma(i)},A{sigma(i+1)});
        G(:,:,k) = input_switch(curr_switches-(tk-ts),A{sigma(i)},A{sigma(i+1)},B);
        i=i+1;
    else
        F(:,:,k) = expm(A{sigma(i)}*ts);
        G(:,:,k) = inv(A{sigma(i)})*(F(:,:,k)-eye(2))*B;
    end
    x_k = reshape(F(:,:,k),[2 2])*x_k+reshape(G(:,:,k),[2 1])*u;
    x(:,k)=x_k;
end


%Monte Carlo test of algorithm
N=10;
t_est=1*ts;
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

SE = zeros(N,2,K);
SE_kf = zeros(N,2,K);

for n=1:N
    x_hat_k = x_0+init_std*randn(2,1);
    x_hat_kf_k = x_hat_k;
    P_k = P_init;
    P_kf_k = P_init;
    
    x_hat = zeros(2,K);
    x_hat_kf = zeros(2,K);
    
    i=1;
    for k=1:K
        tk=ts*k;
        F_true = reshape(F(:,:,k),[2 2]);
        G_true = reshape(G(:,:,k),[2 1]);
        
        curr_switches=switch_times(switch_times>tk-ts & switch_times<=tk);
        if ~isempty(curr_switches);
            t_est = curr_switches-(tk-ts)+ts/100;
            F_hat = c2d_switch(t_est,A{sigma(i)},A{sigma(i+1)});
            G_hat = input_switch(t_est,A{sigma(i)},A{sigma(i+1)},B);
            i=i+1;
            tk;
        else
            F_hat = F_true;%reshape(F(:,:,k),[2 2]);
            G_hat = G_true;%reshape(G(:,:,k),[2 1]);
        end
        
        P_k = F_hat*P_k*F_hat';
        P_kf_k = F_true*P_kf_k*F_true';
        
        K_k = P_k*H'*inv(H*P_k*H'+R);
        K_kf_k = P_kf_k*H'*inv(H*P_kf_k*H'+R);
        
        x_hat_k = F_hat*x_hat_k + G_hat*u;
        x_hat_kf_k = F_true*x_hat_kf_k + G_true*u;
        
        
        y_k = H*x(:,k)+sqrt(R)*randn(1);
        
        x_hat_k = x_hat_k + K_k*(y_k-H*x_hat_k);
        x_hat_kf_k = x_hat_kf_k + K_kf_k*(y_k-H*x_hat_kf_k);
        
        x_hat(:,k) = x_hat_k;
        x_hat_kf(:,k) = x_hat_kf_k;
        
        P_k = (eye(2)-K_k*H)*P_k;
        P_kf_k = (eye(2)-K_kf_k*H)*P_kf_k;
        
    end
        
    error = x_hat-x;
    SE(n,:,:) = error.^2;
    
    error_kf = x_hat_kf-x;
    SE_kf(n,:,:) = error_kf.^2;
end

RMSE = reshape(sqrt(sum(SE,1)/N),[2,K]);
RMSE_kf = reshape(sqrt(sum(SE_kf,1)/N),[2,K]);

figure
plot(ts:ts:T,x(1,:))
hold on
plot(ts:ts:T,x(2,:))
hold off
h=legend('$i_L$ [amps]','$v_0$ [volts]');
set(h,'interpreter','latex');
xlabel('time [s]')
ylabel('magnitude')

figure
%xlim([0,3]);
plot(switch_times,zeros(size(switch_times)),'kx')
hold on
plot(ts:ts:T,RMSE(1,:))
plot(ts:ts:T,RMSE_kf(1,:))
ylim([0 15])
xlabel('time [s]')
ylabel('RMSE [amps]')
h=legend('switch times','RMSE, $\bar{t}_{est}$','RMSE known $\bar{t}$');
set(h,'interpreter','latex');
title('RMS Estimation Error, $i_L$','interpreter','latex')

figure
%xlim([0,3]);
plot(switch_times,zeros(size(switch_times)),'kx')
hold on
plot(ts:ts:T,RMSE(2,:))
plot(ts:ts:T,RMSE_kf(2,:))
ylim([0 15])
xlabel('time [s]')
ylabel('RMSE [volts]')
h=legend('switch times','RMSE, $\bar{t}_{est}$','RMSE known $\bar{t}$');
set(h,'interpreter','latex');
title('RMS Estimation Error, $v_0$','interpreter','latex')
