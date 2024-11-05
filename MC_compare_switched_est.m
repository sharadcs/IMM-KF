close all;

u = 100;
R1 = 2;
L1 = 500e-6;
C0 = 470e-6;
R0 = 50;
n = 2;

%continuous time state dynamics
A1 = [-R1/L1 0; 0 -1/(R0*C0)];
A2 = [-R1/L1 -1/L1; 1/C0 -1/(R0*C0)];

B = [1/L1; 0];

A = {A1,A2};

T=0.05;%0.5;%
ts=0.0005;
K=T/ts;
sigma = [1,2,1,2,1];
switch_dwell = 0.0005;%0.0021;%
switch_periods = switch_dwell:switch_dwell:T-0.001;%0.002:0.002:T;%
switch_times = [switch_periods];% (switch_periods-0.001)];%[0.4444,2.1111,3.8888,5.6666];

init_std = 10;
P_init = init_std^2*eye(2);

%discretization of system
%c2d_switch = @(t,A_beg,A_end) expm(A_end*(ts-t))*expm(A_beg*t);
%input_switch = @(t,A_beg,A_end,B_b) inv(A_end)*(expm(A_end*(ts-t))-eye(2))*B_b ...
%    + expm(A_end*(ts-t))*inv(A_beg)*(expm(A_beg*(t))-eye(2))*B_b;

%discrete time measurements
H = [1 0];%eye(2);%[1 0];
R = 5;%50000;%*eye(2);

%simulated state development
x=zeros(2,K);
y=zeros(1,K);
sigma_k = zeros(1,K);
F=zeros(2,2,K);
G=zeros(2,1,K);
x_0 = [0;0];
x_e = [10;200];
c_e = [-0.1182 0.5726; -1.0506 0.1111]*x_e;
x_k = [x_0];



%% run MC trials of various filtering methods
num_MC_trials = 1;

% errors over trials
Truth_KF = {};
Truth_KF.avg_error = zeros(size(x));
Truth_KF.MSE = zeros(size(x,1),size(x,1),size(x,2));
Truth_KF.sigma_errors = zeros(size(sigma_k));
onestep_KF = {};
onestep_KF.avg_error = zeros(size(x));
onestep_KF.MSE = zeros(size(x,1),size(x,1),size(x,2));
onestep_KF.sigma_errors = zeros(size(sigma_k));
viterbi_KF = {};
viterbi_KF.avg_error = zeros(size(x));
viterbi_KF.MSE = zeros(size(x,1),size(x,1),size(x,2));
viterbi_KF.sigma_errors = zeros(size(sigma_k));
MI_MHE = {};
MI_MHE.avg_error = zeros(size(x));
MI_MHE.MSE = zeros(size(x,1),size(x,1),size(x,2));
MI_MHE.sigma_errors = zeros(size(sigma_k));

for trial = 1:num_MC_trials
    
first_switch = switch_dwell*rand();
switch_times = first_switch+switch_periods;

% generate trajectory and measurements
t_list = [0];
sigma_k = [1];
i = 1;
for k=1:K
    tk = ts*k;
    t_list = [t_list tk];
    %default system
     F(:,:,k) = expm(A{sigma(i)}*ts);
     G(:,:,k) = inv(A{sigma(i)})*(F(:,:,k)-eye(2))*B;
     %c_e'*(F(:,:,k)*x_k+G(:,:,k)*u-x_e)
    
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
    y(:,k)=H*x(:,k)+sqrt(R)*randn(1);
    sigma_k = [sigma_k sigma(i)];
end

% data for various filtering methods
Truth_KF.x = zeros(size(x));
Truth_KF.sigma = zeros(1,size(x,2));
Truth_Kf.P = zeros(n,n,size(x,2));
Truth_KF.t_c = zeros(1,size(x,2));
onestep_KF.x = zeros(size(x));
onestep_KF.sigma = zeros(1,size(x,2));
onestep_Kf.P = zeros(n,n,size(x,2));
onestep_KF.t_c = zeros(1,size(x,2));
viterbi_KF.x = zeros(size(x));
viterbi_KF.sigma = zeros(1,size(x,2));
viterbi_Kf.P = zeros(n,n,size(x,2));
viterbi_KF.t_c = zeros(1,size(x,2));
MI_MHE.x = zeros(size(x));
MI_MHE.sigma = zeros(1,size(x,2));
MI_MHE.P = zeros(n,n,size(x,2));
MI_MHE.t_c = zeros(1,size(x,2));

W = 4;  %window length for windowed methods
y_last_W = zeros(size(y,1),W); %last W measurements

% Initializations for Kalman filters
Truth_KF.x(:,1) = x_0;
Truth_KF.P(:,:,1) = 1000*eye(n);
Truth_KF.sigma(:,1) = sigma_k(1);
onestep_KF.x(:,1) = x_0;
onestep_KF.P(:,:,1) = 1000*eye(n);
onestep_KF.sigma(:,1) = sigma_k(1);
viterbi_KF.x(:,1) = x_0;
viterbi_KF.P(:,:,1) = 1000*eye(n);
viterbi_KF.sigma(:,1) = sigma_k(1);

% Squared error for computing MSE
Truth_KF_SE = zeros(n,n,size(x,2));
onestep_KF_SE = zeros(n,n,size(x,2));
viterbi_KF_SE = zeros(n,n,size(x,2));
MI_MHE_SE = zeros(n,n,size(x,2));

ndivs = 10; % divisions for estimating 
n_states = 2; % number of discrete states

% iterate through time, compute estimates
viterbi_bank.x = {};
viterbi_bank.P = {};
viterbi_use_constraints=false;
for k=1:size(y,2)
    if k>=W
        y_last_W = y(:,k-W+1:k);
    end
    y_true = y(:,k);
    
    % KF with true switch
    pred = F(:,:,k)*Truth_KF.x(:,k)+G(:,:,k)*u;
    pred_cov = F(:,:,k)*Truth_KF.P(:,:,k)*F(:,:,k)';
    y_pred = H*pred;
    y_pred_cov = H*pred_cov*H' + R;
    K_kf = pred_cov*H'*inv(y_pred_cov);
    Truth_KF.x(:,k+1) = pred + K_kf*(y_true-H*pred);
    Truth_KF.P(:,:,k+1) = (eye(2)-K_kf*H)*pred_cov;
    Truth_KF.sigma(:,k+1) = sigma_k(k+1);
    error = Truth_KF.x(:,k+1)-x(:,k);
    Truth_KF_SE(:,:,k) = error*error';
    
    % KF with 1-step optimal switch
    ml_innov = inf;
    mode1 = onestep_KF.sigma(:,k);
    for mode2 = 1:n_states
        if mode1 ~= mode2
            for t_bar = (ts/ndivs)/2:ts/ndivs:ts-(ts/ndivs)/2
                F_hat = c2d_switch(t_bar,A{mode1},A{mode2});
                G_hat = input_switch(t_bar,A{mode1},A{mode2},B);
                pred = F_hat*onestep_KF.x(:,k) + G_hat*u;%F_hat*onestep_KF.x(:,k) + G_hat*u;
                pred_cov = F_hat*onestep_KF.P(:,:,k)*F_hat' + 10*eye(2);%F_hat*onestep_KF.P(:,:,k)*F_hat';
                y_pred = H*pred;
                  y_pred_cov = H*pred_cov*H' + R;
                temp = exp(-(y_pred-y_true)'*y_pred_cov*(y_pred-y_true));
                if temp<ml_innov
                   ml_innov=temp;
                   K_kf = pred_cov*H'*inv(y_pred_cov);
                   onestep_KF.x(:,k+1) = pred + K_kf*(y_true-H*pred);
                   onestep_KF.P(:,:,k+1) = (eye(2)-K_kf*H)*pred_cov;
                   onestep_KF.sigma(:,k+1)=mode2;
                end
            end
        else
            F_hat = expm(A{mode1}*ts);
            G_hat = inv(A{mode1})*(expm(A{mode1}*ts)-eye(2))*B;
            pred = F_hat*onestep_KF.x(:,k) +  G_hat*u;%F_hat*onestep_KF.x(:,k) +  G_hat*u;
            pred_cov = F_hat*onestep_KF.P(:,:,k)*F_hat';%F_hat*onestep_KF.P(:,:,k)*F_hat';
            y_pred = H*pred;
            y_pred_cov = H*pred_cov*H' + R;
            temp = exp(-(y_pred-y_true)'*y_pred_cov*(y_pred-y_true));
            if temp<ml_innov
                ml_innov=temp;
                K_kf = pred_cov*H'*inv(y_pred_cov);
                onestep_KF.x(:,k+1) = pred + K_kf*(y_true-H*pred);
                onestep_KF.P(:,:,k+1) = (eye(2)-K_kf*H)*pred_cov;
                onestep_KF.sigma(:,k+1)=mode2;
             end
        end
    end
    error = onestep_KF.x(:,k+1)-x(:,k);
    onestep_KF_SE(:,:,k) = error*error';
    
    % KF with windowed viterbi
    if k==1
        if ~viterbi_use_constraints
            for i=1:n_states
                viterbi_bank.x{i} = onestep_KF.x(:,1);%Truth_KF.x(:,1);
                viterbi_bank.P{i} = onestep_KF.P(:,:,1);%Truth_KF.P(:,:,1);
            end
        else
            for i=1:n_states*n_states
                viterbi_bank.x{i} = onestep_KF.x(:,1);%Truth_KF.x(:,1);
                viterbi_bank.P{i} = onestep_KF.P(:,:,1);%Truth_KF.P(:,:,1);
            end
        end
    end
    if k<=W
        if ~viterbi_use_constraints
            for i=1:n_states
                viterbi_bank.x{i} = onestep_KF.x(:,k);%Truth_KF.x(:,k);%
                viterbi_bank.P{i} = onestep_KF.P(:,:,k);%Truth_KF.P(:,:,k);%
                map_state=1;
            end
        else
            for i=1:n_states*n_states
                viterbi_bank.x{i} = onestep_KF.x(:,k);%Truth_KF.x(:,k);%
                viterbi_bank.P{i} = onestep_KF.P(:,:,k);%Truth_KF.P(:,:,k);%
                map_state=1;
            end
        end
    else
        if ~viterbi_use_constraints
            for i=1:n_states
                viterbi_init.x{i} = viterbi_KF.x(:,k-W);%Truth_KF.x(:,k-W);%-W+1)% %viterbi_bank.x{i};%
                viterbi_init.P{i} = 100*viterbi_KF.P(:,:,k-W);%+1000*eye(2);%Truth_KF.P(:,:,k-W);%-W+1);%100* %5*viterbi_bank.P{i};%
            end
        else
            for i=1:n_states*n_states
                viterbi_init.x{i} = viterbi_KF.x(:,k-W);%Truth_KF.x(:,k);%-W+1);%
                viterbi_init.P{i} = 5*viterbi_KF.P(:,:,k-W);%Truth_KF.P(:,:,k);%-W+1);%
            end
        end
        
        meas.y = y_last_W;
        meas.H=H;
        meas.R=R;
        meas.ts=ts;
        dyn.A=A;
        dyn.B=B;
        dyn.u=u;
        
        output_viterbi = viterbi_kalman(viterbi_init.x,viterbi_init.P,meas,W,dyn,n_states,viterbi_use_constraints);
        
        viterbi_bank.x=output_viterbi.x;
        viterbi_bank.P=output_viterbi.P;
        map_state=output_viterbi.map_state;;      
    end
    viterbi_KF.x(:,k+1) = viterbi_bank.x{map_state};
    viterbi_KF.P(:,:,k+1) = viterbi_bank.P{map_state};
    viterbi_KF.sigma(:,k+1)=map_state;
    error = viterbi_KF.x(:,k+1)-x(:,k);
    viterbi_KF_SE(:,:,k) = error*error';
    
    
    % MHE for all possible switches
    
    
end

% compute statistics
Truth_KF_error = Truth_KF.x(:,2:end)-x;
Truth_KF.avg_error = ((trial-1)*Truth_KF.avg_error + Truth_KF_error)/trial;
Truth_KF_MSE = ((trial-1)*Truth_KF.MSE + Truth_KF_SE)/trial;
onestep_KF_error = onestep_KF.x(:,2:end)-x;
onestep_KF.avg_error = ((trial-1)*onestep_KF.avg_error + onestep_KF_error)/trial;
onestep_KF_MSE = ((trial-1)*onestep_KF.MSE + onestep_KF_SE)/trial;
viterbi_KF_error = viterbi_KF.x(:,2:end)-x;
viterbi_KF.avg_error = ((trial-1)*viterbi_KF.avg_error + viterbi_KF_error)/trial;
viterbi_KF_MSE = ((trial-1)*viterbi_KF.MSE + viterbi_KF_SE)/trial;
end

%% Plots

% plot of x1 estimation error vs time
figure
plot(t_list(2:end),Truth_KF.avg_error(1,:))
hold on
plot(t_list(2:end),onestep_KF.avg_error(1,:))
plot(t_list(2:end),viterbi_KF.avg_error(1,:))
legend('Truth KF','one step KF','Viterbi KF')
xlabel('time')
ylabel('error in $i_L$','interpreter','Latex')

% plot of x2 estimation error vs time
figure
plot(t_list(2:end),Truth_KF.avg_error(2,:))
hold on
plot(t_list(2:end),onestep_KF.avg_error(2,:))
plot(t_list(2:end),viterbi_KF.avg_error(2,:))
legend('Truth KF','one step KF','Viterbi KF')
xlabel('time')
ylabel('error in $v_0$','interpreter','Latex')
