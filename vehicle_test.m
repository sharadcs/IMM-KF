close all
clear all

A1 = [0 0.1 0 0; 0.01 0 0 0; 0 0 0 0.1; 0 0 0.01 0];
A2 = [0 0.1 0 0; 0.01 0.08 0 0.02; 0 0 0 0.1; 0 0.02 0.01 0.08];
A3 = [0 0.1 0 0; 0.01 0.08 0 0.02; 0 0 0 0.1; 0 0.02 0.01 0.08];
B1 = [0 0.01 0 0.01]';
B2 = [0 10 0 -10]';
B3 = [0 -10 0 10]';
A = {A1,A2,A3};
B = {B1,B2,B3};
num_modes = 3;
N=4;

seed = 195;
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

T = 10;
tau = 0.5; % dwell time constraint
ts = 0.01; % simulation sampling period
K = T/ts;


for trials=1:1
% random switching signal
switch_times = T*rand(1,10);
switch_times = sort(switch_times);
i = 1;
while i<length(switch_times)
    if switch_times(i+1)-switch_times(i)<tau
        if i+1<length(switch_times)
            switch_times = [switch_times(1:i) switch_times(i+2:end)];
        else
            switch_times = switch_times(1:i);
        end
    else
        i=i+1;
    end
end

switch_pattern = ceil(num_modes*rand(length(switch_times)+1,1));
for i=1:length(switch_pattern)-1
    if(switch_pattern(i)==switch_pattern(i+1))
        switch_pattern(i+1)=mod(switch_pattern(i+1)+1,num_modes)+1;
    end
end

% simulate trajectory
x=zeros(4,K);
y=zeros(1,K);
sigma_k = zeros(1,K);
F=zeros(4,4,K);
G=zeros(4,1,K);
x_0 = [0;10;0;10];
x_k = [x_0];
u=1;

t_list = [0];
sigma_k = [1];
i = 1;
sigma = switch_pattern;



for k=1:K
    tk = ts*k;
    t_list = [t_list tk];
    %default system
     F(:,:,k) = expm(A{sigma(i)}*ts);
     G(:,:,k) = inv(A{sigma(i)})*(F(:,:,k)-eye(N))*B{sigma(i)};
     %c_e'*(F(:,:,k)*x_k+G(:,:,k)*u-x_e)
    
    curr_switches=switch_times(switch_times>tk-ts & switch_times<=tk);
    if ~isempty(curr_switches);
        %sigma(i+1)=3-sigma(i);
        F(:,:,k) = c2d_switch(curr_switches-(tk-ts),ts,A{sigma(i)},A{sigma(i+1)});
        G(:,:,k) = input_switch(curr_switches-(tk-ts),ts,A{sigma(i)},A{sigma(i+1)},B{sigma(i)},B{sigma(i+1)});
        i=i+1;
    else
        F(:,:,k) = expm(A{sigma(i)}*ts);
        G(:,:,k) = inv(A{sigma(i)})*(F(:,:,k)-eye(N))*B{sigma(i)};
    end
    x_k = reshape(F(:,:,k),[N N])*x_k+reshape(G(:,:,k),[N 1])*u;
    x(:,k)=x_k;
    %y(:,k)=H*x(:,k)+sqrt(R)*randn(1);
    sigma_k = [sigma_k sigma(i)];
end
% plot trajectory
fig1 = figure
plot(x(1,:),x(3,:))%,'filled')
hold on
%quiver(x(1,:),x(3,:),x(2,:),x(4,:))

down_sample = 10;
tm = down_sample*ts;
down_ind = down_sample:down_sample:length(t_list);
tm_list = t_list(down_ind);
xm = x(:,down_ind);
meas_R = 0.01;
%ym = xm([1,3],:) + sqrt(meas_R)*randn(2,length(down_ind));
ym = xm + sqrt(meas_R)*randn(4,length(down_ind));
plot(ym(1,:),ym(3,:));
hold off

% IMM-EV Kalman Filter implementation
viterbi_KF = {};
viterbi_KF.avg_error = zeros(size(xm));
viterbi_KF.MSE = zeros(size(xm,1),size(xm,1),size(xm,2));
viterbi_KF.sigma_errors = zeros(size(sigma_k));

window_len = 2;
W = window_len;
meas = {};
meas.y = [0;0];
meas.H = eye(4);%[1 0 0 0; 0 0 1 0];
meas.R = meas_R*eye(4);
meas.ts = tm;

dyn = {};
dyn.A = A;
dyn.B = B;
dyn.u = u;
n_states = 3;

viterbi_bank = {};
viterbi_bank.x{i} = [0;0;0;0];
viterbi_bank.P{i} = 100*eye(4);

viterbi_use_constraints = false;

for k = 1:length(tm_list)
    if k==1
        if ~viterbi_use_constraints
            for i=1:n_states
                viterbi_bank.x{i} = xm(:,k);
                viterbi_bank.P{i} = 100*eye(4);
            end
        else
            for i=1:n_states*n_states
                viterbi_bank.x{i} = xm(:,k)
                viterbi_bank.P{i} = eye(4);
            end
        end
    end
    if k<=W
        if ~viterbi_use_constraints
            for i=1:n_states
                viterbi_bank.x{i} = xm(:,k);
                viterbi_bank.P{i} = 100*eye(4);
                map_state=1;
            end
        else
            for i=1:n_states*n_states
                viterbi_bank.x{i} = xm(:,k);
                viterbi_bank.P{i} = eye(4);
                map_state=1;
            end
        end
    else
        if ~viterbi_use_constraints
            for i=1:n_states
                viterbi_init.x{i} = viterbi_KF.x(:,k-W);%Truth_KF.x(:,k-W);%-W+1)% %viterbi_bank.x{i};%
                viterbi_init.P{i} = viterbi_KF.P(:,:,k-W)+0.001*eye(4);%Truth_KF.P(:,:,k-W);%-W+1);%100* %5*viterbi_bank.P{i};%
            end
        else
            for i=1:n_states*n_states
                viterbi_init.x{i} = viterbi_KF.x(:,k-W);%Truth_KF.x(:,k);%-W+1);%
                viterbi_init.P{i} = 5*viterbi_KF.P(:,:,k-W);%Truth_KF.P(:,:,k);%-W+1);%
            end
        end
        
        y_last_W = ym(:,k-W+1:k);
        meas.y = y_last_W;
        
        output_viterbi = viterbi_kalman2(viterbi_init.x,viterbi_init.P,4,meas,window_len,dyn,n_states,viterbi_use_constraints);
        
        viterbi_bank.x=output_viterbi.x;
        viterbi_bank.P=output_viterbi.P;
        map_state=output_viterbi.map_state;      
    end
    %viterbi_bank.x{map_state}
    %viterbi_bank.P{map_state}
    viterbi_KF.x(:,k+1) = viterbi_bank.x{map_state};
    viterbi_KF.P(:,:,k+1) = viterbi_bank.P{map_state};
    viterbi_KF.sigma(:,k+1)=map_state;
    error = viterbi_KF.x(:,k+1)-xm(:,k);
    viterbi_KF_SE(:,:,k) = error*error';
    
end

figure(fig1)
hold on
plot(viterbi_KF.x(1,:),viterbi_KF.x(3,:))
hold off

end

