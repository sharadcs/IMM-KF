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

T=0.05;
ts=0.0005;
K=T/ts;
sigma = [1,2,1,2,1];
switch_periods = 0.0011:0.0011:T;%0.002:0.002:T;%
switch_times = 0.0005+[switch_periods];% (switch_periods-0.001)];%[0.4444,2.1111,3.8888,5.6666];

init_std = 1000;
P_init = init_std^2*eye(2);

%discretization of system
c2d_switch = @(t,A_beg,A_end) expm(A_end*(ts-t))*expm(A_beg*t);
input_switch = @(t,A_beg,A_end,B_b) inv(A_end)*(expm(A_end*(ts-t))-eye(2))*B_b ...
    + expm(A_end*(ts-t))*inv(A_beg)*(expm(A_beg*(t))-eye(2))*B_b;

%discrete time measurements
H = [0 1];%eye(2);%[1 0];
R = 5;%50000;%*eye(2);

%simulate state development
x=zeros(2,K);
F=zeros(2,2,K);
G=zeros(2,1,K);
x_0 = [0;0];
x_e = [10;200];
c_e = [-0.1182 0.5726; -1.0506 0.1111]*x_e;
x_k = [x_0];
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
    sigma_k = [sigma_k sigma(i)];
end


% Viterbi algorithm to estimate switching signal with discretized
% time instants

W = 5; % window length in # of samples with measurements
ndivs = 10;  % # of divisions of sample to check paths
use_constraints=false;

% state definition:
% 1 - system 1, no switch last sample
% 2 - system 1, switched last sample
% 3 - system 2, no switch last sample
% 4 - system 2, switched last sample
if use_constraints
    n_states = 4;%2;%4; % number of discrete states
else
    n_states = 2;
end
    
% matrix of possible state transitions:
if use_constraints
    M_tr = [ 1 0 0 1; 1 0 0 0; 0 1 1 0; 0 0 1 0];%[1 1; 1 1];%[ 1 0 0 1; 1 0 0 0; 0 1 1 0; 0 0 1 0];
else
    M_tr = [1 1; 1 1];
end
    
% trellis defined l-r t-b of transitions
path_probs = zeros(size(M_tr,1)*size(M_tr,2),W); % transition probabilities
ml_transition = zeros(n_states,W+1); % best previous state for each state at each level
state_likelihood = zeros(n_states,W+1); % overall discrete state likelihoods
state_likelihood(:,1)=ones(n_states,1)/n_states;%[0.25 0.25 0.25 0.25]';%0.5*ones(n_states,1);
t_hats = zeros(W+1,1);
KF_bank = {};
KF_hist = {};
for i = 1:n_states
    KF_bank{i} = {};
    KF_bank{i}.x=x_0;
    KF_bank{i}.P=P_init;
end
KF_hist{1} = KF_bank;

% loop through levels of trellis and update probabilities
tic
delay = 10;
for w = 2:W
    x_true = x(:,w+delay);
    y_true = H*x_true+sqrt(R)*randn(1);
    w
    for i = 1:n_states
       x_hat_ml = KF_bank{i}.x;
       P_ml = KF_bank{i}.P;
       for j = 1:n_states
          if M_tr(i,j)
              if use_constraints
                  mode1 = ceil(i/2);%i;%ceil(i/2);
                  mode2 = ceil(j/2);%j;%ceil(j/2);
              else
                  mode1 = i;
                  mode2 = j;
              end
              w;
              j;
              if i~=j
                  'switch';
                  curr_t_l = 0;
                  pred_ml = 0;
                  t_bar_ml = nan;
                  maybe_x = x_hat_ml;
                  maybe_P = P_ml;
                  for t_bar = (ts/ndivs)/2:ts/ndivs:ts-(ts/ndivs)/2
                      F_hat = c2d_switch(t_bar,A{mode1},A{mode2});
                      G_hat = input_switch(t_bar,A{mode1},A{mode2},B);
                      pred = F_hat*KF_bank{mode1}.x + G_hat*u;
                      x(:,w);
                      pred_cov = F_hat*KF_bank{mode1}.P*F_hat';
                      y_pred = H*pred;
                      y_pred_cov = H*pred_cov*H' + R;
                      temp = exp(-(y_pred-y_true)'*inv(y_pred_cov)*(y_pred-y_true));%y_pred_cov
                      pred_ml = pred_ml+temp/ndivs;
                      if temp>=curr_t_l
                          curr_t_l = temp;
                          t_bar_ml = t_bar;
                          K_kf = pred_cov*H'*inv(y_pred_cov);
                          maybe_x = pred + K_kf*(y_true-H*pred);
                          maybe_P = (eye(2)-K_kf*H)*pred_cov;
                      end
                  end
                  t_bar_ml;
                  'switch';
                  pred_ml;
                  if state_likelihood(i,w-1)*pred_ml >= state_likelihood(j,w)
                      state_likelihood(j,w) = state_likelihood(i,w-1)*pred_ml;
                      ml_transition(j,w) = i;
                      x_hat_ml = maybe_x;
                      P_ml = maybe_P;
                  end
                  t_hats(w) = ts*(w-1)+t_bar_ml;
              else
                  'no switch';
                  F_hat = expm(A{mode1}*ts);
                  G_hat = inv(A{mode1})*(expm(A{mode1}*ts)-eye(2))*B;
                  pred = F_hat*KF_bank{mode1}.x +  G_hat*u;
                  x(:,w);
                  pred_cov = F_hat*KF_bank{mode1}.P*F_hat';
                  y_pred = H*pred;
                  y_pred_cov = H*pred_cov*H' + R;
                  temp = exp(-(y_pred-y_true)'*inv(y_pred_cov)*(y_pred-y_true));%y_pred_cov
                  if state_likelihood(i,w-1)*temp >= state_likelihood(j,w)
                      state_likelihood(j,w) = state_likelihood(i,w-1)*temp;
                      ml_transition(j,w) = i;
                      K_kf = pred_cov*H'*inv(y_pred_cov);
                      x_hat_ml = pred + K_kf*(y_true-H*pred);
                      P_ml = (eye(2)-K_kf*H)*pred_cov;
                  end
              end    
          end
       end
       KF_bank{i}.x=x_hat_ml;
       KF_bank{i}.P = P_ml;
    end
    state_likelihood(:,w)=state_likelihood(:,w)/sum(state_likelihood(:,w));
    KF_hist{w} = KF_bank;
end

% compute maximum likelihood paths from trellis
ml_path = zeros(1,W+1);
[max_val,max_ind] = max(state_likelihood(:,W));
ml_path(W+1) = max_ind(1);
for v = 1:W
% $$$     [max_val,max_ind] = max(state_likelihood(:,W+1-v).*M_tr(:,ml_path(W+2-v))');
% $$$     ml_path(W+1-v) = max_ind(1);
    ml_path(W+1-v) = ml_transition(ml_path(W+2-v),W+1-v);
end

for v = 2:W+1
    if use_constraints
        if ceil(ml_path(v)/2)==ceil(ml_path(v-1)/2)%ml_path(v)==ml_path(v-1)
            t_hats(v-1) = 0;
        end
    else
        if ml_path(v)==ml_path(v-1)
            t_hats(v-1) = 0;
        end    
    end
end

toc

KF_opt = {};
KF_opt.x = zeros(n,W-1);
KF_opt.P = zeros(n,n,W-1);
for v = 2:W;
   curr_KFs_k = KF_hist{v};
   KF_opt.x(:,v-1) = curr_KFs_k{ml_path(v)}.x;
   KF_opt.P(:,:,v-1) = curr_KFs_k{ml_path(v)}.P;
end

figure
plot(1:W-1,KF_opt.x(1,:))
hold on
plot(1:W-1,x(1,(2:W)+delay))
xlabel('time [s]')
ylabel('$i_L$ error','interpreter','Latex')

figure
plot(1:W-1,KF_opt.x(2,:))
hold on
plot(1:W-1,x(2,delay+(2:W)))
xlabel('time [s]')
ylabel('$v_0$ error','interpreter','Latex')