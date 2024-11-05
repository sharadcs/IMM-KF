function output = viterbi_kalman(x_0,P_0,N,meas,W,system,n,use_constraints)

y=meas.y;
H=meas.H;
R=meas.R;
ts=meas.ts;

A=system.A;
B=system.B;
u=system.u;

ndivs = 10;  % # of divisions of sample to check paths

% state definition:
% 1 - system 1, no switch last sample
% 2 - system 1, switched last sample
% 3 - system 2, no switch last sample
% 4 - system 2, switched last sample
if use_constraints
    n_states = n*n;%2;%4; % number of discrete states
else
    n_states = n;
end
    
% matrix of possible state transitions:
if use_constraints
    M_tr = [ 1 0 0 1; 1 0 0 0; 0 1 1 0; 0 0 1 0];%[1 1; 1 1];%[ 1 0 0 1; 1 0 0 0; 0 1 1 0; 0 0 1 0];
else
    M_tr = ones(n_states);%[1 1; 1 1];
end

% trellis defined l-r t-b of transitions
path_probs = zeros(size(M_tr,1)*size(M_tr,2),W); % transition probabilities
ml_transition = zeros(n_states,W+1); % best previous state for each state at each level
state_likelihood = zeros(n_states,W+1); % overall discrete state likelihoods
state_likelihood(:,1)=ones(n_states,1)/n_states;
t_hats = zeros(W+1,1);
KF_bank = {};
for i = 1:n_states
    KF_bank{i} = {};
    KF_bank{i}.x=x_0{i};%ceil(i/2)};
    KF_bank{i}.P=P_0{i};%{ceil(i/2)};
end

% loop through levels of trellis and update probabilities
tic
for w = 2:W
    y_true = y(:,w);
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
                      G_hat = input_switch(t_bar,ts,A{mode1},A{mode2},B);
                      pred = F_hat*KF_bank{mode1}.x + G_hat*u;
                      pred_cov = F_hat*KF_bank{mode1}.P*F_hat';
                      y_pred = H*pred;
                      y_pred_cov = H*pred_cov*H' + R;
                      temp = exp(-(y_pred-y_true)'*inv(R)*(y_pred-y_true));%y_pred_cov
                      pred_ml = pred_ml+temp/ndivs;
                      if temp>=curr_t_l
                          curr_t_l = temp;
                          t_bar_ml = t_bar;
                          K_kf = pred_cov*H'*inv(y_pred_cov);
                          maybe_x = pred + K_kf*(y_true-H*pred);
                          maybe_P = (eye(N)-K_kf*H)*pred_cov;
                      end
                  end
                  t_bar_ml;
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
                  G_hat = inv(A{mode1})*(expm(A{mode1}*ts)-eye(N))*B{mode1};
                  pred = F_hat*KF_bank{mode1}.x +  G_hat*u;
                  pred_cov = F_hat*KF_bank{mode1}.P*F_hat';
                  y_pred = H*pred;
                  y_pred_cov = H*pred_cov*H' + R;
                  temp = exp(-(y_pred-y_true)'*inv(R)*(y_pred-y_true));%y_pred_cov
                  if state_likelihood(i,w-1)*temp >= state_likelihood(j,w)
                      state_likelihood(j,w) = state_likelihood(i,w-1)*temp;
                      ml_transition(j,w) = i;
                      K_kf = pred_cov*H'*inv(y_pred_cov);
                      x_hat_ml = pred + K_kf*(y_true-H*pred);
                      P_ml = (eye(N)-K_kf*H)*pred_cov;
                  end
              end    
          end
       end
       KF_bank{i}.x=x_hat_ml;
       KF_bank{i}.P = P_ml;
    end
    state_likelihood(:,w)=state_likelihood(:,w)/sum(state_likelihood(:,w));
end

%state_likelihood

% compute maximum likelihood paths from trellis
ml_path = zeros(1,W);
[max_val,max_ind] = max(state_likelihood(:,W));
ml_path(W) = max_ind(1);
% $$$ for v = 1:W
% $$$ % $$$     [max_val,max_ind] = max(state_likelihood(:,W+1-v).*M_tr(:,ml_path(W+2-v))');
% $$$ % $$$     ml_path(W+1-v) = max_ind(1);
% $$$     v;
% $$$     ml_transition;
% $$$     ml_path;
% $$$     ml_path(W+1-v) = ml_transition(ml_path(W+2-v),W+1-v);
% $$$ end
% $$$ 
% $$$ for v = 2:W+1
% $$$     if use_constraints
% $$$         if ceil(ml_path(v)/2)==ceil(ml_path(v-1)/2)%ml_path(v)==ml_path(v-1)
% $$$             t_hats(v-1) = 0;
% $$$         end
% $$$     elseb
% $$$         if ml_path(v)==ml_path(v-1)
% $$$             t_hats(v-1) = 0;
% $$$         end    
% $$$     end
% $$$ end

map_state = ml_path(W);
for i=1:n_states
    x{i} = KF_bank{i}.x;
    P{i} = KF_bank{i}.P;
end

output.x=x;
output.P=P;
output.map_state=map_state;