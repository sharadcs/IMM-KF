close all
clear all

A1 = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
A2 = [0 1 0 0; 0 1 0 -1; 0 0 0 1; 0 -1 0 1];
A3 = [0 1 0 0; 0 -1 0 1; 0 0 0 1; 0 1 0 -1];
B1 = [0 1 0 1]';
B2 = [0 -1 0 1]';
B3 = [0 1 0 -1]';
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

% loop through trials
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

switch_times = [0.3 0.4];
switch_pattern = [1 2 3];

% simulate trasjectory
x = zeros(N,K);
t = zeros(K);
sigma = zeros(K);
u=0.11;
curr_k = 1;
for i=1:length(switch_pattern)-1
    switch_pattern(i);
    switch_times(i);
    sys_mode = ss(A{switch_pattern(i)},B{switch_pattern(i)},eye(N),zeros(N,1));
    if i==length(switch_pattern)-1
       end_time = T; 
    else
        end_time = ts*floor(switch_times(i)/ts);
    end
    t_mode = curr_k*ts:ts:end_time;
    next_k = floor(end_time/ts);
    x_mode = lsim(sys_mode,u*ones(size(t_mode)),t_mode-ts*curr_k,x(:,curr_k),'zoh');
    x(:,curr_k+1:next_k) = x_mode(2:end,:)';
    curr_k=next_k;
end

figure
plot(x(1,:),x(3,:))

end