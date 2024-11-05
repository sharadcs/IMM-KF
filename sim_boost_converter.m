u = 100;
R1 = 2;
L1 = 500e-6;%1000e-6;
C0 = 470e-6;%1000e-6;
R0 = 50;

% continuous time state dynamics
A1 = [-R1/L1 0; 0 -1/(R0*C0)];
A2 = [-R1/L1 -1/L1; 1/C0 -1/(R0*C0)];

B = [1/L1; 0];
A = {A1,A2};

n_modes = 2;

T=0.5;%1;
ts=0.00001;
sigma = [1,2,1,2,1];
switch_periods = [0.0012:0.0012:0.2496 0.2505:0.0009:0.4998];% 0.5005:0.0007:0.7497 0.7503:0.0006:1];
switch_times = [switch_periods];% (switch_periods-0.001)];%[0.4444,2.1111,3.8888,5.6666];
K=floor(T/ts);

% discretized system
F1 = expm(A1*ts);
F2 = expm(A2*ts);
F = {F1,F2};
G1 = inv(A1)*(F1-eye(2))*B;
G2 = inv(A2)*(F2-eye(2))*B;
G = {G1,G2};

sys1 = ss(A1,B,eye(2),zeros(2,1));
sys2 = ss(A2,B,eye(2),zeros(2,1));
sys = {sys1,sys2};

%simulate state development
x=zeros(2,K);
t=zeros(1,K);
sigma=zeros(1,K);
x_0 = [0;0];
sigma_0 = 1;
x_k = x_0;
sigma_k = sigma_0;
t0 = 0;
tk = t0;
% $$$ for k=1:K
% $$$     tk = ts*k;
% $$$     
% $$$     simed = lsim(sys{sigma_k},[u u]',[0 ts],x_k,'zoh');
% $$$     x_k = simed(2,:)';%F{sigma_k}*x_k+G{sigma_k}*u;
% $$$     
% $$$     if(sum((switch_times<(tk+ts/2)).*(switch_times>(tk-ts/2))))
% $$$         sigma_k = 3-sigma_k;
% $$$     end
% $$$     
% $$$     x(:,k)=x_k;
% $$$     t(k) = tk;
% $$$     sigma(k) = sigma_k;
% $$$ end
curr_ind = 1;
for k = 1:length(switch_periods)
    t_beg = tk;
    t_end = switch_periods(k);
    
    floor((t_end-t_beg)/ts);
    
    curr_window = round((t_end-t_beg)/ts);
    t_window = 0:ts:ts*curr_window;;
    
    simed = lsim(sys{sigma_k},u*ones(size(t_window)),t_window,x_k,'zoh');
    
    x(:,curr_ind:curr_ind+curr_window) = simed';
    t(curr_ind:curr_ind+curr_window) = t_beg:ts:t_end;
    sigma(:,curr_ind:curr_ind+curr_window-1) = sigma_k*ones(1,curr_window);
    
    sigma_k = 3-sigma_k;
    
    curr_ind = curr_ind+curr_window;
    x_k = simed(end,:)';
    tk = t_end;
end

x = x(:,1:end-20);
t = t(1:end-20);
sigma = sigma(1:end-20);


%x = [x_0 x];
%t = [0 t];
%sigma = [sigma_0 sigma];

fig = figure
yyaxis left
plot(t(23000:27000)-t(23000),x(1,23000:27000),'--')
hold on
left_color = [0 0 0];
right_color = [0 0 0];
%set(fig,'defaultAxesColorOrder',[left_color; right_color]);
ylim([-10,130])
yyaxis right
plot(t(23000:27000)-t(23000),x(2,2MC3000:27000))
ylim([-10,130])
xlabel('time [s]')
xlim([0,0.04])
set(gca,'FontSize',12)
yyaxis left
ylabel('$i_L$ [amps]','interpreter','latex','FontSize',16)
yyaxis right
ylabel('$v_0$ [volts]','interpreter','latex','FontSize',16)