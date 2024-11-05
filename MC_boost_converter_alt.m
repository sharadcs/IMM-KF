close all;
clear all;

seed = 195;
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

n_trials = 10;%0;%0;

sim_boost_converter

% sampline rate
ndivs = [1,2,3,4,5];%[1,2,5,10];
downsample = 50;

% discrete time measurements
H = [0 1];%eye(2);
R = 5;%50000;%*eye(2);

collected_MSE = zeros(2,length(downsample),n_trials);

for trial = 1:n_trials % run monte carlo trials for measurements
    
RMSE = zeros(2,length(ndivs));
ME = zeros(2,length(ndivs));
t_d = [];
% generate measurements
x_down = x(:,1:downsample:end);
t_down = t(1:downsample:end);
sigma_down = sigma(1:downsample:end);
y = x_down(2,:) + sqrt(R)*randn(size(x_down(2,:)));
    
for d = 1:length(ndivs)
    
    Ts = downsample*ts;
    t_div = Ts/ndivs(d);
    t_list = [];
    t_d(d) = t_div;%ndivs(d);
    
    % discretization no switch
    F1d = expm(A1*Ts);
    F2d = expm(A2*Ts);
    G1d = inv(A1)*(F1d-eye(2))*B;
    G2d = inv(A2)*(F2d-eye(2))*B;
    
    % discretization for switch
    % matrices for performing MLE on switch time
    A_switch = {};
    B_switch = {};
    F_switch = {};
    G_switch = {};
    for mode = 1:n_modes
        for div=1:ndivs(d)
            t_bar = div*t_div - t_div/2;
            t_list(div) = t_bar;
            A_switch{mode,div} = expm(t_bar*A{mode});
            B_switch{mode,div} = inv(A{mode})*(A_switch{mode,div}-eye(2))*B;
            F_switch{mode,div} = c2d_switch(t_bar,Ts,A{mode},A{3-mode});
            G_switch{mode,div} = input_switch(t_bar,Ts,A{mode},A{3-mode},B);
        end
    end
    
    % Kalman Filtering
    x_hat = zeros(size(x_down));
    x_hat(:,1) = 10*randn(2,1);
    P = zeros(2,2,size(x_down,2));
    P(:,:,1) = 100*eye(2);
    
    
    for k = 1:size(y,2)-1
        % pick state transition matrices
        if(sigma_down(k) == 1 && sigma_down(k+1) == 1)
            curr_F = F1d;
            curr_G = G1d;
        elseif(sigma_down(k) == 1 && sigma_down(k+1) == 2)
            curr_F = F1d;%F12_switch;%
            curr_G = G1d;%G12_switch;%
        elseif(sigma_down(k) == 2 && sigma_down(k+1) == 2)
            curr_F = F2d;%F21_switch;%
            curr_G = G2d;%G21_switch;%
        else
            curr_F = F2d;
            curr_G = G2d;
        end
        
        % perform filtering
        if(sigma_down(k+1)==sigma_down(k))
            x_pred = curr_F*x_hat(:,k) + curr_G*u;
            P_pred = curr_F*P(:,:,k)*curr_F';
        else
            min_val = inf;
            min_ind = 1;
            mode1 = sigma_down(k);
            mode2 = sigma_down(k+1);
            for div = 1:ndivs(d)
                x_pred = A_switch{mode1,div}*x_hat(:,k) + B_switch{mode1,div}*u;
                x_pred = A_switch{mode2,ndivs(d)+1-div}*x_pred + B_switch{mode2,ndivs(d)+1-div}*u;
                %x_pred = F_switch{mode1,div}*x_hat(:,k) + G_switch{mode1,div}*u;
                this_val = (y(:,k+1)-H*x_pred)'*inv(R)*(y(:,k+1)-H*x_pred);
                if this_val<min_val
                   min_val=this_val;
                   min_ind = div;
                   min_x = x_pred;
                end
                
                min_val;
                t_bar = k+1*Ts - switch_times((switch_times>=k*Ts) & (switch_times<(k+1)*Ts));
                t_bar_hat = min_ind*t_div - t_div/2;
                
                t_bar;
                t_bar_hat;
            end
            F_max1 = A_switch{mode1,min_ind};
            F_max2 = A_switch{mode2,ndivs(d)+1-min_ind};
            G_max1 = B_switch{mode1,min_ind};
            G_max2 = B_switch{mode2,ndivs(d)+1-min_ind};
            F_max = F_switch{mode1,min_ind};
            %P_pred = F_max*P_pred*F_max;
            P_pred = F_max2*F_max1*P(:,:,k)*F_max1'*F_max2';
            P_pred = P_pred + 10*eye(2);
            x_pred = min_x;
        end
        K_k = P_pred*H'*inv(H*P_pred*H' + R);
        x_hat(:,k+1) = x_pred + K_k*(y(:,k+1)-H*x_pred);
        P(:,:,k+1) = (eye(2)-K_k*H)*P_pred;
    end
    
    % trial plot
% $$$     figure
% $$$     plot(t_down,x_down(1,:));
% $$$     hold on
% $$$     plot(t_down,x_down(2,:));
% $$$     plot(t_down,x_hat(1,:));
% $$$     plot(t_down,x_hat(2,:));
    
    % compute error statistics
    error = x_down-x_hat;%x_down(:,1:end-1)-x_hat(:,2:end);
    RMSE(:,d) = sqrt(mean(error.^2,2));
    ME(:,d) = mean(error,2);
    
    collected_MSE(:,d,trial) = mean(error.^2,2);
    
% $$$     figure
% $$$     plot(t_down,x_down(1,:));
% $$$     hold on
% $$$     plot(t_down,x_hat(1,:));
% $$$     
% $$$     figure
% $$$     plot(t_down,x_down(2,:));
% $$$     hold on
% $$$     plot(t_down,x_hat(2,:));
end



end

% $$$ % plot data
% $$$ figure
% $$$ plot(t_d,ME(1,:),'b')
% $$$ hold on
% $$$ plot(t_d,ME(1,:)+RMSE(1,:),'b--')
% $$$ plot(t_d,ME(1,:)-RMSE(1,:),'b--')
% $$$ xlabel('Sampling Time [s]')
% $$$ ylabel('error statistics [Amps]')
% $$$ 
% $$$ figure
% $$$ plot(t_d,ME(2,:),'r')
% $$$ hold on
% $$$ plot(t_d,ME(2,:)+RMSE(2,:),'r--')
% $$$ plot(t_d,ME(2,:)-RMSE(2,:),'r--')
% $$$ xlabel('Sampling Time [s]')
% $$$ ylabel('error statistics [Volts]')

RMSE1 = sqrt(mean(collected_MSE(1,:,:),3));
RMSE2 = sqrt(mean(collected_MSE(2,:,:),3));

figure
plot(t_d,RMSE1)
%plot(t_d,ME(1,:),'b')
%hold on
%plot(t_d,ME(1,:)+RMSE(1,:),'b--')
%plot(t_d,ME(1,:)-RMSE(1,:),'b--')
xlabel('Sampling Time [s]')
ylabel('Root Mean-Squared Error [Amps]')

figure
plot(t_d,RMSE2)
%plot(t_d,ME(2,:),'r')
%hold on
%plot(t_d,ME(2,:)+RMSE(2,:),'r--')
%plot(t_d,ME(2,:)-RMSE(2,:),'r--')
xlabel('Sampling Time [s]')
ylabel('error statistics [Volts]')
