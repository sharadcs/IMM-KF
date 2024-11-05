close all;
clear all;

seed = 195;
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

n_trials = 1;

sim_boost_converter

% sampline rate
downsample = 5*[2,3,7,10,15];%[1,2,5,10];

% discrete time measurements
H = [0 1];%eye(2);
R = 5;%50000;%*eye(2);

collected_MSE = zeros(2,length(downsample),n_trials);

for trial = 1:n_trials % run monte carlo trials for measurements
    
RMSE = zeros(2,length(downsample));
ME = zeros(2,length(downsample));
t_d = [];
for d = 1:length(downsample)
    % generate measurements
    x_down = x(:,1:downsample(d):end);
    t_down = t(1:downsample(d):end);
    sigma_down = sigma(1:downsample(d):end);
    y = x_down(2,:) + sqrt(R)*randn(size(x_down(1,:)));
    
    Ts = downsample(d)*ts;
    t_d(d) = Ts;
    t_div = Ts;
    
    % discretization no switch
    F1d = expm(A1*Ts);
    F2d = expm(A2*Ts);
    G1d = inv(A1)*(F1d-eye(2))*B;
    G2d = inv(A2)*(F2d-eye(2))*B;
    
    % discretization for switch
    F12_switch = c2d_switch(Ts/2,Ts,A1,A2);
    F21_switch = c2d_switch(Ts/2,Ts,A2,A1);
    G12_switch = input_switch(Ts/2,Ts,A1,A2,B);
    G21_switch = input_switch(Ts/2,Ts,A2,A1,B);
    
    A_switch = {};
    B_switch = {};
    F_switch = {};
    G_switch = {};
    n_modes = 2;
    for mode = 1:n_modes
        div = 1;
        t_bar = div*t_div - t_div/2;
        t_list(div) = t_bar;
        A_switch{mode,div} = expm(t_bar*A{mode});
        B_switch{mode,div} = inv(A{mode})*(A_switch{mode,div}-eye(2))*B;
        F_switch{mode,div} = c2d_switch(t_bar,Ts,A{mode},A{3-mode});
        G_switch{mode,div} = input_switch(t_bar,Ts,A{mode},A{3-mode},B);
    end
    
    F12 = expm(A1*Ts/2);
    F21 = expm(A2*Ts/2);
    G12 = inv(A1)*(F12-eye(2))*B;
    G21 = inv(A2)*(F21-eye(2))*B;
    
    % Kalman Filtering
    x_hat = zeros(size(x_down));
    x_hat(:,1) = 10*randn(2,1);
    P = zeros(2,2,size(x_down,2));
    P(:,:,1) = 100*eye(2);
    
    for k = 1:size(y,2)-1
% $$$         % pick state transition matrices
% $$$         if((sigma_down(k) == 1) && (sigma_down(k+1) == 1))
% $$$             %curr_F = F1d;
% $$$             %curr_G = G1d;
% $$$             x_pred = F1d*x_hat(:,k) + G1d*u;
% $$$             P_pred = F1d*P(:,:,k)*F1d';
% $$$         elseif((sigma_down(k) == 1) && (sigma_down(k+1) == 2))
% $$$             %curr_F = F1d;%F12_switch;%
% $$$             %curr_G = G1d;%G12_switch;%
% $$$             x_pred = F12*x_hat(:,k) + G12*u;
% $$$             P_pred = F12*P(:,:,k)*F12';
% $$$             x_pred = F21*x_pred + G21*u;
% $$$             P_pred = F21*P_pred*F21';
% $$$         elseif((sigma_down(k) == 2) && (sigma_down(k+1) == 2))
% $$$             %curr_F = F2d;%F21_switch;%
% $$$             %curr_G = G2d;%G21_switch;%
% $$$             x_pred = F21*x_hat(:,k) + G21*u;
% $$$             P_pred = F21*P(:,:,k)*F21';
% $$$             x_pred = F12*x_pred + G12*u;
% $$$             P_pred = F12*P_pred*F12';
% $$$         else
% $$$             %curr_F = F2d;
% $$$             %curr_G = G2d;
% $$$             x_pred = F2d*x_hat(:,k) + G2d*u;
% $$$             P_pred = F2d*P(:,:,k)*F2d';
% $$$         end
        
        % perform filtering
% $$$         x_pred = curr_F*x_hat(:,k) + curr_G*u;
% $$$         P_pred = curr_F*P(:,:,k)*curr_F';
% $$$         if(sigma_down(k+1)~=sigma_down(k))
% $$$             P_pred = P_pred + 10*eye(2);
% $$$         end
        
        if(sigma_down(k)==1)
            curr_F = F1d;
            curr_G = G1d;
        else
            curr_F = F2d;
            curr_G = G2d;
        end
        
                % perform filtering
        if(sigma_down(k+1)==sigma_down(k))
            x_pred = curr_F*x_hat(:,k) + curr_G*u;
            P_pred = curr_F*P(:,:,k)*curr_F';
        else
            mode1 = sigma_down(k);
            mode2 = sigma_down(k+1);
            
            x_pred = A_switch{mode1,1}*x_hat(:,k) + B_switch{mode1,1}*u;
            x_pred = A_switch{mode2,1}*x_pred + B_switch{mode2,1}*u;
            %x_pred = F_switch{mode1,div}*x_hat(:,k) + G_switch{mode1,div}*u;                
           
            F_max1 = A_switch{mode1,1};
            F_max2 = A_switch{mode2,1};
            G_max1 = B_switch{mode1,1};
            G_max2 = B_switch{mode2,1};
            F_max = F_switch{mode1,1};
            %P_pred = F_max*P_pred*F_max;
            P_pred = F_max2*F_max1*P(:,:,k)*F_max1'*F_max2';
            P_pred = P_pred + 10*eye(2);
        end
        
        K_k = P_pred*H'*inv(H*P_pred*H' + R);
        x_hat(:,k+1) = x_pred + K_k*(y(:,k+1)-H*x_pred);
        P(:,:,k+1) = (eye(2)-K_k*H)*P_pred;
    end
    
    % compute error statistics
    error = x_down-x_hat;%x_down(:,1:end-1)-x_hat(:,2:end);%
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

% plot data
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
