classdef IMM_KF_Solver
   properties
      n_modes = 0;
      T = 1;
      t_divs = 0.1;
      n_divs = 10;
      use_constraints = 0; 
       
      A_switch = {};
      A_no_switch = {};
      B_switch = {};
      B_no_switch = {};
      t_list = [];
      
      H = [];
      R = 1;
      
      x = {};
      P = {};
      mode_probs = [];
      t_hats = [];
      
      x_pred = {};
      P_pred = {};
      t_hats_pred = [];
      pred_probs = [];
      
      x_up = {};
      P_up = {};
      updated_probs = [];
   end
   
   methods
       function r = set_vals(this,this_x,this_P,this_probs,this_H,this_R)
           this.x = this_x;
           this.P = this_P;
           this.mode_probs = this_probs;
           this.H = this_H;
           this.R = this_R;
       end
       
       function r = compute_mats(n_d,A,B) % call this first to populate library
           for mode = 1:n_modes
               % discretization for no switch
               A_no_switch{mode} = expm(t_bar*A{mode});
               if(sum(eig(A{mode})<1e-3)>0)
                   B_switch{mode,div} = (t_bar*eye(size(A_no_switch{mode},1)) ...
                                             + t_bar^2*A{mode}/2)*B_no_switch{mode};
               else
                   B_switch{mode,div} = inv(A_{mode})*(A__no_switch{mode,div}-eye(2))*B{mode}
               end
               for div = 1:n_divs
                   t_bar = div*t_divs - t_divs/2;
                   t_list[div] = t_bar;
                   A_switch{mode,div} = expm(t_bar*A{mode});
                   if(sum(eig(A{mode})<1e-3)>0)
                       B_switch{mode,div} = (t_bar*eye(size(A{mode},1)) ...
                                             + t_bar^2*A{mode}/2)*B_no_switch{mode};
                   else
                       B_switch{mode,div} = inv(A_{mode})*(A_switch{mode,div}-eye(2))*B{mode}
                   end
               end
           end
           r = 0;
       end
       
       function r = predict_x_P(init_mode,fin_mode,u,y)
           if(init_mode == fin_mode)
              x_pred{init_mode} = A_no_switch{init_mode}*x{init_mode} + B_no_switch{init_mode}*u;
              P_pred{init_mode} = A_no_switch{init_mode}*P{init_mode}*A_no_switch{init_mode}';
           else
              x_divs = {};
              P_divs = {};
              cost_divs = [];
              for div = 1:n_divs 
                  x_divs{div} = A_switch{init_mode,div}*x{init_mode} + B_switch{div}*u; 
                  P_divs{div} = A_switch{init_mode,div}*P{init_mode}*A_switch{init_mode,div}'
                  x_divs{div} = A_switch{fin_mode,n_divs+1-div}*x_divs{div} + B_switch{fin_mode,n_divs+1-div}*u;
                  P_divs{div} = A_switch{fin_mode,n_divs+1-div}*P_divs{div}*A_switch{fin_mode,n_divs+1-div}';
                  cost_divs[div] = cost(x_divs{div},P_divs{div},y);
                  
              end
              [max_val,max_ind] = max(cost_divs);
              x_pred{init_mode} = x_divs{max_ind};
              P_pred{init_mode} = P_divs{max_ind};
              t_hats[init_mode] = t_list(max_ind);
           end
           r = 0;
       end
       
       function r = update_x_P(init_mode,fin_mode,y)
           K_k = P_pred{init_mode}*H'*inv(H*P_pred{init_mode}*H'+R);
           x_pred{init_mode} = x_pred{init_mod} + K_k(y-H*x_pred{init_mod}); 
           P_pred{init_mode} = (eye(size(A_switch{mode},1))-K_k*H)*P_pred{init_mode};
           pred_probs(init_mode) = cost(x_pred{init_mode},P_pred{init_mode},y)*mode_probs(init_mode);
           r = 0;
       end
       
       function cost = compute_likelihood(this_x,this_P,y)
           cost = (y-H*this_x)'*inv(R)*(y-H*this_x);
       end
       
       function r = compute_IMM_KF(u,y)
           temp_probs = zeros(n_modes);
           for i = 1:n_modes
               for j = 1:n_modes
                   predict_x_P(j,i,u,y);
                   update_x_P(j,i,y);
               end
               [max_val,max_ind] = max(pred_probs);
               x_up{i} = x_pred{max_ind};
               P_up{i} = P_pred{max_ind};
               updated_probs(i) = max_val;
           end
           x = x_up;
           P = P_up;
           mode_probs = updated_probs;
           r = 0;
       end
       
   end
end