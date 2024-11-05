function [Solvers,n_solvers] = generateDubinsMHE(window_len,M,n_switches)
    
% set time steps for optimization
Tvariable tm [];
N=4;

% enumerate possible switching signals
max_switches = n_switches;
n_solvers = window_len; % TODO: set to correct number of solvers needed
%num_switch_signals = factorial(n_modes)/factorial(n_modes-1-max_switches)
%switch_bank = {};
%for i = 1:max_switches+1
%    switch_bank{i} = nchoosek(1:n_modes,i);
%end

%Tvariable active_modes window_len;
Tvariable p [2,window_len];
Tvariable a [2,window_len];
Tvariable pm [2,1];
Tvariable am [2,1];

Tvariable beacon_pos [2,M];
Tvariable y [M,window_len];
Tvariable Rinv [];

% measurement cost function (max likelihood)
pp=reshape(p,[2,1,window_len]);
pp=pp(:,ones(M,1),:);
bb_pos=reshape(beacon_pos,[2,M,1]);
bb_pos=bb_pos(:,:,ones(window_len,1));
rel_pos=pp-bb_pos; %[3,nBeacons,nT]
distance=sqrt(tprod(rel_pos,[-1,1,2],rel_pos,[-1,1,2])); %[nBeacons,nT]
cost_meas = norm2(y-distance);

% State and Input transition matrices
Tvariable ts window_len;

% dynamics contstraints
Tvariable v [];
Tvariable w [];
%Tvariable dm [N,max_switches];
%Tvariable active_const [3*window_len];

%cost_proc = 100*(norm2(d));%+norm2(dm));

constraints = {};
n_solvers;
window_len;
for solv_num = 1:n_solvers
curr_constraints = {};

% constraints on dynamics
    if solv_num == window_len
        dt = ts(2:window_len)-ts(1:window_len-1);
        dt = reshape(dt,[1,window_len-1]);
        dp = p(:,2:end) - p(:,1:end-1);
        da = a(:,2:end) - a(:,1:end-1);
        
        vdt = tprod(repmat(reshape(v,[1,1]),[2,window_len-1]),[1,2],repmat(dt,[2,1]),[1,2]);
        curr_constraints{end+1} = dp == tprod(vdt,[1,2],a(:,2:window_len),[1,2]);
        ww = [w*ones(1,window_len-1); -w*ones(1,window_len-1)];
        wdt = tprod(ww,[1,2],repmat(dt,[2,1]),[1,2]);
        curr_constraints{end+1} = da == tprod(wdt,[1,2],a([2,1],2:window_len),[1,2]);
    else
        p_aug = [p(:,1:solv_num) pm p(:,solv_num+1:window_len)];
        a_aug = [a(:,1:solv_num) am a(:,solv_num+1:window_len)];
        t_aug = [ts(1:solv_num); tm; ts(solv_num+1:window_len)];
        dt = t_aug(2:window_len+1)-t_aug(1:window_len);
        dt = reshape(dt,[1,window_len]);
        dp = p_aug(:,2:end) - p_aug(:,1:end-1);
        da = a_aug(:,2:end) - a_aug(:,1:end-1);

        vdt = tprod(repmat(reshape(v,[1,1]),[2,window_len]),[1,2],repmat(dt,[2,1]),[1,2]);
        curr_constraints{end+1} = dp == tprod(vdt,[1,2],a_aug(:,2:window_len+1),[1,2]);
        ww = [w*ones(1,window_len); -w*ones(1,window_len)];
        wdt = tprod(ww,[1,2],repmat(dt,[2,1]),[1,2]);
        curr_constraints{end+1} = da == tprod(wdt,[1,2],a_aug([2,1],2:window_len+1),[1,2]);
    end

% constraints on switch times
if solv_num<window_len
    tmax_constr = tm <= ts(solv_num+1);
    tmin_constr = tm >= ts(solv_num);
    curr_constraints{end+1} = tmax_constr;
    curr_constraints{end+1} = tmin_constr;
else
    %curr_constraints{end+1} = tmax_constr;
end

% constraints for angle
if solv_num<window_len
    curr_constraints{end+1} = tprod(a_aug,[-1,1],a_aug,[-1,1])==Tones(window_len+1);
else
    curr_constraints{end+1} = tprod(a,[-1,1],a,[-1,1])==Tones(window_len);
end

constraints{end+1} = curr_constraints;
end



% constraints on state
% Tvariable x_max [N,window_len+max_switches];
% Tvariable x_min [N,window_len+max_switches];
% x_max_constr = x_aug <= x_max;
% x_min_constr = x_aug >= x_min;
%          

%constraints{end+1} = x_max_constr;
%constraints{end+1} = x_min_constr;
         
%Tvariable Hess_ [97,97];%[121,121]%(2*N*(window_len+max_switches) + 1 + 2*N*(window_len+max_switches) +2)*[1,1];
%Tvariable dHess_ size(Hess_,1);  % smaller

cost = cost_meas;% + norm2(x);%+cost_proc;

Solvers = {};
for solv_num = 1:n_solvers
this_constraints = constraints{solv_num};
name = sprintf('tmpDubinsFilter_%d',solv_num);
if solv_num<n_solvers
    classname=class2optimizeCS('classname',name,...
                          'objective',cost,...
                          'optimizationVariables',{...
                              p,a,pm,am,tm},...
                          'parameters',{...
                              y,Rinv,ts,v,w,beacon_pos},...
                          'constraints',this_constraints,...
                          'outputExpressions',{...
                              cost,p,a,pm,am,tm},...
                          'addEye2Hessian',false,...
                          'solverVerboseLevel',3);
else
    classname=class2optimizeCS('classname',name,...
                          'objective',cost,...
                          'optimizationVariables',{...
                              p,a},...
                          'parameters',{...
                              y,Rinv,ts,v,w,pm,am,tm,beacon_pos},...
                          'constraints',this_constraints,...
                          'outputExpressions',{...
                              cost,p,a},...
                          'addEye2Hessian',false,...
                          'solverVerboseLevel',3);
end
this_solver = feval(classname);
Solvers{solv_num} = this_solver;
end

end

