function Iris_active()
%% Initialize
clc
close all
if ~isfolder('Iris_active_results')
    mkdir Iris_active_results
end
diary on
fprintf('Active Iris Curve Fitting\nWritten by Babak N. Safa (Ethier Lab)\n')
fprintf('\nMatlab version \t %s \n\n',version)
fprintf('Starting iris_active.m ...\n\n') 
fprintf('%s\n',datetime)

N = 25;
DataFit_Output.N = N;

fprintf('Iris_active.m run, with grid size of N = %d \n',N)

DataFit_Output.Start_date = datetime('now');
%% Read experimental data
raw = csvread('Iris_seg_summary.csv');

time = raw(:,1);
time(1) = 0; %make sure the initial point starts at time zero, this offsetgting does impact the fit, because this is duing the no-deformation phase
r_p = raw(:,2);
% r_s = raw(:,3);
e_r = 1/2*((r_p/mean(r_p(time<15))).^2-1);

%resample time to not overfit the region where there is not light
fprintf('Find the (middle index - 3) of the time vector...\n')
mid_index = ceil(length(time)/2)-3;
fprintf('Downsample the initial idle time of the pupil (1:10), interp1, and moving average (3)...\n')
time_resampled = [linspace(time(1),time(mid_index),5)';
                  linspace(time(mid_index),time(end),15)'];
time_resampled = unique(time_resampled);

%linear interpolation and moving average with three 3 neighboring points
e_r_resampled = interp1(time,e_r,time_resampled,'linear','extrap');
 
figure
hold on
plot(time,e_r,'o')
plot(time_resampled,e_r_resampled,'-*r')
%% Define cost function
%parameters v, tau (sec), T_sphincter
lb    = [1e-6,          1e-6,       1e-6];
ub    = [.49,           100,       100];

if sum(ub<lb)>0
    error('The bounds do not match')
end

fun_par = @(x) [1, x , 0]; %set the matrix modulus to one and dialtor traction to zero

DataFit_Output.par_name       = {'E', '\nu','\tau (sec)','T_sphincter', 'T_diallator'};
DataFit_Output.lb = fun_par(lb);
DataFit_Output.ub = fun_par(ub);

cleanup = true;
fun = @(x) sqrt(sum((FEBio_run_Iris_Active(fun_par(x),fun_par(lb),fun_par(ub),time_resampled,cleanup)-e_r_resampled).^2)/length(e_r_resampled));
%% Test function
v_test = 00;
tau_test = 3; 
T_s_test = 16; %in units of E

x_test = ([v_test, tau_test, T_s_test]-lb)./(ub-lb);
cleanup = false;
test1 = FEBio_run_Iris_Active(fun_par(x_test),fun_par(lb),fun_par(ub),time_resampled,cleanup);
test2 = fun(x_test);
if isempty(test1)||isempty(test2)
    error('FEBio test failed!')
end
%% Curvefitting
s=rng('shuffle');%for reproducing the random number
DataFit_Output.rand_state     = s;
M_par_normal_master = lhsdesign(N,length(lb)); %the list of random guesses the values are picked from
M_par_normal        = M_par_normal_master(1:N,:);

options = optimoptions(@fmincon,'Algorithm','interior-point',...
                      'Display','Iter','StepTolerance',1e-20);

tmp_swap=[];
fprintf('Optimization started ...\n')

parfor i=1:N
    sprintf('Starting cycle number %d/%d',[i,N])
    tic
    x0_temp = M_par_normal(i,:);

    [x_fit,fval,exitflag,~,~,grad,hessian] = fmincon(fun,x0_temp,[],[],[],[],...
            zeros(1,length(lb)),ones(1,length(ub)),[],options);
    t_elapsed = toc;
    fprintf('Optimization for cycle %d ended in %.3f hours\n',[i, t_elapsed/3600])
    %% Evaluate the fit response
    cleanup     = false;
    [e_r_fit, detail] = FEBio_run_Iris_Active(fun_par(x_fit),fun_par(lb),fun_par(ub),...
                        time_resampled,cleanup);
    RMSE_e_r    = sqrt(sum((e_r_fit-e_r_resampled).^2)/length(e_r_resampled));

    elapsed_time = toc
    
    fclose('all'); 
    if ~isfolder('temp')
        mkdir 'temp'    
    end
    movefile('incremental_creep.log',sprintf('temp/temp_%d_iris_active.log',i));
    movefile('incremental_creep.feb',sprintf('temp/temp_%d_iris_active.feb',i));
    movefile('incremental_creep.xplt',sprintf('temp/temp_%d_iris_active.xplt',i));
    
    tmp_swap(i).elapsed_time   = elapsed_time;   
    tmp_swap(i).x_par0         = fun_par(lb + x0_temp.*(ub-lb));
    tmp_swap(i).x_fit          = fun_par(lb + x_fit.*(ub-lb));
    tmp_swap(i).fval           = fval;
    tmp_swap(i).exitflag       = exitflag;
    tmp_swap(i).grad           = grad;
    tmp_swap(i).hessian        = hessian;
    tmp_swap(i).time           = time_resampled; 
    
    tmp_swap(i).e_r_fit        = e_r_fit;
    tmp_swap(i).e_r_0          = e_r_resampled;
    tmp_swap(i).RMSE_e_r       = RMSE_e_r;
    
    tmp_swap(i).detail         = detail;
end  
% %% Plot the resultes
% figure('Units','centimeters','WindowStyle','normal','Position',[1,1,25,8]);
% 
% subplot(1,2,1)
% hold on
% % plot(time,e_r,'ob')
% plot(time_resampled,e_r_resampled,'-ob','linewidth',2)
% plot(time_resampled,e_r_fit,'-*r','linewidth',2)
% 
% set(gca,'xlim',[0,Inf])
% ylabel('Pupil diameter')
% xlabel('Normalized Time')
% set(gca,'FontName','Times','Fontsize',14)
% 
% subplot(1,2,2)
% hold on
% plot(time_resampled,T_s_fit,'-','linewidth',2)
% set(gca,'xlim',[0,Inf],'ylim',[0,Inf])
% 
% ylabel('$T_s$','Interpreter','latex')
% xlabel('Time')
% 
% yyaxis right
% plot(time_resampled,T_d_fit,'-','linewidth',2)
% ylabel('$T_d$','Interpreter','latex')
% xlabel('Normalized Time')
% set(gca,'xlim',[0,Inf],'ylim',[0,2])
% 
% legend({'Sphincter','Dillator'},'Interpreter','latex')
% set(gca,'FontName','Times','Fontsize',14)
%% Save Outputs
DataFit_Output.active_iris = tmp_swap;

fprintf('%s \n',datetime)
fprintf('Analyses concluded! \n')

if isfile('temp/MatlabDiary.txt')
    delete temp/MatlabDiary.txt
end
formatOut = 'mm_dd_yyyy_HH_MM';
date_marker = datestr(now,formatOut);
save([sprintf('temp/DataFit_Output_N%d_',N),date_marker,'.mat'],'DataFit_Output')
diary temp/MatlabDiary.txt
diary off
end