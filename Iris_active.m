function Iris_active()
%% Initialize
clc
close all
diary on

fprintf('Active Iris Curve Fitting\nWritten by Babak N. Safa (Ethier Lab)\n')
fprintf('\nMatlab version \t %s \n\n',version)
fprintf('Starting iris_active.m ...\n\n') 
fprintf('%s\n',datetime)

N = 10;
DataFit_Output.N = N;

fprintf('Iris_active.m run, with grid size of N = %d \n',N)

DataFit_Output.Start_date = datetime('now');
%% Read experimental data
raw = csvread('Iris_seg_summary.csv');

time = raw(:,1);
r_p = raw(:,2);

[~, middle_index ] = min(abs(time-15));
if time(middle_index)<15
    middle_index = middle_index-1;
end
r_p_ref = mean(r_p(1:middle_index));

raw = raw(middle_index+1:end,:);

time = raw(:,1);
time = time - time(1);
r_p = raw(:,2);
e_r = 1/2*((r_p/r_p_ref).^2-1);

time_resampled = [linspace(0,2,30)';
                  linspace(2,15,60)'];
time_resampled = unique(time_resampled);

%linear interpolation
e_r_resampled = interp1(time,e_r,time_resampled,'linear','extrap');
 
figure
hold on
plot(time,e_r,'o')
plot(time_resampled,e_r_resampled,'-*r')
%% Define cost function
%parameters v, tau (sec), beta, T_sphincter
lb    = [1e-6, 1e-6,          1e-6,   1,       1e-6];
ub    = [100, .49,           100,    100,       40];

if sum(ub<lb)>0
    error('The bounds do not match')
end

fun_par = @(x) [x , 0]; %set the matrix dialtor traction to zero

DataFit_Output.par_name       = {'\nu','\tau1 (sec)','\tau2/\tau1','T_sphincter', 'T_diallator'};
DataFit_Output.lb = fun_par(lb);
DataFit_Output.ub = fun_par(ub);

cleanup = true;
fun = @(x) sqrt(sum((FEBio_run_Iris_Active(fun_par(x),fun_par(lb),fun_par(ub),time_resampled,cleanup)-e_r_resampled).^2)/length(e_r_resampled));
%% Test function
E_test = 1;
v_test = 0.4;
tau_test = 1;
beta_test = 5;
T_s_test = 10;

x_test = ([E_test, v_test, tau_test, beta_test, T_s_test]-lb)./(ub-lb);
cleanup = false;
test1 = FEBio_run_Iris_Active(fun_par(x_test),fun_par(lb),fun_par(ub),time_resampled,cleanup);
test2 = fun(x_test);
hold on
plot(time_resampled,test1,'o-')
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
    hold on
    plot(time_resampled,e_r_fit,'-*b')

    RMSE_e_r    = sqrt(sum((e_r_fit-e_r_resampled).^2)/length(e_r_resampled));

    elapsed_time = toc
    
    if ~isfolder('temp')
        mkdir 'temp'    
    end
    
    fclose 'all';
    
    movefile('Iris_active.log',sprintf('./temp/temp_%d_iris_active.log',i));
    movefile('Iris_active.feb',sprintf('./temp/temp_%d_iris_active.feb',i));
    movefile('Iris_active.xplt',sprintf('./temp/temp_%d_iris_active.xplt',i));
    
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