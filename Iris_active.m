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
%% Read experimental data
raw = csvread('Iris_seg_summary.csv');

time = raw(:,1);
time(1) = 0; %make sure the initial point starts at time zero, this offsetgting does impact the fit, because this is duing the no-deformation phase
r_p = raw(:,2);
% r_s = raw(:,3);
e_r = r_p/r_p(1)-1;

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
n = length(time_resampled);
cleanup = true;
fun_par = @(x) [time_resampled, x(1:n,1), x(n+1)*ones(size(time_resampled))];
fun = @(x) sqrt(sum((FEBio_run_Iris_Active(fun_par(x),cleanup)-e_r_resampled).^2)/length(e_r_resampled));
%% Test function
% T_s_test = ones(size(time_resampled));
% T_d_test = ones(size(time_resampled));
% 
% T_test = [time_resampled,T_s_test,T_d_test];
% cleanup = true;
% test1 = FEBio_run_Iris_Active(fun_par(T_test),cleanup);
% test2 = fun(T_test);
% if isempty(test1)||isempty(test2)
%     error('FEBio test failed!')
% end
%% Curvefitting
options = optimoptions('fmincon','Display','iter',...
            'UseParallel',false,...
            'Objectivelimit',1e-3);

tic;
fprintf('Optimization started ...\n')
x_fit = fmincon(fun,ones(n+1,1),[],[],[],[],...
    ones(n+1,1)-ones(n+1,1),10*ones(n+1,1),[],options);
t_elapsed = toc;
fprintf('Optimization ended in %.3f hours\n',t_elapsed/3600)
%% Evaluate the fit response
cleanup = false;
e_r_fit = FEBio_run_Iris_Active(fun_par(x_fit),cleanup);

T_fit = fun_par(x_fit);
T_s_fit = T_fit(:,2);
T_d_fit = T_fit(:,3);

M = [time_resampled,e_r_resampled,e_r_fit,T_s_fit,T_d_fit];
writematrix(M,'Iris_active_results/results.csv');
M_title = {'Time','e_r','e_r_fit','T_s_Fit','T_d_Fit'};
writecell(M_title,'Iris_active_results/readme_DataColumnTitles.txt','Delimiter','tab');
%% Plot the resultes
figure('Units','centimeters','WindowStyle','normal','Position',[1,1,25,8]);

subplot(1,2,1)
hold on
% plot(time,e_r,'ob')
plot(time_resampled,e_r_resampled,'-ob','linewidth',2)
plot(time_resampled,e_r_fit,'-*r','linewidth',2)

set(gca,'xlim',[0,Inf])
ylabel('Pupil diameter')
xlabel('Normalized Time')
set(gca,'FontName','Times','Fontsize',14)

subplot(1,2,2)
hold on
plot(time_resampled,T_s_fit,'-','linewidth',2)
set(gca,'xlim',[0,Inf],'ylim',[0,Inf])

ylabel('$T_s$','Interpreter','latex')
xlabel('Time')

yyaxis right
plot(time_resampled,T_d_fit,'-','linewidth',2)
ylabel('$T_d$','Interpreter','latex')
xlabel('Normalized Time')
set(gca,'xlim',[0,Inf],'ylim',[0,2])

legend({'Sphincter','Dillator'},'Interpreter','latex')
set(gca,'FontName','Times','Fontsize',14)
%% Save Outputs
movefile Iris_active.feb Iris_active_results/Iris_active.feb
movefile Iris_active.log Iris_active_results/Iris_active.log
movefile Iris_active.xplt Iris_active_results/Iris_active.xplt

fprintf('%s \n',datetime)
fprintf('Analyses concluded! \n')

if isfile('Iris_active_results/MatlabDiary.txt')
    delete Iris_active_results/MatlabDiary.txt
end
diary Iris_active_results/MatlabDiary.txt
diary off
end