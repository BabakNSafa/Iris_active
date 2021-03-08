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
%%
time = raw(:,1);
r_p = raw(:,2);
% r_s = raw(:,3);
e_r = r_p/r_p(1)-1;

%resample time to not overfit the region where there is not light
fprintf('Find the middle index-1 of the time vector...\n')
mid_index = ceil(length(time)/2)-1;
fprintf('Down sample the initial idle time of the pupil (1:10), interp1, and moving average (3)...\n')
time_resample = [time(1:10:mid_index);time((mid_index):end)];  

%linear interpolation and moving average
e_r = interp1(time,e_r,time_resample,'linear','extrap');
e_r = movmean(e_r,3);

figure
hold on
plot(time_resample,e_r,'o')
plot(time_resample,'-*r')
%%
n = length(time_resample);
fun_par = @(x) [time_resample, x(1:n,1), x(n+1)*ones(size(time_resample))];
fun = @(x) sqrt(sum((FEBio_run_Iris_Active(fun_par(x))-e_r).^2)/length(e_r));

%test function
T_s_test = ones(size(time_resample));
T_d_test = ones(size(time_resample));

T_test = [time_resample,T_s_test,T_d_test];
test = fun(T_test);
if isempty(test)
    error('FEBio test failed!')
end

options = optimoptions('fmincon','Display','iter',...
            'MaxFunctionEvaluations',10000);

tic;
fprintf('Optimization started ...\n')
x_fit = fmincon(fun,ones(n+1,1),[],[],[],[],...
    ones(n+1,1)-ones(n+1,1),10*ones(n+1,1),[],options);
t_elapsed = toc;
fprintf('Optimization ended in %.3f hours\n',t_elapsed/3600)
%% Evaluate the fit response
e_r_fit = FEBio_run_Iris_Active(fun_par(x_fit));


T_fit = fun_par(x_fit);
T_s_fit = T_fit(:,2);
T_d_fit = T_fit(:,3);

M = [time_resample,e_r,e_r_fit,T_s_fit,T_d_fit];
writematrix(M,'Iris_active_results/results.csv');
M_title = {'Time','e_r','e_r_fit','T_s_Fit','T_d_Fit'};
writecell(M_title,'Iris_active_results/readme_DataColumnTitles.txt','Delimiter','tab');

%% Plot the resultes
figure('Units','centimeters','WindowStyle','normal','Position',[1,1,25,8]);

subplot(1,2,1)
hold on
plot(time_resample,e_r,'--b','linewidth',2)
plot(time_resample,e_r_fit,'--r','linewidth',2)

set(gca,'xlim',[0,Inf])
ylabel('Pupil diameter')
xlabel('Normalized Time')
set(gca,'FontName','Times','Fontsize',14)

subplot(1,2,2)
hold on
plot(time_resample,T_s_fit,'-','linewidth',2)
set(gca,'xlim',[0,Inf],'ylim',[0,Inf])

ylabel('$T_s$','Interpreter','latex')
xlabel('Time')

yyaxis right
plot(time_resample,T_d_fit,'-','linewidth',2)
ylabel('$T_d$','Interpreter','latex')
xlabel('Normalized Time')
set(gca,'xlim',[0,Inf],'ylim',[0,2])

legend({'Sphincter','Dillator'},'Interpreter','latex')
set(gca,'FontName','Times','Fontsize',14)

%% Save Outputs
movefile Iris_active.feb Iris_active_results/Iris_active.feb
movefile Iris_active.log Iris_active_results/Iris_active.log
movefile Iris_active.xplt Iris_active_results/Iris_active.xplt

%%
fprintf('%s \n',datetime)
fprintf('Analyses concluded! \n')

if isfile('Iris_active_results/MatlabDiary.txt')
    delete Iris_active_results/MatlabDiary.txt
end
diary Iris_active_results/MatlabDiary.txt
diary off
end