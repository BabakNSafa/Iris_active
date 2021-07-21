function Iris_active()
%% Initialize
clc
close all
parpool('local',24)

if ~isfolder('temp')
    mkdir 'temp'    
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
    %% Seg
    raw_seg = readtable('Iris_seg_summary.csv');

    time = raw_seg.Time_sec_;
    e_r_p = raw_seg.e_pupil;

    %linear interpolation
    index = (time>17)&(time<20);
    e_r_p_max_0 = mean(e_r_p(index));

    time_resampled = [0;.1; 1];% the .1 step is added to make sure the reference state of elements is logged in the log file
    T_s_lc            = [0;0;1];
    T_d_lc            = [0;0;0];
    load = [time_resampled,T_s_lc,T_d_lc];
    step_size = .1; 
    figure
    hold on
    plot(time,e_r_p,'o')
    plot(time(index),e_r_p(index),'-*r')
%% Define cost function
%parameters E (MPa), v, tau (sec), beta, T_sphincter(MPa), a_sphincter(mm) T_dialator(MPa)
lb    = [0, 0, 1e-10,  2,  0, .1,  0];
ub    = [1, .49, 100, 100, 1, 2, 1];

if sum(ub<lb)>0
    error('The bounds do not match')
end

%set the sphincter width (see
%/home/asixbabak/Dropbox (GaTech)/project-IrisBiomechanics/Iris sphincter
%dimension estimation from literature)
% test cases are a_s = [.4, .7, 1, 1.3]

fun_par_normal = @(x) [x(1),   x(2),   0,  0,  x(3),  (.7-lb(6))/(ub(6)-lb(6)), 0]; %set the time constat, beta, and dialtor traction to the minimum of lb


%in the variable names I have kept the usual FEBio units when using mm as dimension, which are based
%on mm, MPa and sec
DataFit_Output.par_name       = {'E (MPa)','\nu','\tau (sec)','\beta','T_sphincter(MPa)', 'a_sphincter (mm)', 'T_diallator(MPa)'};
DataFit_Output.lb = lb;
DataFit_Output.ub = ub;

%% Define the objective function
    %% parameters for testing
E_test      = .01;% in MPa, which is equivalent of 10 kPa
v_test      = 0.4;
tau_test    = 100;%in sec
beta_test   = 2; 
T_s_test    = .01;% in MPa, which is equivalent of 10 kPa
a_s_test    = 1;% in mm
T_d_test    = .02;% in MPa, which is equivalent of 20 kPa

x_test = [E_test, v_test, tau_test, beta_test T_s_test, a_s_test, T_d_test];
temp = (x_test-lb)./(ub-lb);
x_test_normal = fun_par_normal(temp([1,2,5])); 
cleanup = true;
[e_r_test1, ~, detail] = FEBio_run_Iris_Active(x_test_normal,lb,ub,load,step_size,cleanup);
rho = detail.rho;

    %% Definitions go in here
cleanup = true;

%Fit only maximum constriction pupil strain
fun = @(x) sqrt((min(FEBio_run_Iris_Active(fun_par_normal(x),lb,ub,...
    load,step_size,cleanup))-e_r_p_max_0).^2);

    %% run tests
test1 = sqrt(sum((min(e_r_test1)-e_r_p_max_0).^2)/length(e_r_p_max_0));
test2 = fun(x_test_normal([1,2,5]));

if isempty(test1)||isempty(test2)||abs(test1-test2)>1e-12
    error('FEBio test failed!')
end
%% Curvefitting
s=rng('shuffle');%for reproducing the random number
DataFit_Output.rand_state     = s;
M_par_normal_master = lhsdesign(N,3); %the list of random guesses the values are picked from
M_par_normal        = M_par_normal_master(1:N,:);
lb_normal = zeros(1,size(M_par_normal,2));
ub_normal = ones(1,size(M_par_normal,2));

options = optimoptions(@fmincon,'Algorithm','interior-point',...
                      'Display','Iter','DiffMinChange',.01,...
                      'ObjectiveLimit',1e-5,'StepTolerance',1e-20);

tmp_swap=[];
fprintf('Optimization started ...\n')

parfor i=1:N
    sprintf('Starting cycle number %d/%d',[i,N])
    tic
    x0_temp = M_par_normal(i,:);

    [x_fit,fval,exitflag,~,~,grad,hessian] = fmincon(fun,x0_temp,[],[],[],[],...
            lb_normal,ub_normal,[],options);

    t_elapsed = toc;
    fprintf('Optimization for cycle %d ended in %.3f hours\n',[i, t_elapsed/3600])
    %% Evaluate the fit response
    cleanup     = false;
    [e_r_p_max_fit, ~, detail] = FEBio_run_Iris_Active(fun_par_normal(x_fit),lb,ub,...
                        load,step_size,cleanup);

    RMSE_e_r    = fval;

    elapsed_time = toc
    
    pause(2)
    if isfile('Iris_active.log')
        movefile('Iris_active.log',sprintf('./temp/temp_%d_iris_active.log',i),'f');
    end
    if isfile('Iris_active.feb')
        movefile('Iris_active.feb',sprintf('./temp/temp_%d_iris_active.feb',i),'f');
    end
    if isfile('Iris_active.xplt')
        movefile('Iris_active.xplt',sprintf('./temp/temp_%d_iris_active.xplt',i),'f');
    end
    
    tmp_swap(i).elapsed_time   = elapsed_time;   
    tmp_swap(i).x_par0         = lb + fun_par_normal(x0_temp).*(ub-lb);
    tmp_swap(i).x_fit          = lb + fun_par_normal(x_fit).*(ub-lb);
    tmp_swap(i).fval           = fval;
    tmp_swap(i).exitflag       = exitflag;
    tmp_swap(i).grad           = grad;
    tmp_swap(i).hessian        = hessian;
    tmp_swap(i).time           = time_resampled; 
    
    tmp_swap(i).rho            = rho;
    tmp_swap(i).e_r_p_max_fit  = e_r_p_max_fit;
    tmp_swap(i).e_r_p_max_0    = e_r_p_max_0;
    tmp_swap(i).RMSE_e_r       = RMSE_e_r;
    
    tmp_swap(i).detail         = detail;
    pause(1)
end  
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
diary off

move diary temp/MatlabDiary
end