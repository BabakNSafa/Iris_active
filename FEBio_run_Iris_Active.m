function [e_r_pupil,varargout] = FEBio_run_Iris_Active(x_par_normal,lb,ub,...
                                                       load,step_size,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function reads the jobs/Iris_active.feb and makes
%changes to it
%
%function [dia_pupil,varargout] = FEBio_run_Iris_Active(T,varargin)
% Inputs 
%     x_par_normal (1x5) = [E, v, tau T_Sphicter, T_Diallator]// this
%     variable is normalized between 0 and 1
%     E                 = Matrix modulus
%     v                 = Matrix Poisson's ratio
%     tau               = Matrix kinetics time-constant
%     beta              = time constant power
%     T_Sphicter        = Sphicter muscle stress as a function of time (in units of E)
%     T_Diallator       = Diallator muscle stress as a function of time (in units of E)
%
%     lb                = lower bound for parameters
%     ub                = upper bound for parameters
%
%     load              = [Time (sec), T_s_lc, T_d_lc]; normalized load
%     curve T_s_lc and  T_d_lc between 0 and 1
%     step_size         = the size of time step (sec);
%
%     varargin{1}       = cleanup
%     
%
% Outputs:
%     e_r_pupil(vector)     = raidual strain of pupil as a function of time.
%     detail                = structure that includes element position, stress,
%                     strain, and fluid flux (only xx yy zz elements,
%                     because the model is symmetric avoided all the tensor
%                     componenets, basically a bunch of zeros)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cleanup = true; %default value to cleanup the mess after febio
if ~isempty(varargin)
    cleanup = varargin{1};
end
nOutputs = nargout;

%denormalizing the input model parameters
x = lb + x_par_normal.*(ub-lb);
x(isnan(x)) = 0;%if any of the parameters is NaN make it zero/ this could happen for T_Diallator
E   = x(:,1);
v   = x(:,2);
tau = x(:,3);
beta = x(:,4);
T_Sphicter  = x(:,5);
a_Sphincter = x(:,6);
T_Diallator = x(:,7);

time_resample = load(1,1):step_size:load(end,1);
if time_resample(end)~=load(end,1)
    time_resample = [time_resample; load(end,1)];
end

iris_template = xml2struct('jobs/Iris_active_template.feb');
iris = iris_template;
%% Find the radial positions
rho_flag = true;
kk = 1;
while rho_flag
    A = sscanf(iris.febio_spec.Mesh.Nodes.node{kk}.Text,'%f,%f,%f');
    if A(2) < 10^-6
        rho(kk) = A(1);
        kk = kk +1;
    else
        rho_flag  = false;
    end
end
%% Control
iris.febio_spec.Control.time_steps=sprintf('%d',ceil(load(end,1)/step_size));
iris.febio_spec.Control.step_size=sprintf('%d',step_size);
iris.febio_spec.Control.time_stepper.dtmax=sprintf('%d',step_size);
iris.febio_spec.Control.time_stepper.dtmin=sprintf('%d',step_size/3);
%% Mesh
iris.febio_spec.Mesh = [];
iris.febio_spec.MeshDomains = [];
iris.febio_spec.MeshData = [];

iris.febio_spec.Mesh.Attributes.from='jobs/Iris_active_template.feb';
iris.febio_spec.MeshDomains.Attributes.from='jobs/Iris_active_template.feb';
iris.febio_spec.MeshData.Attributes.from='jobs/Iris_active_template.feb';

%% Material: Iris material
iris.febio_spec.Material.material.solid{1, 1}.Attributes.type='reactive viscoelastic';
iris.febio_spec.Material.material.solid{1, 1}.elastic.Attributes.type='neo-Hookean';
iris.febio_spec.Material.material.solid{1, 1}.elastic.E.Text = sprintf('%.20f',E);
iris.febio_spec.Material.material.solid{1, 1}.elastic.v.Text = sprintf('%.20f',v);

iris.febio_spec.Material.material.solid{1, 1}.bond.Attributes.type='neo-Hookean';
iris.febio_spec.Material.material.solid{1, 1}.bond.E.Text = sprintf('%.20f',E);
iris.febio_spec.Material.material.solid{1, 1}.bond.v.Text = sprintf('%.20f',v);

iris.febio_spec.Material.material.solid{1, 1}.relaxation = [];
iris.febio_spec.Material.material.solid{1, 1}.relaxation.Attributes.type='relaxation-power';
iris.febio_spec.Material.material.solid{1, 1}.relaxation.tau.Text=sprintf('%.20f',tau);
iris.febio_spec.Material.material.solid{1, 1}.relaxation.beta.Text=sprintf('%.20f',beta);

iris.febio_spec.Material.material.solid{1, 2}.T0.Attributes.type = 'math';

%r_Sphincter = 3.4; a_Sphincter = 1;
iris.febio_spec.Material.material.solid{1, 2}.T0.Text = sprintf('H(3.4+%f-(X^2+Y^2)^.5)',a_Sphincter);
%% Output Logfile
% circle.febio_spec.Output=[];
% circle.febio_spec.Output.Attributes.from='output_format.feb';
%% LoadData
for i=1:size(load,1)
    iris.febio_spec.LoadData.load_controller{1, 1}.points.point{i}.Text ...
        = sprintf ('%.20f,%.20f',[load(i,1),T_Sphicter * load(i,2)]);

    iris.febio_spec.LoadData.load_controller{1, 2}.points.point{i}.Text ...
           = sprintf ('%.20f,%.20f',[load(i,1),T_Diallator * load(i,3)]);
end
%% Plot and log file setup (optional & done through the template)
% <var type=?parameter[?fem.material[0].E?]=E?/>
%% Write file
uid = sprintf('%.0f',10^6*cputime);
uid = [uid(randperm(length(uid))),uid(randperm(length(uid)))];
uid_feb = [uid,'.feb'];
uid_log = [uid,'.log'];
uid_xplt = [uid,'.xplt'];

file_name = uid;
struct2xml(iris,file_name);
fclose('all');
movefile([file_name,'.xml'],uid_feb);
%% Run FEBio 3
% FEBio_Path = ':/Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS';
% current_path = getenv('PATH');
% flag_FEBio_path = strfind(current_path,FEBio_Path);
% if isempty(flag_FEBio_path)
%     setenv('PATH',[current_path,FEBio_Path]);
% end
system(['febio3 -i ',uid_feb,' -silent -nosplash']);
%% Visualize the displacement
%% Read output
file=fopen(uid_log);
scan_log=textscan(file,'%s','delimiter','\n','Whitespace','');
scan_log=scan_log{1,1};
fclose(file);
kk=1;
%read the fields, not that if kk=kk+1 was set on Time, there were redundant
%Time readings, which could be deleted with unique()
if strfind(scan_log{end-1},' N O R M A L   T E R M I N A T I O N')
    for i=1:length(scan_log)
        if strfind(scan_log{i},'Number of solid elements')
            %make sure you have the tab character in the beginning
            num_node = sscanf(scan_log{i-1},'	Number of nodes ................................ : %f');                      
            num_elem = sscanf(scan_log{i},'	Number of solid elements ....................... : %f');
        end
        if strfind(scan_log{i},'Data = x;y;z;sx;sy;sz;sxy;syz;sxz;Ex;Ey;Ez;Exy;Eyz;Exz')
            time(kk) = sscanf(scan_log{i-1},'Time = %f');                            
             for j = 1:num_elem           
                temp = sscanf(scan_log{i+j},'%d %f %f %f %f %f %f %f %f %f');
                    element_x(j,kk)  = temp(2);
                    element_y(j,kk)  = temp(3);                
                    element_z(j,kk)  = temp(4); 
                    
                    element_sx(j,kk)  = temp(5); 
                    element_sy(j,kk)  = temp(6);                
                    element_sz(j,kk)  = temp(7);
                    
                    element_sxy(j,kk)  = temp(8); 
                    element_syz(j,kk)  = temp(9);                
                    element_sxz(j,kk)  = temp(10);
                    
                    element_Ex(j,kk)  = temp(11); 
                    element_Ey(j,kk)  = temp(12);                
                    element_Ez(j,kk)  = temp(13); 
                    
                    element_Exy(j,kk)  = temp(14); 
                    element_Eyz(j,kk)  = temp(15);                
                    element_Exz(j,kk)  = temp(16);                     
                    
             end
            for jj = 1:num_node           
                temp = sscanf(scan_log{i+(j+6)+jj},'%d %f %f %f %f %f %f');
                    node_x(jj,kk)  = temp(2);
                    node_y(jj,kk)  = temp(3);                
                    node_z(jj,kk)  = temp(4); 
                    
                    node_ux(jj,kk)  = temp(5); 
                    node_uy(jj,kk)  = temp(6);                
                    node_uz(jj,kk)  = temp(7);
            end
         kk = kk+1;
        end
    end
    %%
    detail.time = time;

    detail.element_x = element_x;
    detail.element_y = element_y;
    detail.element_z = element_z;

    detail.element_sx = element_sx;
    detail.element_sy = element_sy;
    detail.element_sz = element_sz;

    detail.element_sxy = element_sxy;
    detail.element_syz = element_syz;
    detail.element_sxz = element_sxz;
    
    detail.element_Ex = element_Ex;
    detail.element_Ey = element_Ey;
    detail.element_Ez = element_Ez;
    
    detail.element_Exy = element_Exy;
    detail.element_Eyz = element_Eyz;
    detail.element_Exz = element_Exz;
    
    detail.node_x = node_x;
    detail.node_y = node_y;
    detail.node_z = node_z;
    
    detail.node_ux = node_ux;
    detail.node_uy = node_uy;
    detail.node_uz = node_uz;    
    
    detail.rho = rho; %the radial distance of each element
        
    dia_pupil = 2*node_x(1,:);
    e_r_pupil = 1/2*((dia_pupil/(2*min(rho))).^2-1);%lagrangian strain of the pupil
    
    
    x_spatial_element = element_x(1:(length(rho)-1),end);
    e_xx_spatial_element = element_Ex(1:(length(rho)-1),end);
    e_xx_spatial = interp1(x_spatial_element,e_xx_spatial_element,rho,'linear','extrap');
 
    %FEBio has a bug that doesn't write the zero time into log file (im not
    %sure if there is a way around it, but i deal with it like this)
    if time(1) ~= 0
        time = [0, time];
        e_r_pupil = [0, e_r_pupil];
    end
    e_r_pupil = interp1(time, e_r_pupil,time_resample,'linear');    
else
    dia_pupil   = time_resample-time_resample;
    time        = time_resample-time_resample;
    e_r_pupil   = time_resample-time_resample;
    e_xx_spatial = rho-rho;
    detail      = {};
end
    
%% Clean up the mess
if cleanup
    delete(uid_feb);delete(uid_log);delete(uid_xplt);
else
    if isfile(uid_feb)
        movefile(uid_feb,'Iris_active.feb');
    end
    if isfile(uid_log)
        movefile(uid_log,'Iris_active.log');
    end
    if isfile(uid_xplt)
        movefile(uid_xplt,'Iris_active.xplt');
    end
end
%% Outputs
if nOutputs == 2
    varargout = cell(1,nOutputs-1);
    varargout{1} = time_resample;
elseif nOutputs == 3
    varargout = cell(1,nOutputs-1);
    varargout{1} = time_resample;
    varargout{2} = detail;  
elseif nOutputs>2
    error('Too many output variables.')
end
end