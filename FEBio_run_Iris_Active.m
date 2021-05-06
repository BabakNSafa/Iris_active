function [e_r_pupil,varargout] = FEBio_run_Iris_Active(x,time_resample,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function reads the jobs/Iris_active.feb and makes
%changes to it
%
%function [dia_pupil,varargout] = FEBio_run_Iris_Active(T,varargin)
% Inputs 
%     x (1x5) = [E, v, tau T_Sphicter, T_Diallator]
%     E                 = Matrix modulus
%     v                 = Matrix Poisson's ratio
%     tau               = Matrix kinetics time-constant (sec)
%     T_Sphicter        = Sphicter muscle stress as a function of time (in units of E)
%     T_Diallator       = Diallator muscle stress as a function of time (in units of E)
%
%     time_resample     = Time (sec)
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
E   = x(:,1);
v   = x(:,2);
tau = x(:,3);
T_Sphicter  = x(:,4);
T_Diallator = x(:,5);

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
iris.febio_spec.Control.time_steps=sprintf('%d',max(time_resample)*10);
iris.febio_spec.Control.step_size=sprintf('%d',.1);
%% Material: Iris material
% iris.febio_spec.Material.material.solid{1, 1}.kinetics = '1';
% iris.febio_spec.Material.material.solid{1, 1}.trigger = '0';
% 
% iris.febio_spec.Material.material.solid{1, 1} = ...
%     rmfield(iris.febio_spec.Material.material.solid{1, 1},'E');
% iris.febio_spec.Material.material.solid{1, 1} = ...
%     rmfield(iris.febio_spec.Material.material.solid{1, 1},'v');

iris.febio_spec.Material.material.solid{1, 1}.Attributes.type='reactive viscoelastic';
iris.febio_spec.Material.material.solid{1, 1}.elastic.Attributes.type='neo-Hookean';
iris.febio_spec.Material.material.solid{1, 1}.elastic.E.Text = sprintf('%.20f',E);
iris.febio_spec.Material.material.solid{1, 1}.elastic.v.Text = sprintf('%.20f',v);

iris.febio_spec.Material.material.solid{1, 1}.bond.Attributes.type='neo-Hookean';
iris.febio_spec.Material.material.solid{1, 1}.bond.E.Text = sprintf('%.20f',E);
iris.febio_spec.Material.material.solid{1, 1}.bond.v.Text = sprintf('%.20f',v);

iris.febio_spec.Material.material.solid{1, 1}.relaxation.Attributes.type='relaxation-exponential';
iris.febio_spec.Material.material.solid{1, 1}.relaxation.tau.Text=sprintf('%.20f',tau);


iris.febio_spec.Material.material.solid{1, 2}.T0.Attributes.type = 'math';
% iris.febio_spec.Material.material.solid{1, 2}.T0.Text = '(1-H((X^2+Y^2)^.5-4))';
iris.febio_spec.Material.material.solid{1, 2}.T0.Text = '(1-tanh(5*((X^2+Y^2)^.5-4)))*H((X^2+Y^2)^.5-3)';

iris.febio_spec.Material.material.solid{1, 3}.T0.Text = '1';
%% Output Logfile
% circle.febio_spec.Output=[];
% circle.febio_spec.Output.Attributes.from='output_format.feb';
%% LoadData
iris.febio_spec.LoadData.load_controller{1, 2}.points.point{2}.Text ...
    = sprintf ('%.20f,%.20f',[20,T_Sphicter]);

iris.febio_spec.LoadData.load_controller{1, 3}.points.point{2}.Text ...
       = sprintf ('%.20f,%.20f',[30,T_Diallator]);
%% Plot and log file setup (optional & done through the template)
% <var type=?parameter[?fem.material[0].E?]=E?/>
%% Write file
uid = sprintf('%.0f',10^6*cputime);
uid = uid(randperm(length(uid)));
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
        if strfind(scan_log{i},'Data = x;y;z;sx;sy;sz;Ex;Ey')
            time(kk) = sscanf(scan_log{i-1},'Time = %f');                            
             for j = 1:num_elem           
                temp = sscanf(scan_log{i+j},'%d %f %f %f %f %f %f %f %f %f');
                    element_x(j,kk)  = temp(2);
                    element_y(j,kk)  = temp(3);                
                    element_z(j,kk)  = temp(4); 
                    
                    element_sx(j,kk)  = temp(5); 
                    element_sy(j,kk)  = temp(6);                
                    element_sz(j,kk)  = temp(7);
                    
                    element_Ex(j,kk)  = temp(8); 
                    element_Ey(j,kk)  = temp(9);                
                    element_Ez(j,kk)  = temp(10); 
                    
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

    detail.element_Ex = element_Ex;
    detail.element_Ey = element_Ey;
    detail.element_Ez = element_Ez;
    
    detail.node_x = node_x;
    detail.node_y = node_y;
    detail.node_z = node_z;
    
    detail.node_ux = node_ux;
    detail.node_uy = node_uy;
    detail.node_uz = node_uz;    
    
    detail.rho = rho; %the radial distance of each element
        
    dia_pupil = 2*node_x(1,:);
    
    e_r_pupil = 1/2*((dia_pupil/(2*min(rho))).^2-1);%lagrangian strain of the pupil
    e_r_pupil = interp1(time, e_r_pupil,time_resample,'pchip');    
else
    dia_pupil   = time_resample-time_resample;
    time        = time_resample-time_resample;
    e_r_pupil   = time_resample-time_resample;
    detail      = {};
end

%% Clean up the mess
if cleanup
    delete(uid_feb);delete(uid_log);delete(uid_xplt);
else
    movefile(uid_feb,'Iris_active.feb');
    movefile(uid_log,'Iris_active.log');
    movefile(uid_xplt,'Iris_active.xplt');
end
%% Outputs
if nOutputs == 2
    varargout = cell(1,nOutputs-1);
    varargout{1} = time;
elseif nOutputs == 3
    varargout = cell(1,nOutputs-1);
    varargout{1} = time;
    varargout{2} = detail;  
elseif nOutputs>2
    error('Too many output variables.')
end
end