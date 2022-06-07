function [AllSimulations] = f_as_invitro_new_2( theta ,func_mex_model,utility)
%% New SIMULATIONS FUNCTION MAINLY FOR PLOTTING.

ic0       = utility.ic0;
objModel  = utility.objModel; % variables used for indexing
modelName = utility.modelName;
pNames    = utility.pNames ;
p         = utility.p0; % Parmeter vector with values from model files. - Vector used for simulation


id_FFM     = ismember(SBstates(objModel),'FFM') ; 
id_FM      = ismember(SBstates(objModel),'FM');
id_FFMname = ismember(pNames,'FFMinit') ; 
id_FMname  = ismember(pNames,'FMinit');

% Chaning values in vector P, for paramters being optimized. 
if strcmp(modelName,'MLMM_final')
    % parameter and time vector setup
    p(17:37)=[10.^theta(8:19) theta(20)  10.^theta(21:end)];  % GSS to ik1
    if p(ismember(pNames,'ik1')) == 10^theta(end) % bounds check
    else
        disp('check parameter vector')
    end
    
elseif strcmp(modelName,'MLMM_final_2')  ||...
       strcmp(modelName,'MLMM_final_4')  ||...
       strcmp(modelName,'MLMM_final_5')  ||...
       strcmp(modelName,'MLMM_final_6')  ||...
       strcmp(modelName,'MLMM_final_7')  ||...
       strcmp(modelName,'MLMM_final_8') 
   
    
    p(17:36)=[10.^theta(8:18) theta(19)  10.^theta(20:end)];  % GSS to ik1
    
    if p(ismember(pNames,'ik1')) == 10^theta(end) % bounds check
    else
        disp('check parameter vector')
    end
    
elseif strcmp(modelName,'MLMM_final_3')
    p(17:36)=[10.^theta(8:18) theta(19)  10.^theta(20:end)];  % GSS to ik1
    
    if p(ismember(pNames,'ik1')) == 10^theta(end) % bounds check
    else
        disp('check parameter vector')
    end    
    
end



try   
%% Steady state simulation
ic_ss  = ic0;
p_ss   = p;

p_ss( ismember(pNames, 'ss_x')) = 0 ;

ic_ss(id_FFM)   = 1; 
ic_ss(id_FM)    = 1; 
p_ss(id_FFMname)= 1; 
p_ss(id_FMname) = 1;    % no effect of Weight on tissue-level and a i s- models
p_ss(end)       = 1000; % no IPGTT

ss  = func_mex_model(0:100:300,ic_ss, p_ss); % steady state simulation
ic0 = ss.statevalues(end,:); % new ic


%% Hansson et al (2019) 4w - invitro insulin stimulation of primary adipocytes
% Setting up parameters and time-vectors

t_diet     = 28; % 4 weeks
t_fasting  = [ t_diet  t_diet+0.5833 ]  ;  % 14 hours fasting

t_in_vitro   = (0:1:30)/(60*24) + t_fasting(end); % 120 min 

p(end)     =  t_fasting(end) + 1000; % no IPGTT



pchow_vitro = p; 
pchow_vitro(1)  = 10.^theta(5); 
pchow_vitro(2)  = 0.7; 
pchow_vitro(16) = 10.^theta(1); % using same BW parameters as Hanson et al 2018 simulations

phfd_vitro=p;  
phfd_vitro(1)  = 16.5217; 
phfd_vitro(2)  = 0.2; 
phfd_vitro(16) = 10.^theta(2); % using same BW parameters as Hanson et al 2018 simulations
% Here we use Ic0 directley which use a start BW of 23g lean and fat mass ratio 89/11 %

%%% simulations 
chow4w       = func_mex_model( (0:1:t_diet),ic0, pchow_vitro);

pchow_vitro(ismember(pNames,'Insact'))=0; % WASH - remvoing effect of I_E in vivo insulin

ss_chow4w    = func_mex_model( t_fasting, chow4w.statevalues(end,:), pchow_vitro);  

pchow_vitro(ismember(pNames,'ins')) = 0.1; % nM , experiment insulin dose

vitro_chow4w = func_mex_model( t_in_vitro , ss_chow4w.statevalues(end,:), pchow_vitro ); 


%%%

hfd4w       = func_mex_model( (0:1:t_diet),ic0, phfd_vitro );  

phfd_vitro(ismember(pNames,'Insact'))=0; 

ss_hfd4w    = func_mex_model( t_fasting, hfd4w.statevalues(end,:), phfd_vitro); 

phfd_vitro(ismember(pNames,'ins'))=0.1; 

vitro_hfd4w = func_mex_model( t_in_vitro , ss_hfd4w.statevalues(end,:), phfd_vitro ); 

catch err 
    disp(err)
end

%%
AllSimulations.chow4w              = chow4w;
AllSimulations.hfd4w               = hfd4w;
AllSimulations.vitro_chow4w        = vitro_chow4w;
AllSimulations.vitro_hfd4w         = vitro_hfd4w;


end

