function [AS] = f_as_IPGTTs( theta ,func_mex_model,D, utility , time_points_IPGTTs)
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
       strcmp(modelName,'MLMM_final_8')  ||...
       strcmp(modelName,'MLMM_extended_v0')  ||... 
       strcmp(modelName,'MLMM_extended_v0_1') 
   
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


%try   
%% Steady state simulation
ic_ss  = ic0;
p_ss   = p;

p_ss( ismember(pNames, 'ss_x')) = 0 ;

ic_ss(id_FFM)   = 1; 
ic_ss(id_FM)    = 1; 
p_ss(id_FFMname)= 1; 
p_ss(id_FMname) = 1;    % no effect of Weight on tissue-level and a i s- models
p_ss(end)       = 1000; % no IPGTT

ss  = func_mex_model(0:1:300,ic_ss, p_ss); % steady state simulation
ic0 = ss.statevalues(end,:); % new ic



%% Setting up parameters and time-vectors
                 % Time simulate IPGTT. 


%%
% fat and lean mass ratio
fmr = 0.11; lmr = 0.89; % Based on data;

p(end)    = 1000;   

p(id_FFMname)=lmr*D.Hansson2018.BW.Mean(1,1); 
p(id_FMname) =fmr*D.Hansson2018.BW.Mean(1,1);
ic0(id_FFM)=lmr*D.Hansson2018.BW.Mean(1,1); 
ic0(id_FM)= fmr*D.Hansson2018.BW.Mean(1,1);

% 2w chow
pchow    = p;
pchow(1) = 10.^theta(5);   % EI   theta(5) - Chow EI
pchow(2) = 0.7;            % carb fraction 
pchow(16)= 10.^theta(1);   % PA for chow

% 2w hfd
phfd=p;     
phfd(1)= 16.5217;       
phfd(2)=0.2;      
phfd(16)= 10.^theta(2);    % PA for hfd



%%  Simulations                                           
t_diet    = time_points_IPGTTs(end); % (days) 

sim_res = 0.1;

%2w chow
diet_chow             = func_mex_model( (0:sim_res:t_diet),ic0, pchow);
%2w hfd   
diet_hfd              = func_mex_model((0:sim_res:t_diet),ic0, phfd );  
%time_points_IPGTTs = [ 1 7 14 21 28 35 42 49 56];

for kk = 1:length(time_points_IPGTTs)
t_fasting = [  time_points_IPGTTs(kk)   time_points_IPGTTs(kk) + 0.5833 ]  ;      % 14 hours of fasting
t_IPGTT   = (0:1:120)/(60*24) + t_fasting(end); % 120 min 
pchow(end)    = t_fasting(end);   
phfd(end)     = t_fasting(end);   

ss_chow    = func_mex_model( t_fasting, diet_chow.statevalues( time_points_IPGTTs(kk)/sim_res,:), pchow);  % chow

tmp = strcat('AS.IPGTT_chow', int2str(time_points_IPGTTs(kk)),'d' );

eval( strcat ( tmp , '= func_mex_model( t_IPGTT , ss_chow.statevalues(end,:), pchow );' ) ) ;
%2w hfd   
ss_hfd          = func_mex_model( t_fasting, diet_hfd.statevalues( time_points_IPGTTs(kk)/sim_res,:), phfd);  % hfd

tmp = strcat('AS.IPGTT_hfd', int2str(time_points_IPGTTs(kk)),'d' );

eval( strcat ( tmp , '= func_mex_model( t_IPGTT , ss_hfd.statevalues(end,:), phfd );') ) ;


end

%%
AS.ss                  = ss;
AS.diet_chow                  = diet_chow;
AS.diet_hfd                   = diet_hfd;

end

