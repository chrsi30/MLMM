function [AllSimulations] = Function_AllSimulations_t2d( theta ,func_mex_model, D, utility)
ic0 = utility.ic0;
p = utility.p0;
objModel = utility.objModel; % variables used for indexing

pNames = utility.pNames ; 
id_FFM = ismember(SBstates(objModel),'FFM') ; id_FM = ismember(SBstates(objModel),'FM');
id_FFMname = ismember(pNames,'FFMinit') ; id_FMname = ismember(pNames,'FMinit');
% parameter and time vector setup
p(17:37)=[10.^theta(8:19) theta(20)  10.^theta(21:end)];  % GSS to ik1
if p(ismember(pNames,'ik1')) == 10^theta(end) % bounds check
else
    disp('check parameter vector')
end
try   
    
    
%% Steady state simulation
ic_ss = ic0; p_ss = p; p_ss( ismember(pNames, 'ss_x')) = 0 ;
ic_ss(id_FFM)=1 ; ic_ss(id_FM)= 1; p_ss(id_FFMname)=1; p_ss(id_FMname)=1; % no effect of Weight on tissue-level and a i s- models
p_ss(end) = 1000; %no IPGTT
ss=func_mex_model(0:100:300,ic_ss, p_ss); % steady state simulation
cost_ss_gcc = (( ss.variablevalues(end, end-1)  -  D.IPGTTgcc.Mean(1,1)  )^2 )./( D.IPGTTgcc.SEM(1,1) )^2 ;% basal glucose
cost_ss_icc = (( ss.variablevalues(end, end)  -  D.IPGTTicc.Mean(1,1)  )^2 )./( D.IPGTTicc.SEM(1,1) )^2   ;% basal insulin
ss_cost = cost_ss_gcc + cost_ss_icc;
ic0 = ss.statevalues(end,:); % overwrite ic0
%% 4 weeks -chow -hfd - t2d hfd
t_diet=28; 
p(end) = 1000; 

%chow
pchow=p;
pchow(1)=10.^theta(5);   % EI   theta(5) - Chow EI
pchow(2)=0.7;            % carb fraction 
pchow(16)= 10.^theta(1); % PA for chow
pchow(id_FFMname)=0.89*D.BW.mean(1,1); pchow(id_FMname)=0.11*D.BW.mean(1,1);
ic_chow = ic0; ic_chow(id_FFM)=0.89*D.BW.mean(1,1) ; ic_chow(id_FM)= 0.11*D.BW.mean(1,1);
%hfd
phfd=p;     
phfd(1)= 16.5217;       
phfd(2)=0.2;      
phfd(16)      = 10.^theta(2);   
phfd(id_FFMname)=0.89*D.BW.mean(4,1);    
phfd(id_FMname)=0.11*D.BW.mean(4,1);
ic_hfd = ic0; ic_hfd(id_FFM)=0.89*D.BW.mean(4,1) ; ic_hfd(id_FM)= 0.11*D.BW.mean(4,1);
%Scenario 1 
psc1=p;     
psc1(1)= 16.5217;       
psc1(2)=0.2;      
psc1(16)      = 10.^theta(2);   
psc1(id_FFMname)=0.89*D.BW.mean(4,1);    psc1(id_FMname)=0.11*D.BW.mean(4,1);
psc1( ismember( pNames, 't2dswitch'))=1;
ic_sc1 = ic0; ic_sc1(id_FFM)=0.89*D.BW.mean(4,1) ; ic_sc1(id_FM)= 0.11*D.BW.mean(4,1);

chow4w=func_mex_model((0:0.1:t_diet),ic_chow, pchow );   
hfd4w=func_mex_model((0:0.1:t_diet),ic_hfd, phfd );   
hfd4wt2d=func_mex_model((0:0.1:t_diet),ic_sc1, psc1 ); 

%  1    71   141   211   281 index based on time vector
% inital
pchow(end) = 1.001; 
phfd(end) = 1.001; 
psc1(end) = 1.001; 
t_IPGTT=(0:1:300)/(60*24)+ 1 ;
% time of IPGTT (triggers after fasting)
IPGTT_chow0=func_mex_model( t_IPGTT , chow4w.statevalues(1,:), pchow ); 
IPGTT_hfd0=func_mex_model( t_IPGTT , hfd4w.statevalues(1,:), phfd ); 
IPGTT_hfd0t2d=func_mex_model( t_IPGTT , hfd4wt2d.statevalues(1,:), psc1 );      
% 1 week
pchow(end) = 7.001; 
phfd(end) = 7.001; 
psc1(end) = 7.001; 
t_IPGTT=(0:1:300)/(60*24)+ 7 ;
IPGTT_chow1w=func_mex_model( t_IPGTT , chow4w.statevalues(71,:), pchow ); 
IPGTT_hfd1w=func_mex_model( t_IPGTT , hfd4w.statevalues(71,:), phfd ); 
IPGTT_hfd1wt2d=func_mex_model( t_IPGTT , hfd4wt2d.statevalues(71,:), psc1 );      
% 2 week 
pchow(end) = 14.001; 
phfd(end) = 14.001; 
psc1(end) = 14.001; 
t_IPGTT=(0:1:300)/(60*24)+ 14 ;
IPGTT_chow2w=func_mex_model( t_IPGTT , chow4w.statevalues(141,:), pchow ); 
IPGTT_hfd2w=func_mex_model( t_IPGTT , hfd4w.statevalues(141,:), phfd ); 
IPGTT_hfd2wt2d=func_mex_model( t_IPGTT , hfd4wt2d.statevalues(141,:), psc1 );      
% 3 week
pchow(end) = 21.001; 
phfd(end) = 21.001; 
psc1(end) = 21.001; 
t_IPGTT=(0:1:300)/(60*24)+ 21 ;
IPGTT_chow3w=func_mex_model( t_IPGTT , chow4w.statevalues(211,:), pchow ); 
IPGTT_hfd3w=func_mex_model( t_IPGTT , hfd4w.statevalues(211,:), phfd ); 
IPGTT_hfd3wt2d=func_mex_model( t_IPGTT , hfd4wt2d.statevalues(211,:), psc1 );      
%4 week
pchow(end) = 28.001; 
phfd(end) = 28.001; 
psc1(end) = 28.001; 
t_IPGTT=(0:1:300)/(60*24)+ 28 ;
IPGTT_chow4w=func_mex_model( t_IPGTT , chow4w.statevalues(281,:), pchow ); 
IPGTT_hfd4w=func_mex_model( t_IPGTT , hfd4w.statevalues(281,:), phfd ); 
IPGTT_hfd4wt2d=func_mex_model( t_IPGTT , hfd4wt2d.statevalues(281,:), psc1 );    

catch err 
    disp(err)
end

%%
AllSimulations.chow4w             = chow4w;
AllSimulations.hfd4w              = hfd4w;
AllSimulations.hfd4wt2d           = hfd4wt2d;

AllSimulations.IPGTT_chow0       = IPGTT_chow0;
AllSimulations.IPGTT_hfd0        = IPGTT_hfd0;
AllSimulations.IPGTT_hfd0t2d     = IPGTT_hfd0t2d;

AllSimulations.IPGTT_chow1w       = IPGTT_chow1w;
AllSimulations.IPGTT_hfd1w        = IPGTT_hfd1w;
AllSimulations.IPGTT_hfd1wt2d     = IPGTT_hfd1wt2d;

AllSimulations.IPGTT_chow2w       = IPGTT_chow2w;
AllSimulations.IPGTT_hfd2w        = IPGTT_hfd2w;
AllSimulations.IPGTT_hfd2wt2d     = IPGTT_hfd2wt2d;

AllSimulations.IPGTT_chow3w       = IPGTT_chow3w;
AllSimulations.IPGTT_hfd3w        = IPGTT_hfd3w;
AllSimulations.IPGTT_hfd3wt2d     = IPGTT_hfd3wt2d;

AllSimulations.IPGTT_chow4w       = IPGTT_chow4w;
AllSimulations.IPGTT_hfd4w        = IPGTT_hfd4w;
AllSimulations.IPGTT_hfd4wt2d     = IPGTT_hfd4wt2d;

end

