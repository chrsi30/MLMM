function [AllSimulations] = Function_AllSimulations_t2d_prog_invitro_exp( theta ,func_mex_model, D, utility)
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

%try   
     
%% Steady state simulation
ic_ss = ic0; p_ss = p; p_ss( ismember(pNames, 'ss_x')) = 0 ;
ic_ss(id_FFM)=1 ; ic_ss(id_FM)= 1; p_ss(id_FFMname)=1; p_ss(id_FMname)=1; % no effect of Weight on tissue-level and a i s- models
p_ss(end) = 10000; %no IPGTT
ss=func_mex_model(0:100:300,ic_ss, p_ss); % steady state simulation
cost_ss_gcc = (( ss.variablevalues(end, end-1)  -  D.IPGTTgcc.Mean(1,1)  )^2 )./( D.IPGTTgcc.SEM(1,1) )^2 ;% basal glucose
cost_ss_icc = (( ss.variablevalues(end, end)  -  D.IPGTTicc.Mean(1,1)  )^2 )./( D.IPGTTicc.SEM(1,1) )^2   ;% basal insulin
ss_cost = cost_ss_gcc + cost_ss_icc;
ic0 = ss.statevalues(end,:); % overwrite ic0
%% 4 weeks -chow -hfd - t2d hfd
t_diet=28; 
p(end) = 10000; % no IPGTT

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

pchow(ismember(pNames,'Insact'))=0; % remvoing effect of I_E in vivo insulin
phfd(ismember(pNames,'Insact'))=0; 
psc1(ismember(pNames,'Insact'))=0; 

pchow(ismember(pNames,'ins'))=100; % nM , experiment insulin dose
phfd(ismember(pNames,'ins'))=100; 
psc1(ismember(pNames,'ins'))=100; 

t_inc=(0:1:30)/(60*24)+ 1 ;
% time of IPGTT (triggers after fasting)
invitro_chow0=func_mex_model( t_inc , chow4w.statevalues(1,:), pchow ); 
invitro_hfd0=func_mex_model( t_inc , hfd4w.statevalues(1,:), phfd ); 
invitro_hfd0t2d=func_mex_model( t_inc , hfd4wt2d.statevalues(1,:), psc1 );      
% 1 week
t_inc=(0:1:30)/(60*24)+ 7 ;
invitro_chow1w=func_mex_model( t_inc , chow4w.statevalues(71,:), pchow ); 
invitro_hfd1w=func_mex_model( t_inc , hfd4w.statevalues(71,:), phfd ); 
invitro_hfd1wt2d=func_mex_model( t_inc , hfd4wt2d.statevalues(71,:), psc1 );      
% 2 week 
t_inc=(0:1:30)/(60*24)+ 14 ;
invitro_chow2w=func_mex_model( t_inc , chow4w.statevalues(141,:), pchow ); 
invitro_hfd2w=func_mex_model( t_inc , hfd4w.statevalues(141,:), phfd ); 
invitro_hfd2wt2d=func_mex_model( t_inc , hfd4wt2d.statevalues(141,:), psc1 );      
% 3 week
t_inc=(0:1:30)/(60*24)+ 21 ;
invitro_chow3w=func_mex_model( t_inc , chow4w.statevalues(211,:), pchow ); 
invitro_hfd3w=func_mex_model( t_inc , hfd4w.statevalues(211,:), phfd ); 
invitro_hfd3wt2d=func_mex_model( t_inc , hfd4wt2d.statevalues(211,:), psc1 );      
%4 week
t_inc=(0:1:30)/(60*24)+ 28 ;
invitro_chow4w=func_mex_model( t_inc , chow4w.statevalues(281,:), pchow ); 
invitro_hfd4w=func_mex_model( t_inc , hfd4w.statevalues(281,:), phfd ); 
invitro_hfd4wt2d=func_mex_model( t_inc , hfd4wt2d.statevalues(281,:), psc1 );    

% catch err 
%     disp(err)
% end

%%
AllSimulations.chow4w             = chow4w;
AllSimulations.hfd4w              = hfd4w;
AllSimulations.hfd4wt2d           = hfd4wt2d;

AllSimulations.invitro_chow0       = invitro_chow0;
AllSimulations.invitro_hfd0        = invitro_hfd0;
AllSimulations.invitro_hfd0t2d     = invitro_hfd0t2d;

AllSimulations.invitro_chow1w       = invitro_chow1w;
AllSimulations.invitro_hfd1w        = invitro_hfd1w;
AllSimulations.invitro_hfd1wt2d     = invitro_hfd1wt2d;

AllSimulations.invitro_chow2w       = invitro_chow2w;
AllSimulations.invitro_hfd2w        = invitro_hfd2w;
AllSimulations.invitro_hfd2wt2d     = invitro_hfd2wt2d;

AllSimulations.invitro_chow3w       = invitro_chow3w;
AllSimulations.invitro_hfd3w        = invitro_hfd3w;
AllSimulations.invitro_hfd3wt2d     = invitro_hfd3wt2d;

AllSimulations.invitro_chow4w       = invitro_chow4w;
AllSimulations.invitro_hfd4w        = invitro_hfd4w;
AllSimulations.invitro_hfd4wt2d     = invitro_hfd4wt2d;

end

