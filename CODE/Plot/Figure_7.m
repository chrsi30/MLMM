%% Simulations for best cost.
load('./Parameters/theta_best.mat'); theta = Results.xbest; 

%%
id_FFM     = ismember(SBstates(objModel),'FFM') ;
id_FM      = ismember(SBstates(objModel),'FM');
id_FFMname = ismember(pNames,'FFMinit') ; 
id_FMname  = ismember(pNames,'FMinit');

if  strcmp(modelName,'MLMM_extended_v0_1') 
   
    p(17:36)=[10.^theta(8:18) theta(19)  10.^theta(20:end)];  % GSS to ik1
    
    if p(ismember(pNames,'ik1')) == 10^theta(end) % bounds check
    else
        disp('check parameter vector')
    end
          
end
%% Steady state simulation
ic_ss  = ic0;
p_ss   = p;

p_ss( ismember(pNames, 'ss_x')) = 0 ;

ic_ss(id_FFM)   = 1; 
ic_ss(id_FM)    = 1; 
p_ss(id_FFMname)= 1; 
p_ss(id_FMname) = 1;    % no effect of Weight on tissue-level and a i s- models
p_ss(end)       = 1000; % no IPGTT

ss  = func_mex_model(0:1:500,ic_ss, p_ss); % steady state simulation
ic0 = ss.statevalues(end,:); % new ic

%% Setting up parameters and time-vectors
time_points_IPGTTs = 1:1:7*10;

%%
p(ismember(pNames,'xxh')) = 1;
p(ismember(pNames,'r1'))  = 0.3948; 
p(ismember(pNames,'d0'))  = 0.04604149;
p(ismember(pNames,'r2'))  = 0.78 ; 

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
% chow
pchow(ismember(pNames,'xxh')) = 1;
%pchow(ismember(pNames,'xxh')) = 0;
diet_chow             = func_mex_model( (0:sim_res:t_diet),ic0, pchow);

% hfd   - normal
diet_hfd              = func_mex_model((0:sim_res:t_diet),ic0, phfd );  

% hfd   - normal
phfd2 = phfd;
phfd2(ismember(pNames,'r1'))  = 0.3948; 
phfd2(ismember(pNames,'d0'))  = 0.04604149;
phfd2(ismember(pNames,'r2'))  = 1.05  ; 
diet_hfd2             = func_mex_model((0:sim_res:t_diet),ic0, phfd2 );  


% hfd   - humanized?
phfd3 = phfd;
phfd3(ismember(pNames,'r1'))  = 0.55; 
phfd3(ismember(pNames,'d0'))  = 0.029;
phfd3(ismember(pNames,'r2'))  = 1.5  ; 
diet_hfd3             = func_mex_model((0:sim_res:t_diet),ic0, phfd3 );  

            
%% Plotting
label_size       = 16 ;
label_size_ax    = 16 ;
marker_line_size = 1.5;
marker_size      = 5  ;
line_size        = 1.5;
title_size       = 18 ;

control_color   = [ 0,   0,     0 ]; 
HFD_color       = [ 1,   .6,   .1 ]; 
HFD_color2      = [ 1,   .3,    0 ]; 
HFD_color3      = [ 0  , 1,    0 ]; 

id_Beta = ismember(SBstates(objModel),'Beta');
id_G_C  = ismember(SBvariables(objModel),'Gconc');
id_I_C  = ismember(SBvariables(objModel),'Iconc');

id_f_IR  = ismember(SBvariables(objModel),'f_IR_Ins');
id_f_EGP = ismember(SBvariables(objModel),'f_IR_EGP');
id_EGP   = ismember(SBvariables(objModel),'PlotEGP');

fig=0;

%%
t = diet_hfd.time./7;

figure(7)

set(figure(7), 'outerposition',[0 0 2000 500], 'PaperType','a4')


subplot(1,3,1)
plot( t,  diet_hfd.statevalues(:, id_Beta ), 'color', HFD_color  , 'LineWidth',line_size);
hold on
plot( t,  diet_chow.statevalues(:, id_Beta ), 'color', control_color  , 'LineWidth',line_size);
plot( t,  diet_hfd2.statevalues(:, id_Beta ), 'color', HFD_color2  , 'LineWidth',line_size);
plot( t,  diet_hfd3.statevalues(:, id_Beta ), 'color', HFD_color3  , 'LineWidth',line_size);


ylim([ 0 3 ]); yticks([ 0  3]); ylabel('\beta-cell mass (mg)');
xlim([0 10]); xticks([ 0  10]); xlabel('Time (weeks)');


subplot(1,3,2)
plot( t,  diet_hfd.variablevalues(:, id_I_C ), 'color', HFD_color  , 'LineWidth',line_size);
hold on
plot( t,  diet_chow.variablevalues(:, id_I_C ), 'color', control_color  , 'LineWidth',line_size);
plot( t,  diet_hfd2.variablevalues(:, id_I_C ), 'color', HFD_color2  , 'LineWidth',line_size);
plot( t,  diet_hfd3.variablevalues(:, id_I_C ), 'color', HFD_color3  , 'LineWidth',line_size);


ylim([ 0 600 ]); yticks([ 0  600]); ylabel('Fasting plasma insulin levels (pmol/L)');
xlim([0 10]); xticks([ 0  10]); xlabel('Time (weeks)');

subplot(1,3,3)
plot( t,  diet_hfd.variablevalues(:, id_G_C ), 'color', HFD_color  , 'LineWidth',line_size);
hold on
plot( t,  diet_chow.variablevalues(:, id_G_C ), 'color', control_color  , 'LineWidth',line_size);
plot( t,  diet_hfd2.variablevalues(:, id_G_C ), 'color', HFD_color2  , 'LineWidth',line_size);
plot( t,  diet_hfd3.variablevalues(:, id_G_C ), 'color', HFD_color3  , 'LineWidth',line_size);


ylim([ 0 16 ]); yticks([ 0  16]); ylabel('Fasting plasma glucose levels (mmol/L)');
xlim([0 10]); xticks([ 0  10]); xlabel('Time (weeks)');

