
%%
%% UNC
id_PKB_phos473 = ismember(SBvariables(utility.objModel),'totalPKB473'); 
id_PKB         = ismember(SBvariables(utility.objModel),'totalPKB'); 
id_AS          = ismember(SBvariables(utility.objModel),'totalAS160'); 
id_As_phos     = ismember(SBstates(utility.objModel),'AS160p'); 
id_GU          = ismember(SBvariables(utility.objModel),'CL_GIA'); 

for i = 1:nfiles
load( strcat( filedir ,'/',files(i).name ) )
disp(strcat('Processing - ',files(i).name, int2str(i),'/',int2str(nfiles)))
theta=Results.xbest;

[ AS ] = f_as_new( theta ,func_mex_model, D, utility);

ss                   = AS.ss;
chow8w               = AS.chow8w;
hfd8w                = AS.hfd8w;
IPGTT_hfd6w          = AS.IPGTT_hfd6w;
IPGTT_chow6w         = AS.IPGTT_chow6w;
chow2w               = AS.chow2w;
IPGTT_chow2w         = AS.IPGTT_chow2w;
chow12dHFD2d         = AS.chow12dHFD2d;
IPGTT_chow12dHFD2d   = AS.IPGTT_chow12dHFD2d;
chow10dHFD4d         = AS.chow10dHFD4d;
IPGTT_chow10dHFD4d   = AS.IPGTT_chow10dHFD4d;
chow8dHFD6d          = AS.chow8dHFD6d;
IPGTT_chow8dHFD6d    = AS.IPGTT_chow8dHFD6d;
hfd2w                = AS.hfd2w;
IPGTT_hfd2w          = AS.IPGTT_hfd2w;

[ AS3 ] = f_as_invitro_new_2( theta ,func_mex_model, utility);

chow4w       = AS3.chow4w;            
hfd4w        = AS3.hfd4w;               
vitro_chow4w = AS3.vitro_chow4w;       
vitro_hfd4w  = AS3.vitro_hfd4w;        

% MIN and MAX value calculations
sim_id = {'IPGTT_chow2w','IPGTT_chow12dHFD2d','IPGTT_chow10dHFD4d','IPGTT_chow8dHFD6d','IPGTT_hfd2w','IPGTT_hfd6w'};
state_name = {'PKB473', 'PKB308473','AS160p','AS160','RABGTP','PKBda','PMAa'};
for kk = 1:length(sim_id)  
    for jj = 1:length(state_name)
    tmp1 = strcat('Max_', sim_id{kk}, '_', state_name{jj} );
    tmp2 = strcat('Min_', sim_id{kk}, '_', state_name{jj} );
    eval(strcat('id_',  state_name{jj},  "= ismember(SBstates(objModel),'", state_name{jj},"');"));
    eval(strcat('tmp3 = id_',  state_name{jj},';') );  
    
    if i == 1
    tmp4 = strcat( sim_id{kk} , '.statevalues(:, tmp3 );');
    eval(strcat( tmp1 ,'=' ,tmp4)); eval(strcat( tmp2 ,'=' ,tmp4));
    end
    eval(strcat('[ ',tmp2,' ,  ',tmp1,' ]     = f_calc_minmax( ',tmp2,',', tmp1, " , 'state' , ", sim_id{kk} ,' , tmp3);'))
    end
end

sim_id = {'IPGTT_chow2w','IPGTT_chow12dHFD2d','IPGTT_chow10dHFD4d','IPGTT_chow8dHFD6d','IPGTT_hfd2w','IPGTT_hfd6w'};
varb_name = {'insulin_input', 'GLUT4m'};
for kk = 1:length(sim_id)  
    for jj = 1:length(varb_name)
    tmp1 = strcat('Max_', sim_id{kk}, '_', varb_name{jj} );
    tmp2 = strcat('Min_', sim_id{kk}, '_', varb_name{jj} );
    eval(strcat('id_',  varb_name{jj},  "= ismember(SBvariables(objModel),'", varb_name{jj},"');"));
    eval(strcat('tmp3 = id_',  varb_name{jj},';') );  
    
    if i == 1
    tmp4 = strcat( sim_id{kk} , '.variablevalues(:, tmp3 );');
    eval(strcat( tmp1 ,'=' ,tmp4)); eval(strcat( tmp2 ,'=' ,tmp4));
    end
    eval(strcat('[ ',tmp2,' ,  ',tmp1,' ]     = f_calc_minmax( ',tmp2,',', tmp1, " , 'variable' , ", sim_id{kk} ,' , tmp3);'))
    end
end

sim_id2 =  {'vitro_chow4w', 'vitro_hfd4w'    }; %%string(name_2)';
state_name2 = {'totalPKB473', 'totalPKB'} ; %,'U_ii','U_id_A','U_id_M','U_id_H'};

for kk = 1:length(sim_id2)
    for jj = 1:length(state_name2)
        tmp1 = strcat('Max_', sim_id2{kk}, '_', state_name2{jj} );
        tmp2 = strcat('Min_', sim_id2{kk}, '_', state_name2{jj} );
        eval(strcat('id_',  state_name2{jj},  "= ismember(SBstates(objModel),'", state_name2{jj},"');"));
        eval(strcat('tmp3 = id_',  state_name2{jj},';') );
        
        if i == 1
            tmp4 = strcat( sim_id2{kk} , '.statevalues(:, tmp3 );');
            eval(strcat( tmp1 ,'=' ,tmp4)); eval(strcat( tmp2 ,'=' ,tmp4));              
        end
       eval(strcat('[ ',tmp2,' ,  ',tmp1,' ]     = f_calc_minmax( ',tmp2,',', tmp1, " , 'variable' , ", sim_id2{kk} ,' , tmp3);'))
    end
end

PKB_fold_chow = vitro_chow4w.variablevalues(:,id_PKB_phos473)./vitro_chow4w.variablevalues(end,id_PKB);
PKB_fold_hfd  = vitro_hfd4w.variablevalues(:,id_PKB_phos473)./vitro_hfd4w.variablevalues(end,id_PKB);
fold_diet_PKB = PKB_fold_hfd./PKB_fold_chow;% fold between diets ( fold on fold )


AS160_fold_chow = vitro_chow4w.statevalues(:,id_As_phos)./vitro_chow4w.variablevalues(end,id_AS);
AS160_fold_hfd  = vitro_hfd4w.statevalues(:,id_As_phos)./vitro_hfd4w.variablevalues(end,id_AS);
fold_diet_AS160 = AS160_fold_hfd./AS160_fold_chow;% fold between diets ( fold on fold )

AUC_chow4w_GU     = trapz(vitro_chow4w.variablevalues(:,id_GU));
AUC_hfd4w_GU      = trapz( vitro_hfd4w.variablevalues(:,id_GU));
fold_AUC_hfd4w_GU = AUC_hfd4w_GU./AUC_chow4w_GU; % fold on chow


if i == 1
    Max_vitro_hfd4w_fold_diet_PKB = fold_diet_PKB;
    Min_vitro_hfd4w_fold_diet_PKB = fold_diet_PKB;
    
    Max_vitro_hfd4w_fold_diet_AS160 = fold_diet_AS160;
    Min_vitro_hfd4w_fold_diet_AS160 = fold_diet_AS160;
    
    Max_vitro_hfd4w_fold_diet_GU = fold_AUC_hfd4w_GU;
    Min_vitro_hfd4w_fold_diet_GU = fold_AUC_hfd4w_GU;
    
    
    Max_vitro_chow4w_PKB = vitro_chow4w.variablevalues(:,id_PKB_phos473);
    Min_vitro_chow4w_PKB = vitro_chow4w.variablevalues(:,id_PKB_phos473);
    Max_vitro_hfd4w_PKB  = vitro_hfd4w.variablevalues(:,id_PKB_phos473);
    Min_vitro_hfd4w_PKB  = vitro_hfd4w.variablevalues(:,id_PKB_phos473);
    
    Max_vitro_chow4w_AS160 = vitro_chow4w.statevalues(:,id_As_phos);
    Min_vitro_chow4w_AS160 = vitro_chow4w.statevalues(:,id_As_phos);
    Max_vitro_hfd4w_AS160  = vitro_hfd4w.statevalues(:,id_As_phos);
    Min_vitro_hfd4w_AS160  = vitro_hfd4w.statevalues(:,id_As_phos);
    
    Max_vitro_chow4w_GU = vitro_chow4w.variablevalues(:,id_GU);
    Min_vitro_chow4w_GU = vitro_chow4w.variablevalues(:,id_GU);
    Max_vitro_hfd4w_GU  = vitro_hfd4w.variablevalues(:,id_GU);
    Min_vitro_hfd4w_GU  = vitro_hfd4w.variablevalues(:,id_GU);
            
end
%%
[ Min_vitro_hfd4w_fold_diet_PKB ,  Max_vitro_hfd4w_fold_diet_PKB ]     = f_calc_minmax_fold( Min_vitro_hfd4w_fold_diet_PKB,   Max_vitro_hfd4w_fold_diet_PKB  , fold_diet_PKB);

[ Min_vitro_hfd4w_fold_diet_AS160 ,  Max_vitro_hfd4w_fold_diet_AS160 ]     = f_calc_minmax_fold( Min_vitro_hfd4w_fold_diet_AS160,   Max_vitro_hfd4w_fold_diet_AS160  , fold_diet_AS160);

[ Min_vitro_hfd4w_fold_diet_GU ,  Max_vitro_hfd4w_fold_diet_GU ]     = f_calc_minmax_fold( Min_vitro_hfd4w_fold_diet_GU,   Max_vitro_hfd4w_fold_diet_GU  , fold_AUC_hfd4w_GU);

%%
[ Min_vitro_chow4w_PKB ,  Max_vitro_chow4w_PKB ]     = f_calc_minmax_fold( Min_vitro_chow4w_PKB,   Max_vitro_chow4w_PKB  , vitro_chow4w.variablevalues(:,id_PKB_phos473));
[ Min_vitro_hfd4w_PKB ,  Max_vitro_hfd4w_PKB ]     = f_calc_minmax_fold( Min_vitro_hfd4w_PKB,   Max_vitro_hfd4w_PKB  , vitro_hfd4w.variablevalues(:,id_PKB_phos473));

[ Min_vitro_chow4w_AS160 ,  Max_vitro_chow4w_AS160 ]     = f_calc_minmax_fold( Min_vitro_chow4w_AS160,   Max_vitro_chow4w_AS160  , vitro_chow4w.statevalues(:,id_As_phos));
[ Min_vitro_hfd4w_AS160 ,  Max_vitro_hfd4w_AS160 ]     = f_calc_minmax_fold( Min_vitro_hfd4w_AS160,   Max_vitro_hfd4w_AS160  , vitro_hfd4w.statevalues(:,id_As_phos));

[ Min_vitro_chow4w_GU ,  Max_vitro_chow4w_GU ]     = f_calc_minmax_fold( Min_vitro_chow4w_GU,   Max_vitro_chow4w_GU  , vitro_chow4w.variablevalues(:,id_GU));
[ Min_vitro_hfd4w_GU ,  Max_vitro_hfd4w_GU ]     = f_calc_minmax_fold( Min_vitro_hfd4w_GU,   Max_vitro_hfd4w_GU  , vitro_hfd4w.variablevalues(:,id_GU));


end
%%

%% Plotting

label_size       = 16 ;
label_size_ax    = 16 ;
marker_line_size = 1.5;
marker_size      = 5  ;
line_size        = 1.5;
title_size       = 18 ;

simulation_TrainColor = [ 0.5,0.5,0.5];
training_dataColor   =  [ 0,  0,  0 ]; 

simulation_predColor =  [0,.9,  .5]; 
validation_dataColor =  [0,.7,  .3 ]; 


%% Simulations for best cost.
load('./Parameters/theta_best.mat'); theta2 = Results.xbest; 

[ AS_train ] = f_as_new( theta2 ,func_mex_model, D, utility);
[ AS3_train ] = f_as_invitro_new_2( theta ,func_mex_model, utility);


% %% TRAINING DATA
PKB_fold_chow = AS3_train.vitro_chow4w.variablevalues(end,id_PKB_phos473)./AS3_train.vitro_chow4w.variablevalues(end,id_PKB);
PKB_fold_hfd  = AS3_train.vitro_hfd4w.variablevalues(end,id_PKB_phos473)./AS3_train.vitro_hfd4w.variablevalues(end,id_PKB);
fold_diet_PKB = PKB_fold_hfd./PKB_fold_chow;% fold between diets ( fold on fold )


AS160_fold_chow = AS3_train.vitro_chow4w.statevalues(end,id_As_phos)./AS3_train.vitro_chow4w.variablevalues(end,id_AS);
AS160_fold_hfd  = AS3_train.vitro_hfd4w.statevalues(end,id_As_phos)./AS3_train.vitro_hfd4w.variablevalues(end,id_AS);
fold_diet_AS160 = AS160_fold_hfd./AS160_fold_chow;% fold between diets ( fold on fold )

AUC_chow4w_GU     = trapz(AS3_train.vitro_chow4w.variablevalues(:,id_GU));
AUC_hfd4w_GU      = trapz(AS3_train.vitro_hfd4w.variablevalues( :,id_GU));
fold_AUC_hfd4w_GU = AUC_hfd4w_GU./AUC_chow4w_GU; % fold on chow


%% Figure 5 - Main text

figure (51)
set(figure(51), 'outerposition',[0 0 1500 500], 'PaperType','a4')

%%% Figure 5 E(i) - PKB
subplot(1,3,1) 
hold on
y = D.Hansson2019.invitro.PKBS473.Mean(:,2);     %[1 D.Hansson2019.invitro.PKBS473.Mean(:,2) 1 fold_diet_PKB ];
b = bar(y, 'grouped','FaceColor','flat' ,'HandleVisibility','off');
b(1).CData(1,:) = [1 1 1]; 
errorbar(1 ,D.Hansson2019.invitro.PKBS473.Mean(:,2) ,D.Hansson2019.invitro.PKBS473.SEM(:,2),'k*' ,'linewidth',1.5)
xlim([0 3]); xticks([ 1  2]);
errhigh = abs(Max_vitro_hfd4w_fold_diet_PKB(end) - fold_diet_PKB   );
errlow  = abs(Min_vitro_hfd4w_fold_diet_PKB(end) - fold_diet_PKB   );
rectangle('Position',[ 1.625  Min_vitro_hfd4w_fold_diet_PKB(end) 0.75  (errhigh + errlow) ],'FaceColor',simulation_TrainColor,'EdgeColor','none',...
    'LineWidth',0.5)
plot( [1.625 1.625+ 0.75 ], [ 1  1 ].*fold_diet_PKB,  'color', simulation_TrainColor - 0.2 ,'linewidth',1.5)
ylabel({'PKB pS473 / total PKB';'(Fold of avg. control)'});


%%% Figure 5 E(ii) - AS160
subplot(1,3,2) 
hold on
y = D.Hansson2019.invitro.AS160.Mean(:,2);     %[1 D.Hansson2019.invitro.PKBS473.Mean(:,2) 1 fold_diet_PKB ];
b = bar(y, 'grouped','FaceColor','flat' ,'HandleVisibility','off');
b(1).CData(1,:) = [1 1 1]; 
errorbar(1 ,D.Hansson2019.invitro.AS160.Mean(:,2) ,D.Hansson2019.invitro.AS160.SEM(:,2),'k*' ,'linewidth',1.5)
xlim([0 3]); xticks([ 1  2]);
errhigh = abs(Max_vitro_hfd4w_fold_diet_AS160(end) - fold_diet_AS160   );
errlow  = abs(Min_vitro_hfd4w_fold_diet_AS160(end) - fold_diet_AS160   );
rectangle('Position',[ 1.625  Min_vitro_hfd4w_fold_diet_AS160(end) 0.75  (errhigh + errlow) ],'FaceColor',simulation_TrainColor,'EdgeColor','none',...
    'LineWidth',0.5)
plot( [1.625 1.625+ 0.75 ], [ 1  1 ].*fold_diet_AS160,  'color', simulation_TrainColor - 0.2 ,'linewidth',1.5)
ylabel({'AS160 pT642 / total AS160';'(Fold of avg. control)'});

%%% Figure 5 E(iii) - GLUT4 - Glucose uptake
subplot(1,3,3) 
hold on
y = D.Hansson2019.invitro.GU.Mean(:,2);     %[1 D.Hansson2019.invitro.PKBS473.Mean(:,2) 1 fold_diet_PKB ];
b = bar(y, 'grouped','FaceColor','flat' ,'HandleVisibility','off');
b(1).CData(1,:) = [1 1 1]; 
errorbar(1 ,D.Hansson2019.invitro.GU.Mean(:,2) ,D.Hansson2019.invitro.GU.SEM(:,2),'k*' ,'linewidth',1.5)
xlim([0 3]); xticks([ 1  2]);
errhigh = abs(Max_vitro_hfd4w_fold_diet_GU(end) - fold_AUC_hfd4w_GU   );
errlow  = abs(Min_vitro_hfd4w_fold_diet_GU(end) - fold_AUC_hfd4w_GU   );
rectangle('Position',[ 1.625  Min_vitro_hfd4w_fold_diet_GU(end)   0.75  (errhigh + errlow) ],'FaceColor',simulation_TrainColor,'EdgeColor','none',...
    'LineWidth',0.5)
plot( [1.625 1.625+ 0.75 ], [ 1  1 ].*fold_AUC_hfd4w_GU,  'color', simulation_TrainColor - 0.2 ,'linewidth',1.5 )
ylabel({'AUC adipocse glucose uptake ';'(Fold of avg. control'});
 

%% Figure 
figure (52)
set(figure(52), 'outerposition',[0 0 1200 1200], 'PaperType','a4')

t_IPGTT=(0:1:120); 
fill_time_IPGTT = [ t_IPGTT fliplr(t_IPGTT) ];

sim_id = {'IPGTT_chow2w','IPGTT_chow12dHFD2d','IPGTT_chow10dHFD4d','IPGTT_chow8dHFD6d','IPGTT_hfd2w','IPGTT_hfd6w'};
state_name = {'PKB473' ,'PKB308473','AS160p','AS160','RABGTP','PKBda','PMAa'};
clear tmp2; clear tmp3;
for ii = 1:length(state_name)
    for jj = 1:length(sim_id)
        values1=[];
        eval(strcat('values1(:,1)=Max_',sim_id{jj}, '_', state_name{ii} ,';'));
        eval(strcat('values1(:,2)=Min_',sim_id{jj}, '_', state_name{ii} ,';'));
        if  jj == 1
        tmp2(ii) = max(values1(:,1)) ;
        tmp3(ii) = min(values1(:,2)) ;
        end
        if tmp2(ii) < max(values1(:,1)) ; tmp2(ii) = max(values1(:,1)); end
        if tmp3(ii) > min(values1(:,2)) ; tmp3(ii) = min(values1(:,2)); end
        
    end
end; clear ii; clear jj;

for ii = 1:length(state_name)
    for jj = 1:length(sim_id)
        
        subplot( length(state_name) + 1, length(sim_id), length(sim_id)*(ii-1) + jj  )
        hold on
        
        
        time_fill = fill_time_IPGTT;
        
        values1=[];
        eval(strcat('values1(:,1)=Max_',sim_id{jj}, '_', state_name{ii} ,';'));
        eval(strcat('values1(:,2)=Min_',sim_id{jj}, '_', state_name{ii} ,';'));
        f=fill(time_fill, [ values1(:,2)' fliplr( values1(:,1)' ) ],simulation_predColor, 'EdgeColor', 'none' , 'HandleVisibility','off');
        set(f,'facealpha',.25)
        hold on
        eval(strcat('plot( t_IPGTT, AS_train.',sim_id{jj}, ".statevalues(:,id_", state_name{ii} ,"),'color', simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', 2);"))
        
        box off
        xlim([0 120]);
        set(gca,'FontSize',(label_size_ax-6))
        if ii == length(state_name)
            xlabel('Time (min)','FontSize', label_size-6 )
        end
        
        if jj == 1
            ylabel( strcat( state_name{ii}, ' (a.u.)'),'FontSize', label_size-6 )           
        end
        ylim( [  (tmp3(ii) - tmp3(ii)*0.05)   (tmp2(ii) + tmp2(ii)*0.05) ])
        
    end
end

sim_id = {'IPGTT_chow2w','IPGTT_chow12dHFD2d','IPGTT_chow10dHFD4d','IPGTT_chow8dHFD6d','IPGTT_hfd2w','IPGTT_hfd6w'};
varb_name = { 'GLUT4m'}; % 'insulin_input',

for ii = 1:length(varb_name)
    for jj = 1:length(sim_id)
        subplot( length(state_name) + 1, length(sim_id), length(sim_id)*(length(state_name)) + jj  )
        hold on
        values1=[];
        time_fill = fill_time_IPGTT;
        
        eval(strcat('values1(:,1)=Max_',sim_id{jj}, '_', varb_name{ii} ,';'));
        eval(strcat('values1(:,2)=Min_',sim_id{jj}, '_', varb_name{ii} ,';' ));
        f=fill(time_fill, [ values1(:,2)' fliplr( values1(:,1)' ) ], simulation_predColor, 'EdgeColor', 'none' , 'HandleVisibility','off');
        set(f,'facealpha',.25)
        hold on
        eval(strcat('plot( t_IPGTT, AS_train.',sim_id{jj}, ".variablevalues(:,id_", varb_name{ii} ,"),'color', simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', 2);"))
        
        box off
        set(gca,'FontSize',label_size_ax-6)
        xlabel('Time (min)','FontSize', label_size -6)
        xlim([0 120]);
        %if ii == 1 ylim([0 1.5]); elseif ii == 2 ylim([0 0.5]); end
        
        if jj == 1 && ii == 1
%             ylabel( 'Insulin (nM)' ,'FontSize', label_size-6 )
%         elseif jj == 1 && ii == 2
            ylabel( 'Membrane GLUT4 (a.u.)','FontSize', label_size-6 )
            ylim([0 0.5])
        end
        
    end
end



%%
function[min ,max ] = f_calc_minmax_fold( min, max , fold)
sim_string = 'fold';
eval(strcat('max( max < ', sim_string,')=', sim_string,'( max <', sim_string ,');'))
eval(strcat('min( min > ', sim_string,')= ',sim_string,'( min >', sim_string ,');'))

end
