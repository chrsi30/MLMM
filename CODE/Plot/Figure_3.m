
%%

id_BW=ismember(SBvariables(objModel),'BW');
id_Glu = ismember(SBvariables(objModel),'Gconc');
id_Ins = ismember(SBvariables(objModel),'Iconc');

%% UNC
for i = 1:nfiles
load( strcat( filedir ,'/',files(i).name ) )
disp(strcat('Processing - ',files(i).name, int2str(i),'/',int2str(nfiles)))
theta=Results.xbest;

[ AS ] = f_as_new( theta ,func_mex_model, D, utility);

ss                   = AS.ss;
chow8w               = AS.chow8w;
hfd8w                = AS.hfd8w;
IPGTT_hfd6w          = AS.IPGTT_hfd6w;
IPGTT_chow6w          = AS.IPGTT_chow6w;
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

%% MIN and MAX value calculations

%% Predictions
if i == 1 % Setting up vectors
Max_chow8w_BW  =  chow8w.variablevalues(:,id_BW);  Min_chow8w_BW =  chow8w.variablevalues(:,id_BW);     
Max_hfd8w_BW   =  hfd8w.variablevalues(:,id_BW);  Min_hfd8w_BW =  hfd8w.variablevalues(:,id_BW);      
Max_hfd6w_Glu  =  IPGTT_hfd6w.variablevalues(:,id_Glu);  Min_hfd6w_Glu =  IPGTT_hfd6w.variablevalues(:,id_Glu); 
Max_chow6w_Glu =  IPGTT_chow6w.variablevalues(:,id_Glu);  Min_chow6w_Glu =  IPGTT_chow6w.variablevalues(:,id_Glu); 
Max_hfd6w_Ins  =  IPGTT_hfd6w.variablevalues(:,id_Ins);  Min_hfd6w_Ins =  IPGTT_hfd6w.variablevalues(:,id_Ins); 
Max_chow6w_Ins =  IPGTT_chow6w.variablevalues(:,id_Ins);  Min_chow6w_Ins =  IPGTT_chow6w.variablevalues(:,id_Ins); 
end

[ Min_chow8w_BW ,  Max_chow8w_BW ]     = f_calc_minmax( Min_chow8w_BW,   Max_chow8w_BW  , 'variable' , chow8w , id_BW);
[ Min_hfd8w_BW ,  Max_hfd8w_BW ]       = f_calc_minmax( Min_hfd8w_BW,    Max_hfd8w_BW   , 'variable' , hfd8w  , id_BW);

[ Min_hfd6w_Glu ,  Max_hfd6w_Glu ]     = f_calc_minmax( Min_hfd6w_Glu,   Max_hfd6w_Glu  , 'variable' , IPGTT_hfd6w  , id_Glu);
[ Min_chow6w_Glu ,  Max_chow6w_Glu ]   = f_calc_minmax( Min_chow6w_Glu,  Max_chow6w_Glu , 'variable' , IPGTT_chow6w  , id_Glu);

[ Min_chow6w_Ins ,  Max_chow6w_Ins]   = f_calc_minmax( Min_chow6w_Ins,  Max_chow6w_Ins , 'variable' , IPGTT_chow6w  , id_Ins);
[ Min_hfd6w_Ins ,  Max_hfd6w_Ins]     = f_calc_minmax( Min_hfd6w_Ins,  Max_hfd6w_Ins , 'variable' , IPGTT_hfd6w  , id_Ins);


%% 
id_tp_bw_8w_2 = [ 22 29 36 42 51 57];
bw_cost8 = nansum (  ( (chow8w.variablevalues(id_tp_bw_8w_2, id_BW)    -  D.NEW.BW.Mean(4:end,1) ).^2 )./( D.NEW.BW.SEM(4:end,1).^2 ));
bw_cost9 = nansum (  ( (hfd8w.variablevalues(id_tp_bw_8w_2,  id_BW)    -  D.NEW.BW.Mean(4:end,2) ).^2 )./( D.NEW.BW.SEM(4:end,2).^2 ));
bw_cost_up_pred = bw_cost8 + bw_cost9;

id_tp_IPGTT = [ 1 16 31 61 121 ];
Glu_prediction_cost = nansum (  ( (IPGTT_hfd6w.variablevalues(id_tp_IPGTT, id_Glu)    -  D.NEW.IPGTT.Mean(:,2) ).^2 )./( D.NEW.IPGTT.SEM(:,2).^2 ));
IPGTT_glucose_cost_6w =  Glu_prediction_cost;

id_tp_bw_8w = [ 1 8 15];
bw_cost6 = nansum (  ( (chow8w.variablevalues(id_tp_bw_8w, id_BW)    -  D.NEW.BW.Mean(1:3,1) ).^2 )./( D.NEW.BW.SEM(1:3,1).^2 ));
bw_cost7 = nansum (  ( (hfd8w.variablevalues(id_tp_bw_8w,  id_BW)    -  D.NEW.BW.Mean(1:3,2) ).^2 )./( D.NEW.BW.SEM(1:3,2).^2 )); 

bw_cost_up = bw_cost6 + bw_cost7;

%% Hansson et al (2018)
if i == 1  % Setting up vectors
% 2w chow    
Max_chow2w_BW =  chow2w.variablevalues(:,id_BW);          Min_chow2w_BW =  chow2w.variablevalues(:,id_BW);     
Max_chow2w_Glu =  IPGTT_chow2w.variablevalues(:,id_Glu);  Min_chow2w_Glu =  IPGTT_chow2w.variablevalues(:,id_Glu); 
Max_chow2w_Ins = IPGTT_chow2w.variablevalues(:, id_Ins);  Min_chow2w_Ins =  IPGTT_chow2w.variablevalues(:, id_Ins); 

%12d chow + 2d HFD  
Max_chow12dHFD2d_BW =  chow12dHFD2d.variablevalues(:,id_BW);          Min_chow12dHFD2d_BW =  chow12dHFD2d.variablevalues(:,id_BW);     
Max_chow12dHFD2d_Glu =  IPGTT_chow12dHFD2d.variablevalues(:,id_Glu);  Min_chow12dHFD2d_Glu =  IPGTT_chow12dHFD2d.variablevalues(:,id_Glu); 
Max_chow12dHFD2d_Ins = IPGTT_chow12dHFD2d.variablevalues(:, id_Ins);  Min_chow12dHFD2d_Ins =  IPGTT_chow12dHFD2d.variablevalues(:, id_Ins);

%10d chow + 4d HFD  
Max_chow10dHFD4d_BW =  chow10dHFD4d.variablevalues(:,id_BW);          Min_chow10dHFD4d_BW =  chow10dHFD4d.variablevalues(:,id_BW); 
Max_chow10dHFD4d_Glu =  IPGTT_chow10dHFD4d.variablevalues(:,id_Glu);  Min_chow10dHFD4d_Glu =  IPGTT_chow10dHFD4d.variablevalues(:,id_Glu); 
Max_chow10dHFD4d_Ins = IPGTT_chow10dHFD4d.variablevalues(:, id_Ins);  Min_chow10dHFD4d_Ins =  IPGTT_chow10dHFD4d.variablevalues(:, id_Ins);

%8d chow + 6d HFD  
Max_chow8dHFD6d_BW =  chow8dHFD6d.variablevalues(:,id_BW);          Min_chow8dHFD6d_BW =  chow8dHFD6d.variablevalues(:,id_BW);     
Max_chow8dHFD6d_Glu =  IPGTT_chow8dHFD6d.variablevalues(:,id_Glu);  Min_chow8dHFD6d_Glu =  IPGTT_chow8dHFD6d.variablevalues(:,id_Glu); 
Max_chow8dHFD6d_Ins = IPGTT_chow8dHFD6d.variablevalues(:, id_Ins);  Min_chow8dHFD6d_Ins =  IPGTT_chow8dHFD6d.variablevalues(:, id_Ins); 

% 2w hfd
Max_hfd2w_BW =  hfd2w.variablevalues(:,id_BW);          Min_hfd2w_BW =  hfd2w.variablevalues(:,id_BW);     
Max_hfd2w_Glu =  IPGTT_hfd2w.variablevalues(:,id_Glu);  Min_hfd2w_Glu =  IPGTT_hfd2w.variablevalues(:,id_Glu); 
Max_hfd2w_Ins = IPGTT_hfd2w.variablevalues(:, id_Ins);  Min_hfd2w_Ins =  IPGTT_hfd2w.variablevalues(:, id_Ins); 
end 

[ Min_chow2w_BW ,  Max_chow2w_BW ]     = f_calc_minmax( Min_chow2w_BW,   Max_chow2w_BW  , 'variable' , chow2w , id_BW);
[ Min_chow2w_Glu ,  Max_chow2w_Glu ]   = f_calc_minmax( Min_chow2w_Glu,  Max_chow2w_Glu , 'variable' , IPGTT_chow2w , id_Glu);
[ Min_chow2w_Ins ,  Max_chow2w_Ins ]   = f_calc_minmax( Min_chow2w_Ins,  Max_chow2w_Ins , 'variable' , IPGTT_chow2w , id_Ins);

[ Min_chow12dHFD2d_BW ,  Max_chow12dHFD2d_BW ]     = f_calc_minmax( Min_chow12dHFD2d_BW,   Max_chow12dHFD2d_BW  , 'variable' , chow12dHFD2d , id_BW);
[ Min_chow12dHFD2d_Glu ,  Max_chow12dHFD2d_Glu ]   = f_calc_minmax( Min_chow12dHFD2d_Glu,  Max_chow12dHFD2d_Glu , 'variable' , IPGTT_chow12dHFD2d , id_Glu);
[ Min_chow12dHFD2d_Ins ,  Max_chow12dHFD2d_Ins ]   = f_calc_minmax( Min_chow12dHFD2d_Ins,  Max_chow12dHFD2d_Ins , 'variable' , IPGTT_chow12dHFD2d , id_Ins);

[ Min_chow10dHFD4d_BW ,  Max_chow10dHFD4d_BW ]     = f_calc_minmax( Min_chow10dHFD4d_BW,   Max_chow10dHFD4d_BW  , 'variable' , chow10dHFD4d , id_BW);
[  Min_chow10dHFD4d_Glu ,  Max_chow10dHFD4d_Glu ]   = f_calc_minmax( Min_chow10dHFD4d_Glu,  Max_chow10dHFD4d_Glu , 'variable' , IPGTT_chow10dHFD4d , id_Glu);
[ Min_chow10dHFD4d_Ins  ,  Max_chow10dHFD4d_Ins ]   = f_calc_minmax( Min_chow10dHFD4d_Ins,  Max_chow10dHFD4d_Ins , 'variable' , IPGTT_chow10dHFD4d , id_Ins);

[ Min_chow8dHFD6d_BW ,  Max_chow8dHFD6d_BW ]     = f_calc_minmax( Min_chow8dHFD6d_BW,   Max_chow8dHFD6d_BW  , 'variable' , chow8dHFD6d , id_BW);
[ Min_chow8dHFD6d_Glu ,  Max_chow8dHFD6d_Glu ]   = f_calc_minmax( Min_chow8dHFD6d_Glu,  Max_chow8dHFD6d_Glu , 'variable' , IPGTT_chow8dHFD6d , id_Glu);
[ Min_chow8dHFD6d_Ins ,  Max_chow8dHFD6d_Ins ]   = f_calc_minmax( Min_chow8dHFD6d_Ins,  Max_chow8dHFD6d_Ins , 'variable' , IPGTT_chow8dHFD6d , id_Ins);

[ Min_hfd2w_BW ,   Max_hfd2w_BW ]     = f_calc_minmax( Min_hfd2w_BW,   Max_hfd2w_BW  , 'variable' , hfd2w , id_BW);
[ Min_hfd2w_Glu ,  Max_hfd2w_Glu ]   = f_calc_minmax( Min_hfd2w_Glu,  Max_hfd2w_Glu , 'variable' , IPGTT_hfd2w , id_Glu);
[ Min_hfd2w_Ins ,  Max_hfd2w_Ins ]   = f_calc_minmax( Min_hfd2w_Ins,  Max_hfd2w_Ins , 'variable' , IPGTT_hfd2w , id_Ins);

%% Cost calculations

if i == 1 
    x_pred = bw_cost_up_pred + IPGTT_glucose_cost_6w;
    x_training = cost_function_v6( theta ,func_mex_model, D, utility);
    best_param_train = theta;
    best_param_pred = theta;
end

cost_training = cost_function_v6( theta ,func_mex_model, D, utility);
cost_pred = bw_cost_up_pred + IPGTT_glucose_cost_6w;

if x_training > cost_training
x_training = cost_training;
best_param_train = theta;
end

if x_pred > cost_pred
    x_pred = cost_pred;
    best_bw_cost_up_pred=bw_cost_up_pred;
    best_IPGTT_glucose_cost_6w = IPGTT_glucose_cost_6w;
    best_param_pred = theta;
end
end

%% Simulations for best cost.
load('./Parameters/theta_best.mat'); theta2 = Results.xbest; 

[ AS_train ] = f_as_new( theta2 ,func_mex_model, D, utility);


theta1 = best_param_pred;

[ AS_pred ] = f_as_new( theta1 ,func_mex_model, D, utility);


%% Figure 3 - B

simulation_predColor =  [0,.9,  .5]; 
validation_dataColor =  [0,.7,  .3 ]; 

simulation_TrainColor = [ 0.5,0.5,0.5];
training_dataColor   =  [0, 0,  0 ]; 

label_size       = 16 ;
label_size_ax    = 16 ;
marker_line_size = 1.5;
marker_size      = 5  ;
line_size        = 1.5;
title_size       = 18 ;
%%
figure(31)
set(figure(31), 'outerposition',[0 0 2000 1000], 'PaperType','a4')

%%% time resolution
t_diet=0:1:14;
fill_time_2w = [ t_diet fliplr(t_diet) ];

%%% Figure 3b first row -  BW - model training to Hansson 2018 et al data
Experiments_to_plot = {'chow2w', 'chow12dHFD2d', 'chow10dHFD4d','chow8dHFD6d','hfd2w'};
for kk = 1:length(Experiments_to_plot)
    
    
    values=[];
    eval(strcat('values(:,1)=Max_',Experiments_to_plot{kk},'_BW;'));
    eval(strcat('values(:,2)=Min_',Experiments_to_plot{kk},'_BW;'));
    
    
    
    if kk == 1 % "
        time_plot = t_diet;
        time_fill = fill_time_2w;
        tmp1 = 1;
        tstring = '2 weeks control' ;
        
    elseif kk == 2
        time_plot =12:1:t_diet(end);
        time_fill = [ time_plot fliplr(time_plot) ];
        tstring = '2 days HFD' ;
        
    elseif kk == 3
        time_plot = 10:1:t_diet(end);
        time_fill = [ time_plot fliplr(time_plot) ];
        tmp1 = 2;
        tstring = '4 days HFD' ;
        
    elseif kk == 4
        time_plot =8:1:t_diet(end);
        time_fill = [ time_plot fliplr(time_plot) ];
        tmp1 = 3;
        tstring = '6 days HFD' ;
        
    elseif kk == 5
        time_plot = t_diet;
        time_fill = fill_time_2w;
        tmp1 = 4;
        
        tstring = '2 weeks HFD' ;
        
    end
    
    
    subplot(3,5,kk)
    
    title(tstring,'FontSize', title_size )
    
    hold on
    
    
    if kk == 2
        f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], simulation_predColor, 'EdgeColor', 'none' );
        set(f,'facealpha',.5)
        hold on
        
        eval(strcat('plot( time_plot, AS_train.',Experiments_to_plot{kk}, ".variablevalues(:,id_BW),'color', simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', line_size);"))
        
        values=[];
        values(:,1)=Max_chow2w_BW;
        values(:,2)=Min_chow2w_BW;
        time_plot =0:1:12;
        time_fill = [ time_plot fliplr(time_plot) ];
        f=fill(time_fill, [ values(1:13,2)' fliplr( values(1:13,1)' ) ], simulation_predColor, 'EdgeColor', 'none' );
        set(f,'facealpha',.5)
        hold on
        plot( time_plot, AS_train.chow2w.variablevalues(1:13,id_BW),'color', simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', line_size);
        
    else
        f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], simulation_TrainColor, 'EdgeColor', 'none' );
        set(f,'facealpha',.5)
        hold on
        
        eval(strcat('plot( time_plot, AS_train.',Experiments_to_plot{kk}, ".variablevalues(:,id_BW),'color', simulation_TrainColor - 0.2 , 'linewidth', line_size);"))
        
        data.time = D.Hansson2018.BW.time(:,tmp1);
        data.mean = D.Hansson2018.BW.Mean(:,tmp1);
        data.SEM  = D.Hansson2018.BW.SEM(:,tmp1);
        errorbar(data.time,data.mean, data.SEM, '*','color', training_dataColor ,...
        'MarkerSize',marker_size,'linewidth',marker_line_size ) 
        
    end

 
    set(gca,'FontSize',label_size_ax)
    xlabel('Time (days)','FontSize', label_size )
    ylabel('Body weight (g)','FontSize', label_size )
    box off
    xlim([0 14]);
    xticks([0 ,14])
    xticklabels({'0','14'})
    ylim([20 34]);
    
end


%%% Figure 3b second row -  Gluocse IPGTT - model training to Hansson 2018 et al data
t_diet=14;
t_fasting=[ t_diet  t_diet+0.5833 ]  ;  %t_diet:0.0001:t_diet+0.5833;
t_IPGTT=(0:1:120)/(60*24)+t_fasting(end);

fill_time_IPGTT = [ t_IPGTT-t_fasting(end) fliplr(t_IPGTT-t_fasting(end)) ]*24*60;
tp_IPGTT=[0 15 30 60 120];
time_best_cost = (t_IPGTT-t_fasting(end) )*24*60; 

for kk = 1:length(Experiments_to_plot)
    
    
    values=[];
    eval(strcat('values(:,1)=Max_',Experiments_to_plot{kk},'_Glu;'));
    eval(strcat('values(:,2)=Min_',Experiments_to_plot{kk},'_Glu;'));
    
    
    if kk == 1      % 2w control
        tmp1 = 1;
    elseif kk == 2  % 2d hfd
        tmp1 = 2;
    elseif kk == 3  % 4d hfd
        
    elseif kk == 4  % 6 hfd
        tmp1 = 3;
    elseif kk == 5  % 2w hfd
        tmp1 = 4;
    end
    
    
    subplot(3,5,5 + kk)

    time_fill = fill_time_IPGTT;
    time_plot = (t_IPGTT-t_fasting(end))*24*60;
    
    if kk == 3
    f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], simulation_predColor, 'EdgeColor', 'none' );
    else
    f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], simulation_TrainColor, 'EdgeColor', 'none' );
    end
    
    
    set(f,'facealpha',.5)
    hold on
    
    if kk == 3
    else
        data.time = D.Hansson2018.IPGTT.Glucose.time;
        data.mean =  D.Hansson2018.IPGTT.Glucose.Mean(:,tmp1) ;
        data.SEM  =   D.Hansson2018.IPGTT.Glucose.SEM(:,tmp1);
        errorbar(data.time,data.mean, data.SEM, '*','color', training_dataColor ,...
        'MarkerSize',marker_size,'linewidth',marker_line_size ) 
    end
    
    if kk == 3
        eval(strcat('plot( time_plot, AS_train.IPGTT_',Experiments_to_plot{kk}, ".variablevalues(:,id_Glu),'color',     simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', line_size);"))  
    else
        eval(strcat('plot( time_plot, AS_train.IPGTT_',Experiments_to_plot{kk}, ".variablevalues(:,id_Glu),'color', simulation_TrainColor - 0.2 , 'linewidth', line_size);"))
    end
        
    box off
    set(gca,'FontSize',label_size_ax)
    ylabel(' Plasma glucose (mM)','FontSize', label_size)
    xlabel('Time (min)','FontSize', label_size )
    xlim([0 120]);
    ylim([4 32]);
    
   if kk == 1
       legend({'1','2','3'})
   end 
   
end

%%% Figure 3b third row -  Insulin IPGTT - model training to Hansson 2018 et al data

for kk = 1:length(Experiments_to_plot)
    
    if kk == 1      % 2w control
        tmp1 = 1;
    elseif kk == 2  % 2d hfd
        tmp1 = 2;
    elseif kk == 3  % 4d hfd
        
    elseif kk == 4  % 6 hfd
        tmp1 = 3;
    elseif kk == 5  % 2w hfd
        tmp1 = 4;
    end
    
    subplot(3,5,10 + kk)
    hold on
    
    time_fill = fill_time_IPGTT;
    time_plot = (t_IPGTT-t_fasting(end))*24*60;
    
    values=[];
    eval(strcat('values(:,1)=Max_',Experiments_to_plot{kk},'_Ins;'));
    eval(strcat('values(:,2)=Min_',Experiments_to_plot{kk},'_Ins;'));
    
    f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], simulation_predColor, 'EdgeColor', 'none' );
    set(f,'facealpha',.5)
    hold on
    
    eval(strcat('plot( time_plot, AS_train.IPGTT_',Experiments_to_plot{kk}, ".variablevalues(:,id_Ins),'color', simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', line_size);"))
    
    if kk == 3
    else
        data.time =  D.Hansson2018.IPGTT.Insulin.time;
        data.mean =  D.Hansson2018.IPGTT.Insulin.Mean(:,tmp1);
        data.SEM  =  D.Hansson2018.IPGTT.Insulin.SEM(:,tmp1);
        errorbar(data.time,data.mean, data.SEM, '^','color', validation_dataColor ,...
            'MarkerSize',marker_size,'linewidth',marker_line_size )
    end
    
    
    box off
    set(gca,'FontSize',label_size_ax)
    ylabel(' Plasma Insulin (mM)','FontSize', label_size)
    xlabel('Time (min)','FontSize', label_size )
    xlim([0 120]);
    ylim([0 900]); yticks([ 0 400 800])
    if kk == 1
        legend({'1','2','3'})
    end
end

%%
figure (32)
set(figure(32), 'outerposition',[0 0 2000 1000], 'PaperType','a4')

Experiments_to_plot = {'chow8w','hfd8w'};
for kk = 1:length(Experiments_to_plot)
    
    values=[];
    eval(strcat('values(:,1)=Max_',Experiments_to_plot{kk},'_BW;'));
    eval(strcat('values(:,2)=Min_',Experiments_to_plot{kk},'_BW;'));
    
    if kk == 1
        data.time=D.NEW.BW.time;
        data.mean= D.NEW.BW.Mean(:,1);
        data.SEM=D.NEW.BW.SEM(:,1) ;
        tstring = '8 weeks control';
    elseif kk == 2
        data.time=D.NEW.BW.time;
        data.mean= D.NEW.BW.Mean(:,2);
        data.SEM=D.NEW.BW.SEM(:,2) ;
        tstring = '8 weeks HFD';

    end
    
    subplot(3,5,kk)
    title(tstring,'FontSize', title_size )
    hold on

    time_fill = fill_time_2w;
    time_plot = 0:1:56;

    f=fill(time_fill, [ values(1:15,2)' fliplr( values(1:15,1)' ) ], simulation_TrainColor, 'EdgeColor', 'none' );
    set(f,'facealpha',.5)
    hold on
    fill_time_8w = [ 14:1:56 fliplr(14:1:56) ];
    f=fill(fill_time_8w, [ values(15:end,2)' fliplr( values(15:end,1)' ) ], simulation_predColor, 'EdgeColor', 'none' );
    set(f,'facealpha',.5)
    
    % ploting best sim for training and pred
    eval(strcat('plot( time_plot, AS_train.',Experiments_to_plot{kk}, ".variablevalues(:,id_BW),'color', simulation_TrainColor - 0.2 , 'linewidth', line_size);"))
    eval(strcat('plot( time_plot, AS_pred.',Experiments_to_plot{kk}, ".variablevalues(:,id_BW),'color', simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', line_size);"))
    
    
    errorbar(data.time(1:3),data.mean(1:3), data.SEM(1:3), '*','color', training_dataColor ,...
        'MarkerSize',marker_size,'linewidth',marker_line_size )
    
    errorbar(data.time(4:end),data.mean(4:end), data.SEM(4:end), '^','color', validation_dataColor ,...
        'MarkerSize',marker_size,'linewidth',marker_line_size )

    set(gca,'FontSize',label_size_ax)
    xlabel('Time (days)','FontSize', label_size )
    ylabel('Body weight (g)','FontSize', label_size )
    box off
    xlim([0 56]);
    xticks([0 ,56])
    xticklabels({'0','56'})
    ylim([20 60]);
    if kk == 1
        legend({'1','2','3'})
    end
    
end

%%% 
time_fill = fill_time_IPGTT;
time_plot = (t_IPGTT-t_fasting(end))*24*60;


tstring = 'IPGTT - 6 weeks control';
subplot(3,5,6)
title(tstring,'FontSize', title_size )
hold on

hold on
data.time=D.NEW.IPGTT.time;
data.mean=D.NEW.IPGTT.Mean(:,1);
data.SEM =D.NEW.IPGTT.SEM(:,1);

values=[];
values(:,1)=Max_chow6w_Glu; values(:,2)=Min_chow6w_Glu; 

f=fill(fill_time_IPGTT, [ values(:,2)' fliplr( values(:,1)' ) ], simulation_predColor ,'EdgeColor', 'none' );
set(f,'facealpha',.5)
hold on

plot( time_plot, AS_train.IPGTT_chow6w.variablevalues(:,id_Glu),'color', simulation_TrainColor - 0.2 , 'linewidth', 2);
plot( time_plot, AS_pred.IPGTT_chow6w.variablevalues(:,id_Glu),'color',simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', 2);

errorbar(data.time,data.mean, data.SEM, '^','color', validation_dataColor ,...
'MarkerSize',marker_size,'linewidth',marker_line_size )

set(gca,'FontSize',label_size_ax)
ylabel(' Plasma glucose (mM)','FontSize', label_size ,'FontWeight','bold')
xlabel('Time (min)','FontSize', label_size )
xlim([0 120]);
ylim([5 45]);

%%%
subplot(3,5,7)
tstring = 'IPGTT - 6 weeks HFD';
title(tstring,'FontSize', title_size )
hold on
    
hold on
data.time=D.NEW.IPGTT.time;
data.mean=D.NEW.IPGTT.Mean(:,2);
data.SEM =D.NEW.IPGTT.SEM(:,2);

values=[];
values(:,1)=Max_hfd6w_Glu; values(:,2)=Min_hfd6w_Glu; 

f=fill(fill_time_IPGTT, [ values(:,2)' fliplr( values(:,1)' ) ], simulation_predColor ,'EdgeColor', 'none' );
set(f,'facealpha',.5)
hold on

plot( time_plot, AS_train.IPGTT_hfd6w.variablevalues(:,id_Glu),'color', simulation_TrainColor - 0.2 , 'linewidth', line_size);
plot( time_plot, AS_pred.IPGTT_hfd6w.variablevalues(:,id_Glu),'color',simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', line_size);

errorbar(data.time,data.mean, data.SEM, '^','color', validation_dataColor ,...
'MarkerSize',marker_size,'linewidth',marker_line_size )
    
set(gca,'FontSize',label_size_ax)
ylabel(' Plasma glucose (mM)','FontSize', label_size ,'FontWeight','bold')
xlabel('Time (min)','FontSize', label_size )
    xlim([0 120]);
    ylim([5 45]);

%%%

subplot(3,5,11)
title(tstring,'FontSize', title_size )
hold on

values=[];
values(:,1)=Max_chow6w_Ins; values(:,2)=Min_chow6w_Ins; 

f=fill(fill_time_IPGTT, [ values(:,2)' fliplr( values(:,1)' ) ], simulation_predColor ,'EdgeColor', 'none' );
set(f,'facealpha',.5)
hold on

plot( time_plot, AS_train.IPGTT_chow6w.variablevalues(:,id_Ins),'color', simulation_TrainColor - 0.2 , 'linewidth', 2);
plot( time_plot, AS_pred.IPGTT_chow6w.variablevalues(:,id_Ins),'color',simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', 2);


set(gca,'FontSize',label_size_ax)
ylabel(' Plasma insulin (pM)','FontSize', label_size ,'FontWeight','bold')
xlabel('Time (min)','FontSize', label_size )
xlim([0 120]);
ylim([0 1500]); yticks([ 0 400 800 1200])

%%%
subplot(3,5,12)
tstring = 'IPGTT - 6 weeks control';
title(tstring,'FontSize', title_size )
hold on
    
values=[];
values(:,1)=Max_hfd6w_Ins; values(:,2)=Min_hfd6w_Ins; 

f=fill(fill_time_IPGTT, [ values(:,2)' fliplr( values(:,1)' ) ], simulation_predColor ,'EdgeColor', 'none' );
set(f,'facealpha',.5)
hold on

plot( time_plot, AS_train.IPGTT_hfd6w.variablevalues(:,id_Ins),'color', simulation_TrainColor - 0.2 , 'linewidth', line_size);
plot( time_plot, AS_pred.IPGTT_hfd6w.variablevalues(:,id_Ins),'color',simulation_predColor - [ 0 , 0.2, 0.2 ] , 'linewidth', line_size);

set(gca,'FontSize',label_size_ax)
ylabel(' Plasma insulin (pM)','FontSize', label_size ,'FontWeight','bold')
xlabel('Time (min)','FontSize', label_size )
xlim([0 120]);
ylim([0 1500]); yticks([ 0 400 800 1200])





