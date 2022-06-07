
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

time_points_IPGTTs = 56; % number of time point to calc unc for

[AS2] = f_as_IPGTTs( Results.xbest ,func_mex_model,D, utility, time_points_IPGTTs) ;

for kk = 1:length(time_points_IPGTTs)  
    tmp = strcat('IPGTT_chow', int2str(time_points_IPGTTs(kk)),'d  = AS2.IPGTT_chow', int2str(time_points_IPGTTs(kk)),'d ;');
    eval(tmp);
    name_1(kk,:) = strcat('IPGTT_chow', int2str(time_points_IPGTTs(kk)),'d' );
    
    tmp = strcat('IPGTT_hfd', int2str(time_points_IPGTTs(kk)),'d  = AS2.IPGTT_hfd', int2str(time_points_IPGTTs(kk)),'d ;');
    eval(tmp);
    name_2(kk,:) = strcat('IPGTT_hfd', int2str(time_points_IPGTTs(kk)),'d' );
    
end


%% MIN and MAX value calculations
sim_id = {'IPGTT_chow2w','IPGTT_chow12dHFD2d','IPGTT_chow10dHFD4d','IPGTT_chow8dHFD6d','IPGTT_hfd2w','IPGTT_chow6w','IPGTT_hfd6w'};
varb_name = {'U_id','U_ii','U_id_A','U_id_M','U_id_H', 'PlotEGP', 'PlotGdiff' , 'PlotGIP'   ,'PlotU'       };

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

sim_id2 = string(name_2)';
varb_name2 = {'U_id_A'} ; %%,'U_ii','U_id_A','U_id_M','U_id_H'};

for kk = 1:length(sim_id2)
    for jj = 1:length(varb_name2)
        tmp1 = strcat('Max_', sim_id2{kk}, '_', varb_name2{jj} );
        tmp2 = strcat('Min_', sim_id2{kk}, '_', varb_name2{jj} );
        eval(strcat('id_',  varb_name2{jj},  "= ismember(SBvariables(objModel),'", varb_name2{jj},"');"));
        eval(strcat('tmp3 = id_',  varb_name2{jj},';') );
        
        if i == 1
            tmp4 = strcat( sim_id2{kk} , '.variablevalues(:, tmp3 );');
            eval(strcat( tmp1 ,'=' ,tmp4)); eval(strcat( tmp2 ,'=' ,tmp4));
        end
        eval(strcat('[ ',tmp2,' ,  ',tmp1,' ]     = f_calc_minmax( ',tmp2,',', tmp1, " , 'variable' , ", sim_id2{kk} ,' , tmp3);'))
    end
end


reac_name = {'EGP'};
for kk = 1:length(sim_id)  
    for jj = 1:length(reac_name)
    tmp1 = strcat('Max_', sim_id{kk}, '_', reac_name{jj} );
    tmp2 = strcat('Min_', sim_id{kk}, '_', reac_name{jj} );
    eval(strcat('id_',  reac_name{jj},  "= ismember(SBreactions(objModel),'", reac_name{jj},"');"));
    eval(strcat('tmp3 = id_',  reac_name{jj},';') );  
    
    if i == 1
    tmp4 = strcat( sim_id{kk} , '.reactionvalues(:, tmp3 );');
    eval(strcat( tmp1 ,'=' ,tmp4)); eval(strcat( tmp2 ,'=' ,tmp4));
    end
    eval(strcat('[ ',tmp2,' ,  ',tmp1,' ]     = f_calc_minmax( ',tmp2,',', tmp1, " , 'reaction' , ", sim_id{kk} ,' , tmp3);'))
    end
end
end

%% Plotting

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

%% Simulations for best cost.
load('./Parameters/theta_best.mat'); theta2 = Results.xbest; 

[ AS_train ] = f_as_new( theta2 ,func_mex_model, D, utility);

time_points_IPGTTs = [ 1:56 ];
[AS2_train] = f_as_IPGTTs( theta2 ,func_mex_model,D, utility, time_points_IPGTTs) ;

%% Figure4 A-C
t_IPGTT=(0:1:120);
fill_time_IPGTT = [ t_IPGTT fliplr(t_IPGTT) ];

time_fill = fill_time_IPGTT;

UidColor =  [0.9,.6,0.9];
UiiColor = [ 0.5,0.5,0.5];
%
figure (41)
set(figure(41), 'outerposition',[0 0 2500 500], 'PaperType','a4')

sim_id = {'IPGTT_chow2w','IPGTT_chow12dHFD2d','IPGTT_chow10dHFD4d','IPGTT_chow8dHFD6d','IPGTT_hfd2w','IPGTT_chow6w','IPGTT_hfd6w'};
varb_name = {'U_id','U_ii','U_id_A','U_id_M','U_id_H'};

for kk = 1:length(sim_id)
    subplot(1,7,kk)
    %title(strcat(sim_id{kk}, '-Uid and Uii'),'fontsize',title_size);
    hold on

    values=[]; values1=[];   
    eval(strcat('values(:,1)=Max_',sim_id{kk},'_U_ii;'));
    eval(strcat('values(:,2)=Min_',sim_id{kk},'_U_ii;'));
    f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], UiiColor, 'EdgeColor', 'none' );
    set(f,'facealpha',.5)
    hold on
    eval(strcat('plot( t_IPGTT, AS_train.',sim_id{kk}, ".variablevalues(:,id_U_ii),'color', UiiColor - 0.3, 'linewidth', marker_line_size,'HandleVisibility','off');"))
    
    eval(strcat('values1(:,1)=Max_',sim_id{kk},'_U_id;'));
    eval(strcat('values1(:,2)=Min_',sim_id{kk},'_U_id;'));
    f=fill(time_fill, [ values1(:,2)' fliplr( values1(:,1)' ) ], UidColor, 'EdgeColor', 'none' );
    set(f,'facealpha',.5)
    hold on
    eval(strcat('plot( t_IPGTT, AS_train.',sim_id{kk}, ".variablevalues(:,id_U_id),'color', UidColor , 'linewidth', marker_line_size,'HandleVisibility','off');"))
    
    box off
    set(gca,'FontSize',label_size_ax)
    ylabel('x','FontSize', label_size )
    xlabel('Time (min)','FontSize', label_size )
    xlim([0 120]);
    ylim([0 0.3]);
    
    if kk == 1
        legend({'insulin independent uptake', 'insulin dependent uptake'});
    end
    
end

figure (42)
set(figure(42), 'outerposition',[0 0 2500 500], 'PaperType','a4')

organ_id = {'H','M'};
UidColor2  = [1 0 0];
UidColor3  = [0 1 0];
for ii = 1:length(organ_id)
    tmp1 = cell2mat(organ_id(ii));
    if ii == 1; col = UidColor2  ; elseif ii == 2 ; col = UidColor3  ; elseif ii == 3 ;  col = UidColor4 ;  end
    for kk = 1:length(sim_id)
        values=[];values1=[];
        
        subplot(1,7,kk)      
        hold on
        values=[];
        time_fill = fill_time_IPGTT;        
        eval(strcat('values1(:,1)=Max_',sim_id{kk},'_U_id_',tmp1,';'));
        eval(strcat('values1(:,2)=Min_',sim_id{kk},'_U_id_',tmp1,';'));
        f=fill(time_fill, [ values1(:,2)' fliplr( values1(:,1)' ) ], col - col.*0.2 , 'EdgeColor', 'none' );
        set(f,'facealpha',.3)
        hold on
        eval(strcat('plot( t_IPGTT, AS_train.',sim_id{kk}, ".variablevalues(:,id_U_id_",tmp1,"),'color', col - col.*0.2  , 'linewidth', marker_line_size,'HandleVisibility','off');"))
        box off
        set(gca,'FontSize',label_size_ax)
        ylabel('X','FontSize', label_size )
        xlabel('Time (min)','FontSize', label_size )
        xlim([0 120]);
        ylim([0 0.3]);
        
        if kk == 1
            legend({'Hepatic', 'Muscle'});
        end
        
    end
    
end


figure (43)
set(figure(43), 'outerposition',[0 0 2500 500], 'PaperType','a4')

organ_id = {'A'};
UidColor  = [0 0 1]; col = UidColor;

for ii = 1:length(organ_id)
    tmp1 = cell2mat(organ_id(ii));
    for kk = 1:length(sim_id)
        values=[];values1=[];
        
        subplot(1,7,kk)
        hold on
        values=[];
        time_fill = fill_time_IPGTT;
        
        eval(strcat('values1(:,1)=Max_',sim_id{kk},'_U_id_',tmp1,';'));
        eval(strcat('values1(:,2)=Min_',sim_id{kk},'_U_id_',tmp1,';'));
        f=fill(time_fill, [ values1(:,2)' fliplr( values1(:,1)' ) ], col - col.*0.2 , 'EdgeColor', 'none' );
        set(f,'facealpha',.3)
        hold on
        eval(strcat('plot( t_IPGTT, AS_train.',sim_id{kk}, ".variablevalues(:,id_U_id_",tmp1,"),'color', col - col.*0.2  , 'linewidth', marker_line_size,'HandleVisibility','off');"))
        box off
        set(gca,'FontSize',label_size_ax)
        ylabel('X','FontSize', label_size )
        xlabel('Time (min)','FontSize', label_size )
        xlim([0 120]);
        ylim([0 0.3]);
        if kk == 1 ;         legend({'Adipose'}); end
        
        
    end
end

%%
%% Fat mass gain -for each HFD 

simulation_predColor =  [0,.9,  .5]; 
simulation_TrainColor = [ 0.5,0.5,0.5];

fmr = 0.11; lmr = 0.89; % Based on data;


figure (44)
set(figure(44), 'outerposition',[0 0 2500 1000], 'PaperType','a4')
title('Predication Fat mass gain');
subplot(1,2,1)

basal = AS_train.chow2w.statevalues(1, ismember(SBstates(objModel),'FM'));

y1  =  AS_train.IPGTT_chow2w.statevalues(1, ismember(SBstates(objModel),'FM'));% /(fmr*D.Hansson2018.BW.Mean(1,1)) ; %- fmr*D.Hansson2018.BW.Mean(1,1);

y2 =  AS_train.IPGTT_chow12dHFD2d.statevalues(1, ismember(SBstates(objModel),'FM')); %/AS_train.chow2w.statevalues(12, ismember(SBstates(objModel),'FM')) ; %- AS_train.chow2w.statevalues(12, ismember(SBstates(objModel),'FM'));

y3  =  AS_train.IPGTT_chow10dHFD4d.statevalues(1, ismember(SBstates(objModel),'FM')); %/(fmr*D.Hansson2018.BW.Mean(1,2)) ; %- fmr*D.Hansson2018.BW.Mean(1,2);

y4 =  AS_train.IPGTT_chow8dHFD6d.statevalues(1, ismember(SBstates(objModel),'FM')); %/(fmr*D.Hansson2018.BW.Mean(1,3)) ; %- fmr*D.Hansson2018.BW.Mean(1,3);

y5 =  AS_train.IPGTT_hfd2w.statevalues(1, ismember(SBstates(objModel),'FM')); %/(fmr*D.Hansson2018.BW.Mean(1,4)) ; %- fmr*D.Hansson2018.BW.Mean(1,4);

y6 =  AS_train.IPGTT_chow6w.statevalues(1, ismember(SBstates(objModel),'FM')); %/(fmr*D.NEW.BW.Mean(1,2)) ; %- fmr*D.NEW.BW.Mean(1,2);

y7 =  AS_train.IPGTT_hfd6w.statevalues(1, ismember(SBstates(objModel),'FM')); %/(fmr*D.NEW.BW.Mean(1,2)) ; %- fmr*D.NEW.BW.Mean(1,2);

yy = [ 
      y1  ;...
      y2  ;...    
      y3  ;... 
      y4  ;... 
      y5  ;... 
      NaN  ;... 
      y6  ;... 
      y7  ];  
b = bar(yy,'FaceColor','flat' ); % ,'HandleVisibility','off'
ColorFM  = [0  0  1];
b(1).CData(1,:) = simulation_predColor;
b(1).CData(2,:) = simulation_predColor;
b(1).CData(3,:) = simulation_predColor;
b(1).CData(4,:) = simulation_predColor;
b(1).CData(5,:) = simulation_predColor;
b(1).CData(6,:) = simulation_predColor;
b(1).CData(7,:) = simulation_predColor;
b(1).CData(8,:) = simulation_predColor;

set(gca,'FontSize',label_size_ax)
ylabel('Fat mass (g)','FontSize', label_size )

subplot(1,2,2)
x  =  trapz(AS_train.IPGTT_chow2w.variablevalues(:, ismember(SBvariables(objModel),'U_id_A')));
x1 =  trapz(AS_train.IPGTT_chow2w.variablevalues(:, ismember(SBvariables(objModel),'U_id_M')));
x2 =  trapz(AS_train.IPGTT_chow2w.variablevalues(:, ismember(SBvariables(objModel),'U_id_H')));
y1 = ( [ x x1 x2 ] ) * 100;

x  =  trapz(AS_train.IPGTT_chow12dHFD2d.variablevalues(:, ismember(SBvariables(objModel),'U_id_A')));
x1 =  trapz(AS_train.IPGTT_chow12dHFD2d.variablevalues(:, ismember(SBvariables(objModel),'U_id_M')));
x2 =  trapz(AS_train.IPGTT_chow12dHFD2d.variablevalues(:, ismember(SBvariables(objModel),'U_id_H')));
y2 = ( [ x x1 x2 ] ) * 100;

x  =  trapz(AS_train.IPGTT_chow10dHFD4d.variablevalues(:, ismember(SBvariables(objModel),'U_id_A')));
x1 =  trapz(AS_train.IPGTT_chow10dHFD4d.variablevalues(:, ismember(SBvariables(objModel),'U_id_M')));
x2 =  trapz(AS_train.IPGTT_chow10dHFD4d.variablevalues(:, ismember(SBvariables(objModel),'U_id_H')));
y3 = ( [ x x1 x2 ] ) * 100;

x  =  trapz(AS_train.IPGTT_chow8dHFD6d.variablevalues(:, ismember(SBvariables(objModel),'U_id_A')));
x1 =  trapz(AS_train.IPGTT_chow8dHFD6d.variablevalues(:, ismember(SBvariables(objModel),'U_id_M')));
x2 =  trapz(AS_train.IPGTT_chow8dHFD6d.variablevalues(:, ismember(SBvariables(objModel),'U_id_H')));
y4 = ( [ x x1 x2 ] ) * 100;

x  =  trapz(AS_train.IPGTT_hfd2w.variablevalues(:, ismember(SBvariables(objModel),'U_id_A')));
x1 =  trapz(AS_train.IPGTT_hfd2w.variablevalues(:, ismember(SBvariables(objModel),'U_id_M')));
x2 =  trapz(AS_train.IPGTT_hfd2w.variablevalues(:, ismember(SBvariables(objModel),'U_id_H')));
x_total = x + x1 + x2;
y5 = ( [ x x1 x2 ]) * 100;

x  =  trapz(AS_train.IPGTT_chow6w.variablevalues(:, ismember(SBvariables(objModel),'U_id_A')));
x1 =  trapz(AS_train.IPGTT_chow6w.variablevalues(:, ismember(SBvariables(objModel),'U_id_M')));
x2 =  trapz(AS_train.IPGTT_chow6w.variablevalues(:, ismember(SBvariables(objModel),'U_id_H')));
y6 = ( [ x x1 x2 ] ) * 100;

x  =  trapz(AS_train.IPGTT_hfd6w.variablevalues(:, ismember(SBvariables(objModel),'U_id_A')));
x1 =  trapz(AS_train.IPGTT_hfd6w.variablevalues(:, ismember(SBvariables(objModel),'U_id_M')));
x2 =  trapz(AS_train.IPGTT_hfd6w.variablevalues(:, ismember(SBvariables(objModel),'U_id_H')));
y7 = ( [ x x1 x2 ] ) * 100;

yy = [ y1 ; y2 ; y3 ; y4; y5 ; [ NaN  NaN NaN]; y6 ; y7 ] ;
b = bar(yy, 'stacked','FaceColor','flat' ); % ,'HandleVisibility','off'

ColorA  = [0  0  1]; 
ColorM  = [0  1  0];
ColorH  = [1  0  0];

for ii = 1:size(yy,1)
b(1).CData(ii,:) = ColorA; 
b(2).CData(ii,:) = ColorM; 
b(3).CData(ii,:) = ColorH; 
end
hold on

set(gca,'FontSize',label_size_ax)
ylabel('AUC glucose uptake (mmol)','FontSize', label_size )
box off
legend({' Adipose', 'Muscle', 'Liver'}, 'Location','southwest'  );
xticklabels({'2 w control','2 d HFD','4 d HFD','6 d HFD','2 w HFD', '', '6 w control','6 w HFD'})
xtickangle(90)

%% ADIPOSE TISSUE - change over 8 weeks
figure (45)
set(figure(45), 'outerposition',[0 0 1500 1000], 'PaperType','a4')

fmr=0.11;

for kk = 1:length(time_points_IPGTTs)
    eval(strcat('xxx1(kk,:)= AS2_train.IPGTT_hfd', int2str(time_points_IPGTTs(kk)) ,'d.statevalues(1, ismember(SBstates(objModel), "FM")) - (fmr*D.Hansson2018.BW.Mean(1,4));'));    
    eval(strcat('xxx2(kk,:)= AS2_train.IPGTT_chow', int2str(time_points_IPGTTs(kk)) ,'d.statevalues(1, ismember(SBstates(objModel), "FM")) - (fmr*D.Hansson2018.BW.Mean(1,1));'));    
    if kk == 1
    eval(strcat('basal1= trapz(AS2_train.IPGTT_hfd', int2str(time_points_IPGTTs(kk)) ,'d.variablevalues(:, ismember(SBvariables(objModel), "U_id_A")));')); 
    %eval(strcat('basal2= trapz(AS2_train.IPGTT_hfd', int2str(time_points_IPGTTs(kk)) ,'d.variablevalues(:, ismember(SBvariables(objModel), "U_id_A_no_fat")));')); 
    eval(strcat('basal3= trapz(AS2_train.IPGTT_hfd', int2str(time_points_IPGTTs(kk)) ,'d.variablevalues(:, ismember(SBvariables(objModel), "U_id_A_no_glucose")));'));        
    eval(strcat('basal1_chow= trapz(AS2_train.IPGTT_chow', int2str(time_points_IPGTTs(kk)) ,'d.variablevalues(:, ismember(SBvariables(objModel), "U_id_A")));')); 
    end       
    eval(strcat('yyy1(kk,:)= trapz(AS2_train.IPGTT_hfd', int2str(time_points_IPGTTs(kk)) ,'d.variablevalues(:, ismember(SBvariables(objModel), "U_id_A")))./basal1;'));         
    %eval(strcat('yyy2(kk,:)= trapz(AS2_train.IPGTT_hfd', int2str(time_points_IPGTTs(kk)) ,'d.variablevalues(:, ismember(SBvariables(objModel), "U_id_A_no_fat")))/basal2;'));     
    eval(strcat('yyy3(kk,:)= trapz(AS2_train.IPGTT_hfd', int2str(time_points_IPGTTs(kk)) ,'d.variablevalues(:, ismember(SBvariables(objModel), "U_id_A_no_glucose")))/basal3;'));    
    eval(strcat('zzz(kk,:)= trapz(AS2_train.IPGTT_chow', int2str(time_points_IPGTTs(kk)) ,'d.variablevalues(:, ismember(SBvariables(objModel), "U_id_A")))./basal1_chow;'));         
end

Color1  = [1  1  1];
Color2  = [0  .4 1];
Color3  = [.1 .7 1];
Color4  = [0 1 0];

tx = 1:1:length(time_points_IPGTTs);

plot( tx, ones(1,length(tx)), 'k','LineWidth',1) % , 'HandleVisibility','off'
hold on
sz = 100;

scatter(tx,yyy1,sz,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Color1,'LineWidth',1)
scatter(tx,yyy3,sz,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Color3,'LineWidth',1)

legend('1 = no change',...
       'HFD - Uptake in model',...
       'HFD - Isolated Adipose glucose uptake ', 'Location','northwest'    );
xlim([ 0 56]); xticks([ 0 56]); xticklabels([ 0 8]);
ylim([ 0  6]); yticks([ 0 1 2 3 4 5 6]); 
set(gca,'FontSize',label_size_ax)
ylabel('AUC adipose glucose uptake (fold change)','FontSize', label_size )
xlabel('Time (weeks)','FontSize', label_size )
box off

%% SA Figure
figure (46)
set(figure(46), 'outerposition',[0 0 2000 1500], 'PaperType','a4')
varb_name = { 'PlotEGP', 'PlotGdiff' , 'PlotGIP'   ,'PlotU'       };

sim_id = {'IPGTT_chow2w','IPGTT_hfd2w','IPGTT_hfd6w'};

for kk = 1:length(sim_id)
    subplot(4,3,kk)
    title(strcat(sim_id{kk}),'fontsize',title_size);
    hold on

    values=[]; values1=[];   
    eval(strcat('values(:,1)=Max_',sim_id{kk},'_PlotEGP;'));
    eval(strcat('values(:,2)=Min_',sim_id{kk},'_PlotEGP;'));
    f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], UiiColor, 'EdgeColor', 'none' );
    set(f,'facealpha',.5)
    hold on
    eval(strcat('plot( t_IPGTT, AS_train.',sim_id{kk}, ".variablevalues(:,id_PlotEGP),'color', UiiColor - 0.3, 'linewidth', marker_line_size,'HandleVisibility','off');"))
    
    box off
    set(gca,'FontSize',label_size_ax)
    ylabel('EGP (mol/h)','FontSize', label_size )
    xlabel('Time (min)','FontSize', label_size )
    xlim([0 120]);
    ylim([0 0.14]);

end

for kk = 1:length(sim_id)
    subplot(4,3,kk + 3)
    hold on

    values=[]; values1=[];   
    eval(strcat('values(:,1)=Max_',sim_id{kk},'_PlotU;'));
    eval(strcat('values(:,2)=Min_',sim_id{kk},'_PlotU;'));
    f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], UiiColor, 'EdgeColor', 'none' );
    set(f,'facealpha',.5)
    hold on
    eval(strcat('plot( t_IPGTT, AS_train.',sim_id{kk}, ".variablevalues(:,id_PlotU),'color', UiiColor - 0.3, 'linewidth', marker_line_size,'HandleVisibility','off');"))
    
    box off
    set(gca,'FontSize',label_size_ax)
    ylabel('U (mol/h)','FontSize', label_size )
    xlabel('Time (min)','FontSize', label_size )
    xlim([0 120]);
    ylim([0 0.4]);

    %
end


for kk = 1:length(sim_id)
    subplot(4,3,kk + 6)
    hold on

    values=[]; values1=[];   
    eval(strcat('values(:,1)=Max_',sim_id{kk},'_PlotGdiff;'));
    eval(strcat('values(:,2)=Min_',sim_id{kk},'_PlotGdiff;'));
    f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], UiiColor, 'EdgeColor', 'none' );
    set(f,'facealpha',.5)
    hold on
    eval(strcat('plot( t_IPGTT, AS_train.',sim_id{kk}, ".variablevalues(:,id_PlotGdiff),'color', UiiColor - 0.3, 'linewidth', marker_line_size,'HandleVisibility','off');"))
    
    box off
    set(gca,'FontSize',label_size_ax)
    ylabel('G_Diff (mol/h)','FontSize', label_size )
    xlabel('Time (min)','FontSize', label_size )
    xlim([0 120]);
    %ylim([0 14]);

end


for kk = 1:length(sim_id)
    subplot(4,3,kk + 9)
    hold on

    values=[]; values1=[];   
    eval(strcat('values(:,1)=Max_',sim_id{kk},'_PlotGIP;'));
    eval(strcat('values(:,2)=Min_',sim_id{kk},'_PlotGIP;'));
    f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], UiiColor, 'EdgeColor', 'none' );
    set(f,'facealpha',.5)
    hold on
    eval(strcat('plot( t_IPGTT, AS_train.',sim_id{kk}, ".variablevalues(:,id_PlotGIP),'color', UiiColor - 0.3, 'linewidth', marker_line_size,'HandleVisibility','off');"))
    
    box off
    set(gca,'FontSize',label_size_ax)
    ylabel('G_I_P (mol/h)','FontSize', label_size )
    xlabel('Time (min)','FontSize', label_size )
    xlim([0 120]);
    %ylim([0 14]);

end

