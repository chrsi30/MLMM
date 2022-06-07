
for i = 1:nfiles
load( strcat( filedir ,'/',files(i).name ) )
disp(strcat('Processing - ',files(i).name, int2str(i),'/',int2str(nfiles)))
theta=Results.xbest;

time_points_IPGTTs = 1:1:7*4;
[AS] = f_as_IPGTTs( Results.xbest ,func_mex_model,D, utility, time_points_IPGTTs) ;

for kk = 1:length(time_points_IPGTTs)  
    tmp = strcat('IPGTT_chow', int2str(time_points_IPGTTs(kk)),'d  = AS.IPGTT_chow', int2str(time_points_IPGTTs(kk)),'d ;');
    eval(tmp);    
    tmp1 = strcat('IPGTT_hfd', int2str(time_points_IPGTTs(kk)),'d  = AS.IPGTT_hfd', int2str(time_points_IPGTTs(kk)),'d ;');
    eval(tmp1);
end


diet_chow = AS.diet_chow;
diet_hfd  =AS.diet_hfd;                


%% MIN and MAX value calculations

%%

varb_name1 = {'FatMass',...
              'Gconc',...
              'Iconc',...
              'PlotEGP', ...
              'ISEC',...
              'U_id_A',...
              'U_id_H',...
              'U_id_M',...
              'PlotRABGTP',...
              'PlotPMAa',...
              'PlotPKBda',...
              'PlotPKBda',...
              'f_IR_CLGI', ...
              'f_IR_EGP', ...
              'f_IR_Ins', ...
              'f_IR_PMA' , ...
              'f_IR_IRS_tot_prot'};

sim_id1 = {'diet_chow' , 'diet_hfd' };

for kk = 1:length(sim_id1)
    for jj = 1:length(varb_name1)
        tmp1 = strcat('Max_', sim_id1{kk}, '_', varb_name1{jj} );
        tmp2 = strcat('Min_', sim_id1{kk}, '_', varb_name1{jj} );
        eval(strcat('id_',  varb_name1{jj},  "= ismember(SBvariables(objModel),'", varb_name1{jj},"');"));
        eval(strcat('tmp3 = id_',  varb_name1{jj},';') );
        
        if i == 1
            tmp4 = strcat( sim_id1{kk} , '.variablevalues(:, tmp3 );');
            eval(strcat( tmp1 ,'=' ,tmp4)); eval(strcat( tmp2 ,'=' ,tmp4));
        end
        eval(strcat('[ ',tmp2,' ,  ',tmp1,' ]     = f_calc_minmax( ',tmp2,',', tmp1, " , 'variable' , ", sim_id1{kk} ,' , tmp3);'))
    end
end          
          
          
%% Short term

varb_name2 = {'Gconc',...
              'Iconc',...
              'U_id',...
              'PlotEGP'};
sim_id2 = { 'IPGTT_chow1d', 'IPGTT_chow7d', 'IPGTT_chow14d', 'IPGTT_chow21d', 'IPGTT_chow28d' };
          
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

sim_id3 = { 'IPGTT_hfd1d', 'IPGTT_hfd7d', 'IPGTT_hfd14d', 'IPGTT_hfd21d', 'IPGTT_hfd28d' };

for kk = 1:length(sim_id3)
    for jj = 1:length(varb_name2)
        tmp1 = strcat('Max_', sim_id3{kk}, '_', varb_name2{jj} );
        tmp2 = strcat('Min_', sim_id3{kk}, '_', varb_name2{jj} );
        eval(strcat('id_',  varb_name2{jj},  "= ismember(SBvariables(objModel),'", varb_name2{jj},"');"));
        eval(strcat('tmp3 = id_',  varb_name2{jj},';') );
        
        if i == 1
            tmp4 = strcat( sim_id3{kk} , '.variablevalues(:, tmp3 );');
            eval(strcat( tmp1 ,'=' ,tmp4)); eval(strcat( tmp2 ,'=' ,tmp4));
        end
        eval(strcat('[ ',tmp2,' ,  ',tmp1,' ]     = f_calc_minmax( ',tmp2,',', tmp1, " , 'variable' , ", sim_id3{kk} ,' , tmp3);'))
    end
end



end

%% Simulations for best cost.
load('./Parameters/theta_best.mat'); theta2 = Results.xbest; 

time_points_IPGTTs = 1:1:7*4;
[AS2] = f_as_IPGTTs( Results.xbest ,func_mex_model,D, utility, time_points_IPGTTs) ;

%% Plotting

label_size       = 16 ;
label_size_ax    = 16 ;
marker_line_size = 1.5;
marker_size      = 5  ;
line_size        = 1.5;
title_size       = 18 ;

control_color   = [ 0,   .5,  1]; 
HFD_color       = [ 1,  .6,   .1 ]; 


fig=0;
%% Long term


varb_name1 = {'FatMass',...
    'Gconc',...
    'Iconc',...
    'PlotEGP', ...
    'ISEC',...
    'U_id_A',...
    'U_id_H',...
    'U_id_M',...
    'PlotRABGTP',...
    'PlotPMAa',...
    'PlotPKBda',...
    'PlotPKBda',...
    'f_IR_CLGI', ...
    'f_IR_EGP', ...
    'f_IR_Ins'}; % 'f_IR_CLGI', ...              'f_IR_PMA' , ... 'f_IR_IRS_tot_prot'
              

sim_id1 = {'diet_chow' , 'diet_hfd' };

sim_res = 0.1;

time_fill = [ 0:sim_res:time_points_IPGTTs(end) fliplr(0:sim_res:time_points_IPGTTs(end)) ] ;

figure (61)
set(figure(61), 'outerposition',[0 0 1500 1000], 'PaperType','a4')

for jj = 1:length(varb_name1)
     subplot(5,3,jj)
    for kk = 1:length(sim_id1)
        
        if kk == 1 ; col = control_color; elseif  kk == 2; col = HFD_color;  end
        
        values=[]; values1=[];
        eval(strcat('values(:,1)=Max_', sim_id1{kk},'_',varb_name1{jj},';'));
        eval(strcat('values(:,2)=Min_', sim_id1{kk},'_',varb_name1{jj},';'));
        f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], col, 'EdgeColor', 'none' );
        set(f,'facealpha',.5)
        hold on
        
        eval(strcat('plot( AS2.', sim_id1{kk},'.time, AS2.',sim_id1{kk}, ".variablevalues(:,id_", varb_name1{jj} ,"),'color', col - col.*0.1, 'linewidth', marker_line_size,'HandleVisibility','off');"))
        
        set(gca,'FontSize',label_size_ax)
        ylabel(varb_name1{jj},'FontSize', label_size )
        xlabel('Time (min)','FontSize', label_size )
        
    end
    
end




%%

varb_name1 = {'f_IR_CLGI', ...
              'f_IR_PMA' , ...
              'f_IR_IRS_tot_prot'};

sim_id1 = {'diet_chow' , 'diet_hfd' };

sim_res = 0.1;


figure (62)
set(figure(62), 'outerposition',[0 0 1000 1000], 'PaperType','a4')


for jj = 1:length(varb_name1)
    subplot(3,1,jj)
    
    for kk = 1:length(sim_id1)
        
        if kk == 1 ; col = control_color; elseif  kk == 2; col = HFD_color;  end
        
        values=[]; values1=[];
        eval(strcat('values(:,1)=Max_', sim_id1{kk},'_',varb_name1{jj},';'));
        eval(strcat('values(:,2)=Min_', sim_id1{kk},'_',varb_name1{jj},';'));
        f=fill(time_fill, [ 1./values(:,2)' fliplr( 1./values(:,1)' ) ], col, 'EdgeColor', 'none' );
        set(f,'facealpha',.5)
        hold on
        
        %eval(strcat('plot( AS2.', sim_id1{kk},'.time, AS2.',sim_id1{kk}, ".variablevalues(:,id_", varb_name1{jj} ,"),'color', col - col.*0.1, 'linewidth', marker_line_size,'HandleVisibility','off');"))
        
        set(gca,'FontSize',label_size_ax)
        ylabel(varb_name1{jj},'FontSize', label_size )
        xlabel('Time (min)','FontSize', label_size )
        
    end

end









%% Short term
varb_name2 = {'Gconc',...
              'Iconc',...
              'U_id',...
              'PlotEGP'};
          
sim_id2 = { 'IPGTT_chow1d', 'IPGTT_chow7d', 'IPGTT_chow14d', 'IPGTT_chow21d', 'IPGTT_chow28d' };

sim_id3 = { 'IPGTT_hfd1d', 'IPGTT_hfd7d', 'IPGTT_hfd14d', 'IPGTT_hfd21d', 'IPGTT_hfd28d' };


t_IPGTT   = 0:1:120 ; % min

time_fill = [ t_IPGTT fliplr(t_IPGTT) ] ;

figure (63)
set(figure(63), 'outerposition',[0 0 1500 1000], 'PaperType','a4')


for jj = 1:length(varb_name2)
    for kk = 1:length(sim_id2)
        subplot(length(varb_name2),length(sim_id2), length(sim_id2)*(jj-1) + kk )
        
        hold on 
        values=[]; values1=[];
        eval(strcat('values(:,1)=Max_', sim_id2{kk},'_',varb_name2{jj},';'));
        eval(strcat('values(:,2)=Min_', sim_id2{kk},'_',varb_name2{jj},';'));
        f=fill(time_fill, [ values(:,2)' fliplr( values(:,1)' ) ], control_color, 'EdgeColor', 'none' );
        set(f,'facealpha',.5)
        hold on
        
        eval(strcat('plot( t_IPGTT, AS2.',sim_id2{kk}, ".variablevalues(:,id_", varb_name2{jj} ,"),'color', control_color - control_color.*0.1, 'linewidth', marker_line_size,'HandleVisibility','off');"))
        
        
        values1=[];
        eval(strcat('values1(:,1)=Max_', sim_id3{kk},'_',varb_name2{jj},';'));
        eval(strcat('values1(:,2)=Min_', sim_id3{kk},'_',varb_name2{jj},';'));
        f=fill(time_fill, [ values1(:,2)' fliplr( values1(:,1)' ) ], HFD_color, 'EdgeColor', 'none' );
        set(f,'facealpha',.5)
        hold on
        
        eval(strcat('plot( t_IPGTT, AS2.',sim_id3{kk}, ".variablevalues(:,id_", varb_name2{jj} ,"),'color', HFD_color - HFD_color.*0.1, 'linewidth', marker_line_size,'HandleVisibility','off');"))
               
        set(gca,'FontSize',label_size_ax)
        ylabel(varb_name2{jj},'FontSize', label_size )
        xlabel('Time (min)','FontSize', label_size )
        if jj == 1
            ylim([ 5 36]);
            ylabel('Glucose')
        elseif jj == 2
            ylim([ 50 1100]);
            ylabel('Inuslin')
        elseif jj == 3
            ylim([ 0 0.35]);
            ylabel('Glucose uptake')
        elseif jj == 4
            ylim([ 0 0.12]);
            ylabel('EGP')
        end

    
    end
end



