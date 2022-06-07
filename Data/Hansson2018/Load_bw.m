% BW data from Hansson 2018.
% Excel file C:\Users\chrsi30\Documents\System biologi\Projekt\MLMM\Manus1\Data\RAW DATA\Hansson 2018

%% 2w chow
chow2w_start=[23	23	25	25	25	22	28	25	22	21	22.5	22.6	22.9	28.3	22.5	22.5	23.5	20.8	23.8	22.7	24.7	23.5	23.7	22.2	23.5	22.7	21	26.7	25	22.8 ];
chow2w_start_MEAN = mean(chow2w_start);
chow2w_start_SEM = std(chow2w_start)/sqrt(length(chow2w_start));

chow2w_end   = [ 25.1	23.6	26.6	27.2	23.5	23.7	28.2	22.4	22.5	24.7	29.2	24.2	23.6	23.2	25.2	21.8	23.5	25.1	22.5	23.8	24.7	25.7	25	26.6	22.1	25.2	22.9	24.4	29.4];
chow2w_end_MEAN = mean(chow2w_end);
chow2w_end_SEM = std(chow2w_end)/sqrt(length(chow2w_end));

chow_2w_bw_MEAN = [ chow2w_start_MEAN ; chow2w_end_MEAN];
chow_2w_bw_SEM = [ chow2w_start_SEM ; chow2w_end_SEM];


%% 10 d chow & 4 d HFD 
chow_10d_hfd_4d_start      = [ 23.1	24.4	22.3	21.1	25.7	24.6	23.6	22	25.7	23.8	25.2	23.3	24.3	23.7	22.3	22.3	24.4 ];
chow_10d_hfd_4d_start_MEAN = mean(chow_10d_hfd_4d_start);
chow_10d_hfd_4d_start_SEM  = std(chow_10d_hfd_4d_start)/sqrt(length(chow_10d_hfd_4d_start)) ;

chow_10d_hfd_4d_end = [26	23.4	27.9	26.1	27.5	27	27.8	30.4	25	27.3	32.5	27.1	27	26.8	24.4	25.3	26	27.1];
chow_10d_hfd_4d_end_MEAN = mean(chow_10d_hfd_4d_end);
chow_10d_hfd_4d_end_SEM = std(chow_10d_hfd_4d_end)/sqrt(length(chow_10d_hfd_4d_end));

chow_10d_hfd_4d_bw_MEAN = [ chow_10d_hfd_4d_start_MEAN ; chow_10d_hfd_4d_end_MEAN];
chow_10d_hfd_4d_bw_SEM = [ chow_10d_hfd_4d_start_SEM ; chow_10d_hfd_4d_end_SEM];


%% 8 d chow & 6 d HFD 
chow_8d_hfd_6d_start = [ 27.1	21.6	24.7	26.8	26	24.3	23.2	22.9	23.3	22.6	24.4	23.3	23.5	24	24.7 ];
chow_8d_hfd_6d_start_MEAN = mean(chow_8d_hfd_6d_start);
chow_8d_hfd_6d_start_SEM = std(chow_8d_hfd_6d_start)/sqrt(length(chow_8d_hfd_6d_start)) ;

chow_8d_hfd_6d_end = [33.2	30.8	23.7	32	28.1	26.7	26.5	26.8	26.4	29	26.5	27	28.5	28.5	27.1 ];
chow_8d_hfd_6d_end_MEAN = mean(chow_8d_hfd_6d_end);
chow_8d_hfd_6d_end_SEM = std(chow_8d_hfd_6d_end)/sqrt(length(chow_8d_hfd_6d_end));

chow_8d_hfd_6d_bw_MEAN = [ chow_8d_hfd_6d_start_MEAN ; chow_8d_hfd_6d_end_MEAN];
chow_8d_hfd_6d_bw_SEM = [ chow_8d_hfd_6d_start_SEM ; chow_8d_hfd_6d_end_SEM];

%% 2 w hfd

hfd2w_start =[ 23	22	23	24	21.7	22.6	24.4	26.5	22	21.5	22 ];
hfd2w_start_MEAN = mean(hfd2w_start);
hfd2w_start_SEM = std(hfd2w_start)/sqrt(length(hfd2w_start)) ;

hfd2w_end = [29.2	30.9	31.2	28.5	29	30	31.5	30.4	26.9	24	26.2 ];
hfd2w_end_MEAN = mean(hfd2w_end);
hfd2w_end_SEM = std(hfd2w_end)/sqrt(length(hfd2w_end));

hfd2w_bw_MEAN = [ hfd2w_start_MEAN ; hfd2w_end_MEAN];
hfd2w_bw_SEM = [ hfd2w_start_SEM ; hfd2w_end_SEM];

%%

bw_MEAN = [chow_2w_bw_MEAN , chow_10d_hfd_4d_bw_MEAN , chow_8d_hfd_6d_bw_MEAN ,  hfd2w_bw_MEAN ];
bw_SEM = [chow_2w_bw_SEM , chow_10d_hfd_4d_bw_SEM , chow_8d_hfd_6d_bw_SEM ,  hfd2w_bw_SEM ];

Hanson_2018_BW.Mean = bw_MEAN;
Hanson_2018_BW.SEM = bw_SEM;
Hanson_2018_BW.time = [ [0 14 ]' , [10 14]', [ 8 14]', [ 0 14]'  ];

save('Hanson_2018_BW','Hanson_2018_BW')

% Hanson_2018_BW = table( [ 0 14;10 14; 8 14;0 14],...
%     [ chow2w_start_MEAN  chow2w_end_MEAN; hfd4d_start_MEAN hfd4d_end_MEAN; chow_8d_hfd_6d_start_MEAN chow_8d_hfd_6d_end_MEAN ; hfd2w_start_MEAN hfd2w_end_MEAN],...
%     [ chow2w_start_SEM  chow2w_end_SEM;...
%       hfd4d_start_SEM hfd4d_end_SEM;...
%       chow_8d_hfd_6d_start_SEM chow_8d_hfd_6d_end_SEM ; hfd2w_start_SEM hfd2w_end_SEM],...
%     'VariableNames', {' time', 'mean', 'SEM' } ) ;
% 
% save('Hanson_2018_BW','Hanson_2018_BW')
% figure(1)
% errorbar( Hanson_2018_BW.time(1,:),  Hanson_2018_BW.mean(1,:),  Hanson_2018_BW.SEM(1,:) )
% hold on
% errorbar( Hanson_2018_BW.time(2,:),  Hanson_2018_BW.mean(2,:),  Hanson_2018_BW.SEM(2,:) )
% errorbar( Hanson_2018_BW.time(3,:),  Hanson_2018_BW.mean(3,:),  Hanson_2018_BW.SEM(3,:) )
% errorbar( Hanson_2018_BW.time(4,:),  Hanson_2018_BW.mean(4,:),  Hanson_2018_BW.SEM(4,:) )
% legend({'chow','4 d hfd', '6 d hfd' , '14 d hfd'})



