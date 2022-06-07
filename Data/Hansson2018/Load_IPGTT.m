%%% Hansson 2018 et al
% From excel file 

chow_2w_IPGTT_glucose = [6       6.1     5.7     5.6;...
    20.3	19.2	18.9	20.9;...
    19.2	16.6	17.5	16.3;...
    10.6	10.7	10.1	10.2;...
    8.6     7.3     5.5     7 ];

chow_2w_IPGTT_glucose_MEAN = mean(chow_2w_IPGTT_glucose,2);
n = size(chow_2w_IPGTT_glucose,2) ; 
chow_2w_IPGTT_glucose_SEM = std(chow_2w_IPGTT_glucose,0,2)./sqrt(n) ;

chow_12d_hfd_2d_IPGTT_glucose = [7.6     7.5     8.7     8.3     7.1;...
    28.5	22.8	24.7	22.7	22.6;...
    29.4	16.2	17.6	20.8	19.9;...
    21.7	13      11.8	14.1	11.7;...
    8.3	9.6     7.3     6.6     7.1 ];

chow_12d_hfd_2d_IPGTT_glucose_MEAN = mean(chow_12d_hfd_2d_IPGTT_glucose,2);
n = size(chow_12d_hfd_2d_IPGTT_glucose,2) ; 
chow_12d_hfd_2d_IPGTT_glucose_SEM = std(chow_12d_hfd_2d_IPGTT_glucose,0,2)./sqrt(n) ;


chow_8d_hfd_6d_IPGTT_glucose  = [ 9      7.8     8       8       6.7;...
    25.1	33.2	24.9	31.7	22.8;...
    25.8	30.8	24.7	28      19.4;...
    25.4	13.9	18.7	15.4	14.1;...
    14.1	7.6     7.7 	7.7 	7.4 ];

chow_8d_hfd_6d_IPGTT_glucose_MEAN = mean(chow_8d_hfd_6d_IPGTT_glucose,2);
n = size(chow_8d_hfd_6d_IPGTT_glucose,2) ; 
chow_8d_hfd_6d_IPGTT_glucose_SEM = std(chow_8d_hfd_6d_IPGTT_glucose,0,2)./sqrt(n) ;


hfd_2w_IPGTT_glucose  = [   9       9.2     9.1     11.8     9.1 ;...
    25.8    27.9    23.2    25.3     25.5 ;...
    27.3    30.2    24.9    26       26.1 ;...
    25.7    28.4    19.8    23.8     23.7 ;...
    15.3    14      10.6    10.6     13.1 ];

hfd_2w_IPGTT_glucose_MEAN = mean(hfd_2w_IPGTT_glucose,2);
n = size(hfd_2w_IPGTT_glucose,2) ; 
hfd_2w_IPGTT_glucose_SEM = std(hfd_2w_IPGTT_glucose,0,2)./sqrt(n) ;

glucose_mean = [chow_2w_IPGTT_glucose_MEAN , chow_12d_hfd_2d_IPGTT_glucose_MEAN , chow_8d_hfd_6d_IPGTT_glucose_MEAN ,  hfd_2w_IPGTT_glucose_MEAN ];
glucose_SEM = [chow_2w_IPGTT_glucose_SEM , chow_12d_hfd_2d_IPGTT_glucose_SEM , chow_8d_hfd_6d_IPGTT_glucose_SEM ,  hfd_2w_IPGTT_glucose_SEM];

%Normalize SEM
glucose_SEM_n = glucose_SEM;
glucose_SEM_n( glucose_SEM_n(:,1) < mean(glucose_SEM_n(:,1)) , 1 ) = mean(glucose_SEM_n(:,1)); %14d chow
glucose_SEM_n( glucose_SEM_n(:,2) < mean(glucose_SEM_n(:,2)) , 2 ) = mean(glucose_SEM_n(:,2)); %2d HFD
glucose_SEM_n( glucose_SEM_n(:,3) < mean(glucose_SEM_n(:,3)) , 3 ) = mean(glucose_SEM_n(:,3));%6d HFD
glucose_SEM_n( glucose_SEM_n(:,4) < mean(glucose_SEM_n(:,4)) , 4 ) = mean(glucose_SEM_n(:,4)); %14d
glucose_SEM = glucose_SEM_n;


Hanson_2018_IPGTT.Glucose.Mean = glucose_mean;
Hanson_2018_IPGTT.Glucose.SEM  = glucose_SEM;
Hanson_2018_IPGTT.Glucose.time = [ 0 15 30 60 120];

% 
% % figure(3)
% subplot(1,4,1)
% errorbar( [1:5] , glucose_mean(:,1),  glucose_SEM(:,1), 'k*' , 'linewidth', 1.5 )
% hold on
% errorbar( [2:6] ,glucose_mean(:,1),  glucose_SEM_n(:,1), 'bO' , 'linewidth', 1.5 )
% 
% subplot(1,4,2)
% errorbar( [1:5] , glucose_mean(:,2),  glucose_SEM(:,2), 'k*' , 'linewidth', 1.5 )
% hold on
% errorbar( [2:6] ,glucose_mean(:,2),  glucose_SEM_n(:,2), 'bO' , 'linewidth', 1.5 )
% 
% subplot(1,4,3)
% errorbar( [1:5] , glucose_mean(:,3),  glucose_SEM(:,3), 'k*' , 'linewidth', 1.5 )
% hold on
% errorbar( [2:6] ,glucose_mean(:,3),  glucose_SEM_n(:,3), 'bO' , 'linewidth', 1.5 )
% 
% subplot(1,4,4)
% errorbar( [1:5] , glucose_mean(:,4),  glucose_SEM(:,4), 'k*' , 'linewidth', 1.5 )
% hold on
% errorbar( [2:6] ,glucose_mean(:,4),  glucose_SEM_n(:,4), 'bO' , 'linewidth', 1.5 )

%%
%%%% insulin
conversion_Factor=0.005;%  plasma insulin pmol/L to micro µg/L    -  1pmol/L = 0.005µg/L

chow_2w_IPGTT_insulin         = [0.80833	0.491823	0.2         0.2         0.2 ;...
    2.035878  NaN         0.9966778   0.95901     1.180612 ;...
    0.811651  0.65601     1.012866    NaN         1.059292 ]./conversion_Factor;

chow_2w_IPGTT_insulin_MEAN = nanmean(chow_2w_IPGTT_insulin,2);
n = size(chow_2w_IPGTT_insulin,2) ;
chow_2w_IPGTT_insulin_SEM = nanstd(chow_2w_IPGTT_insulin,0,2)./sqrt(n) ;

chow_12d_hfd_2d_IPGTT_insulin = [ 1.404895	2.19        1.117579	0.81523     0.250913 ;...
    NaN      3.530257    2.668546    1.958832    1.569577  ;...
    0.94861  1.605593    2.668546    1.958832     1.569577 ]./conversion_Factor; % <-- OBS! Måste vara fel

chow_12d_hfd_2d_IPGTT_insulin_MEAN = nanmean(chow_12d_hfd_2d_IPGTT_insulin,2);
n = size(chow_12d_hfd_2d_IPGTT_insulin,2) ;
chow_12d_hfd_2d_IPGTT_insulin_SEM = nanstd(chow_12d_hfd_2d_IPGTT_insulin,0,2)./sqrt(n) ;

chow_8d_hfd_6d_IPGTT_insulin  = [1.839022  2.661749	2.40482     1.063534	0.400395 ;...
    2.150062  2.319392    2.413649    2.479174    2.291834 ;...
    3.641315  2.177129    1.83        1.645       1.08     ]./conversion_Factor;

chow_8d_hfd_6d_IPGTT_insulin_MEAN = nanmean(chow_8d_hfd_6d_IPGTT_insulin,2);
n = size(chow_8d_hfd_6d_IPGTT_insulin,2) ;
chow_8d_hfd_6d_IPGTT_insulin_SEM = nanstd(chow_8d_hfd_6d_IPGTT_insulin,0,2)./sqrt(n) ;


hfd_2w_IPGTT_insulin         = [2.922635	2.318197	1.361193	1.089486	2.489375 ;...
    3.043741  2.551524    2.422465    1.755598    3.851014 ;...
    6.121039  3.067287    2.614557    2.196725    4.52     ]./conversion_Factor;

hfd_2w_IPGTT_insulin_MEAN = nanmean(hfd_2w_IPGTT_insulin,2);
n = size(hfd_2w_IPGTT_insulin,2) ;
hfd_2w_IPGTT_insulin_SEM = nanstd(chow_8d_hfd_6d_IPGTT_insulin,0,2)./sqrt(n) ;


insulin_mean = [chow_2w_IPGTT_insulin_MEAN , chow_12d_hfd_2d_IPGTT_insulin_MEAN , chow_8d_hfd_6d_IPGTT_insulin_MEAN ,  hfd_2w_IPGTT_insulin_MEAN ];
insulin_SEM = [chow_2w_IPGTT_insulin_SEM , chow_12d_hfd_2d_IPGTT_insulin_SEM , chow_8d_hfd_6d_IPGTT_insulin_SEM ,  hfd_2w_IPGTT_insulin_SEM];

Hanson_2018_IPGTT.Insulin.Mean = insulin_mean;
Hanson_2018_IPGTT.Insulin.SEM  = insulin_SEM;
Hanson_2018_IPGTT.Insulin.time = [ 0 15 120];

% figure(2)
% errorbar( [1:4] , Stenkula_IPGTT_icc_Mean(1,:),  Stenkula_IPGTT_icc_SEM(1,:), 'r*', 'linewidth', 1.5 )
% hold on
% errorbar( [1.1:4.1] , insulin_mean(1,:),  insulin_SEM(1,:), 'k*' , 'linewidth', 1.5 )
% legend({'old','new'});
% xticklabels({'2wChow', '2dHFD','6dHFD','2wHFD'});
% xticks([1:5]);
% 
% 
% figure(3)
% subplot(1,4,1)
% errorbar( [1:3] , insulin_mean(:,1),  insulin_SEM(:,1), 'k*' , 'linewidth', 1.5 )
% subplot(1,4,2)
% errorbar( [1:3] , insulin_mean(:,2),  insulin_SEM(:,2), 'k*' , 'linewidth', 1.5 )
% subplot(1,4,3)
% errorbar( [1:3] , insulin_mean(:,3),  insulin_SEM(:,3), 'k*' , 'linewidth', 1.5 )
% subplot(1,4,4)
% errorbar( [1:3] , insulin_mean(:,4),  insulin_SEM(:,4), 'k*' , 'linewidth', 1.5 )
% 
% 
save('Hanson_2018_IPGTT','Hanson_2018_IPGTT');


