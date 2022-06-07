%%% NEW DATA 6w IPGTT

                              
chow_6w_IPGTT  =  [5.8	6.3     6.8     6.7     5       6.4     4.7     6       7.9     6.8 ;...
                                  26.1	31.3	25.8	22.4	25.2	22.6	21.9	24.8	20.1	20.3;...
                                  23.3	28.5	21.7	NaN     NaN 	16.6	20.3	22.9	17      17.5;...
                                  13	17.1	17.2	16.5	11.4	9.9     9.9     15.9	12.7	13.9;...
                                  7.5	9.7     8       8.7     7.7     7.2     6.9     8.7     8       9.2 ];

chow_6w_IPGTT_glucose_MEAN = nanmean(chow_6w_IPGTT,2);
n = size(chow_6w_IPGTT,2) ; 
chow_6w_IPGTT_glucose_SEM = nanstd(chow_6w_IPGTT,0,2)./sqrt(n) ;
           


hfd_6w_IPGTT  =  [ 8.3	9.6     11.4	9.7     9.6     8       9.6     10.1	12.2	13   ;...
                                  28.7	29.2	27.8	26.7	30.2	25.7	24.4	26.6	33.2	32.5 ;...
                                  29.8	33.3    NaN     NaN     NaN     29.7	26.4	30.6	33.1	33.3 ;...
                                  29.9	33.3	23      22.9	30.2	25.5	20      29      31.5	28.2 ;...
                                  20.6	30.1	12.9	11.4	20.4	20      10.9	19.9	20.7	15.5 ];

hfd_6w_IPGTT_glucose_MEAN = nanmean(hfd_6w_IPGTT,2);
n = size(hfd_6w_IPGTT,2) ; 
hfd_6w_IPGTT_glucose_SEM = nanstd(hfd_6w_IPGTT,0,2)./sqrt(n) ;
        
NEWDATA_IPGTT.Mean = [chow_6w_IPGTT_glucose_MEAN, hfd_6w_IPGTT_glucose_MEAN]  ;
NEWDATA_IPGTT.SEM  = [chow_6w_IPGTT_glucose_SEM, hfd_6w_IPGTT_glucose_SEM]  ;
NEWDATA_IPGTT.time = [ 0 15 30 60 120];

save('NEWDATA_IPGTT','NEWDATA_IPGTT')

% load('Stenkula_IPGTT_6w.mat')
% subplot(2,1,1)
% errorbar(Stenkula_IPGTT_6w.time, Stenkula_IPGTT_6w.ChowMean,Stenkula_IPGTT_6w.ChowSEM, 'k*')
% hold on
% errorbar(NEWDATA_IPGTT.time +1, NEWDATA_IPGTT.Mean(:,1) ,NEWDATA_IPGTT.SEM(:,1), 'g*')
% 
% subplot(2,1,2)
% errorbar(Stenkula_IPGTT_6w.time, Stenkula_IPGTT_6w.HFDMean,Stenkula_IPGTT_6w.HFDSEM, 'k*')
% hold on
% errorbar(NEWDATA_IPGTT.time +1, NEWDATA_IPGTT.Mean(:,2) ,NEWDATA_IPGTT.SEM(:,2), 'g*')