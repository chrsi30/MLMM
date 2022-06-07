
%%% Hansson 2019 et al
%%PKB
pPKB_S473_PKB_tot_chow = [ 1.059257 0.940743 1.030086 0.969914 ];
pPKB_S473_PKB_tot_chow_mean = mean(pPKB_S473_PKB_tot_chow);
pPKB_S473_PKB_tot_chow_SEM = std(pPKB_S473_PKB_tot_chow)/sqrt(length(pPKB_S473_PKB_tot_chow));

pPKB_S473_PKB_tot_hfd= [ 0.471271 0.460685 0.546091 0.448691 0.413751 0.836023 ];
pPKB_S473_PKB_tot_hfd_mean = mean(pPKB_S473_PKB_tot_hfd);
pPKB_S473_PKB_tot_hfd_SEM = std(pPKB_S473_PKB_tot_hfd)/sqrt(length(pPKB_S473_PKB_tot_hfd));


Hansson2019_invitro.PKBS473.Mean = [ pPKB_S473_PKB_tot_chow_mean , pPKB_S473_PKB_tot_hfd_mean ];
Hansson2019_invitro.PKBS473.SEM = [ pPKB_S473_PKB_tot_chow_SEM , pPKB_S473_PKB_tot_hfd_SEM ];

%% AS160
pAS160_T612_AS160tot_chow = [ 1.07774 0.92226 1.210462 0.789538 ];
pAS160_T612_AS160tot_chow_mean = mean( pAS160_T612_AS160tot_chow );
pAS160_T612_AS160tot_chow_SEM = std(pAS160_T612_AS160tot_chow)/sqrt(length(pAS160_T612_AS160tot_chow));

pAS160_T612_AS160tot_hfd = [ 0.406311 0.449394 0.333613 0.069322 0.189767 0.269715 ];
pAS160_T612_AS160tot_hfd_mean = mean(pAS160_T612_AS160tot_hfd);
pAS160_T612_AS160tot_hfd_SEM = std(pAS160_T612_AS160tot_hfd)/sqrt(length(pAS160_T612_AS160tot_hfd));

Hansson2019_invitro.AS160.Mean = [ pAS160_T612_AS160tot_chow_mean , pAS160_T612_AS160tot_hfd_mean ];
Hansson2019_invitro.AS160.SEM = [ pAS160_T612_AS160tot_chow_SEM , pAS160_T612_AS160tot_hfd_SEM ];

%% 14D Glucose uptake - Figure 6b)
% Epididymal (left graph) and Inguinal (right graph) adipocytes isolated after chow, HFD or reversed feeding were subjected to tracer glucose uptake assay,
% either non-stimulated (basal), or insulin-stimulated (0.1 nM). FOR 30 min
% (?)
% ? USING ONLY insulin stimulated Inguinal adipocyte data 
% FOLD , Digitized data

Hansson2019_invitro.GU.Mean = [ 1413  322 ]./1413;
Hansson2019_invitro.GU.SEM  = [ 173    67 ]./1413;





save('Hansson2019_invitro','Hansson2019_invitro')