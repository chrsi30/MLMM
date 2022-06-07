function[ out, cost ]= cost_function_v6_pred( theta ,func_mex_model, D, utility, pol, e,x )

if nargin <6, e='IPGTT_hfd6w_gcc'; end %?
if nargin <7, x=120/(60*24)+56+0.5833; end

try   
    
[ AS ] = f_as_new( theta ,func_mex_model, D, utility);

ss                   = AS.ss;
chow8w               = AS.chow8w;
hfd8w                = AS.hfd8w;
IPGTT_hfd6w          = AS.IPGTT_hfd6w;
chow2w               = AS.chow2w;
IPGTT_chow2w         = AS.IPGTT_chow2w;
IPGTT_chow12dHFD2d   = AS.IPGTT_chow12dHFD2d;
chow10dHFD4d         = AS.chow10dHFD4d;
chow8dHFD6d          = AS.chow8dHFD6d;
IPGTT_chow8dHFD6d    = AS.IPGTT_chow8dHFD6d;
hfd2w                = AS.hfd2w;
IPGTT_hfd2w          = AS.IPGTT_hfd2w;


%
id_BW    =  ismember(SBvariables(utility.objModel),'BW');
id_Iconc =  ismember(SBvariables(utility.objModel),'Iconc');
id_Gconc =  ismember(SBvariables(utility.objModel),'Gconc');

%% Steady state 

cost_ss_gcc = (( ss.variablevalues(end, id_Gconc )  -  D.Hansson2018.IPGTT.Glucose.Mean(1,1)  )^2 ) ./( D.Hansson2018.IPGTT.Glucose.SEM(1,1) )^2 ;% basal glucose
cost_ss_icc = (( ss.variablevalues(end, id_Iconc )   -  D.Hansson2018.IPGTT.Insulin.Mean(1,1)  )^2 )./(  D.Hansson2018.IPGTT.Insulin.SEM(1,1) )^2 ;% basal insulin

ss_cost = cost_ss_gcc + cost_ss_icc;
    
%% COST -  Hansson et al (2018) 2w 

%%% Body weight
bw_cost1 = nansum (  ( (chow2w.variablevalues([1,end],  id_BW) -  D.Hansson2018.BW.Mean(:,1) ).^2 ) ./( D.Hansson2018.BW.SEM(:,1).^2 ));    

bw_cost3 = nansum (  ( (chow10dHFD4d.variablevalues(end,id_BW) -  D.Hansson2018.BW.Mean(end,2)).^2 )./(  D.Hansson2018.BW.SEM(end,2).^2 )); 

bw_cost4 = nansum (  ( (chow8dHFD6d.variablevalues(end, id_BW) -  D.Hansson2018.BW.Mean(end,3)).^2 )./( D.Hansson2018.BW.SEM(end ,3).^2 )); 

bw_cost5 = nansum (  ( (hfd2w.variablevalues([1,end],   id_BW) -  D.Hansson2018.BW.Mean(:,4) ).^2 )./(  D.Hansson2018.BW.SEM(:,4).^2 ));    

cost_bw = bw_cost1 + bw_cost3  + bw_cost4 + bw_cost5 ;

%%% IPGTT cost 
id_tp_IPGTT = [ 1 16 31 61 121 ];

IPGTT_glucose_cost1 =  nansum (  ( (IPGTT_chow2w.variablevalues(id_tp_IPGTT, id_Gconc)         - D.Hansson2018.IPGTT.Glucose.Mean(:,1)     ).^2 )./(  D.Hansson2018.IPGTT.Glucose.SEM(:,1) .^2 ));  
IPGTT_glucose_cost2 =  nansum (  ( (IPGTT_chow12dHFD2d.variablevalues(id_tp_IPGTT, id_Gconc)   - D.Hansson2018.IPGTT.Glucose.Mean(:,2)     ).^2 )./(  D.Hansson2018.IPGTT.Glucose.SEM(:,2) .^2 ));  
IPGTT_glucose_cost3 =  nansum (  ( (IPGTT_chow8dHFD6d.variablevalues(id_tp_IPGTT, id_Gconc)    - D.Hansson2018.IPGTT.Glucose.Mean(:,3)     ).^2 )./(  D.Hansson2018.IPGTT.Glucose.SEM(:,3) .^2 ));  
IPGTT_glucose_cost4 =  nansum (  ( (IPGTT_hfd2w.variablevalues(id_tp_IPGTT, id_Gconc)          - D.Hansson2018.IPGTT.Glucose.Mean(:,4)     ).^2 )./(  D.Hansson2018.IPGTT.Glucose.SEM(:,4) .^2 ));  

tp_icc = 1; % [ 1 16 121 ]; % only first data point

IPGTT_insulin_cost1 =  nansum (  ( (IPGTT_chow2w.variablevalues(tp_icc, id_Iconc)         - D.Hansson2018.IPGTT.Insulin.Mean(1,1)     ).^2 )./(  D.Hansson2018.IPGTT.Insulin.SEM(1,1) .^2 ));  
IPGTT_insulin_cost2 =  nansum (  ( (IPGTT_chow12dHFD2d.variablevalues(tp_icc, id_Iconc)   - D.Hansson2018.IPGTT.Insulin.Mean(1,2)     ).^2 )./(  D.Hansson2018.IPGTT.Insulin.SEM(1,2) .^2 ));  
IPGTT_insulin_cost3 =  nansum (  ( (IPGTT_chow8dHFD6d.variablevalues(tp_icc, id_Iconc)    - D.Hansson2018.IPGTT.Insulin.Mean(1,3)     ).^2 )./(  D.Hansson2018.IPGTT.Insulin.SEM(1,3) .^2 ));  
IPGTT_insulin_cost4 =  nansum (  ( (IPGTT_hfd2w.variablevalues(tp_icc, id_Iconc)          - D.Hansson2018.IPGTT.Insulin.Mean(1,4)     ).^2 )./(  D.Hansson2018.IPGTT.Insulin.SEM(1,4) .^2 ));  

cost_IPGTT =  (IPGTT_glucose_cost1 + IPGTT_insulin_cost1) + (IPGTT_glucose_cost2 + IPGTT_insulin_cost2) + (IPGTT_glucose_cost3 + IPGTT_insulin_cost3)  + (IPGTT_glucose_cost4 + IPGTT_insulin_cost4);
              
%% NEW DATA - 8w chow/hfd -  Stenkula et al 
% Training on the first 2 weeks
id_tp_bw_8w = [ 1 8 15];

bw_cost6 = nansum (  ( (chow8w.variablevalues(id_tp_bw_8w, id_BW)    -  D.NEW.BW.Mean(1:3,1) ).^2 )./( D.NEW.BW.SEM(1:3,1).^2 ));
bw_cost7 = nansum (  ( (hfd8w.variablevalues(id_tp_bw_8w,  id_BW)    -  D.NEW.BW.Mean(1:3,2) ).^2 )./( D.NEW.BW.SEM(1:3,2).^2 )); 

bw_cost_new = bw_cost6 + bw_cost7;


%% Adhoc 
limit = chi2inv(0.95, (numel(D.Hansson2018.BW.Mean)-4 + numel(D.Hansson2018.IPGTT.Glucose.Mean) + 4 + 2 + 6));

adhoc_cost = 0;

% U_A
IPGTT_chow2w_UA = trapz(IPGTT_chow2w.variablevalues(:, ismember(SBvariables(utility.objModel),'U_id_A')) ); 
IPGTT_chow2w_UM = trapz(IPGTT_chow2w.variablevalues(:, ismember(SBvariables(utility.objModel),'U_id_M')) ); 

adhoc_cost = adhoc_cost + Constraint(IPGTT_chow2w_UA, IPGTT_chow2w_UM*0.5, limit);
adhoc_cost = adhoc_cost + Constraint(IPGTT_chow2w_UM*0.15, IPGTT_chow2w_UA, limit);

% U_I  
IPGTT_chow2w_U_id = trapz(IPGTT_chow2w.variablevalues(:, ismember(SBvariables(utility.objModel),'U_id')) ); 
IPGTT_chow2w_U_ii = trapz(IPGTT_chow2w.variablevalues(:, ismember(SBvariables(utility.objModel),'U_ii')) ); 
adhoc_cost = adhoc_cost + Constraint(IPGTT_chow2w_U_ii, IPGTT_chow2w_U_id, limit);
adhoc_cost = adhoc_cost + Constraint(IPGTT_chow2w_U_id*0.05, IPGTT_chow2w_U_ii, limit);

% Insulin  

max_insulin_level = max(D.Hansson2018.IPGTT.Insulin.Mean(:,4))  + 100;

IPGTT_chow2w_max_insulin= max(IPGTT_hfd2w.variablevalues(:, ismember(SBvariables(utility.objModel),'Iconc')) ); 

adhoc_cost_insulin = adhoc_cost + Constraint(IPGTT_chow2w_max_insulin, max_insulin_level, limit);

%% in vitro - primary adipcoytes - insulin signaling


[ AS2 ] = f_as_invitro_new( theta ,func_mex_model,utility);
        
vitro_chow4w = AS2.vitro_chow4w;       
vitro_hfd4w  = AS2.vitro_hfd4w;         

id_PKB_phos473 = ismember(SBvariables(utility.objModel),'totalPKB473'); 
id_PKB         = ismember(SBvariables(utility.objModel),'totalPKB'); 
id_AS          = ismember(SBvariables(utility.objModel),'totalAS160'); 
id_As_phos     = ismember(SBstates(utility.objModel),'AS160p'); 

PKB_fold_chow = vitro_chow4w.variablevalues(end,id_PKB_phos473)/vitro_chow4w.variablevalues(end,id_PKB);
PKB_fold_hfd  = vitro_hfd4w.variablevalues(end,id_PKB_phos473)/vitro_hfd4w.variablevalues(end,id_PKB);
fold_diet_PKB = PKB_fold_hfd/PKB_fold_chow;% fold between diets ( fold on fold )

AS160_fold_chow = vitro_chow4w.statevalues(end,id_As_phos)/vitro_chow4w.variablevalues(end,id_AS);
AS160_fold_hfd  = vitro_hfd4w.statevalues(end,id_As_phos)/vitro_hfd4w.variablevalues(end,id_AS);
fold_diet_AS160 = AS160_fold_hfd/AS160_fold_chow;% fold between diets ( fold on fold )

cost_PKB   =   ( (  fold_diet_PKB   - D.Hansson2019.invitro.PKBS473.Mean(1,2) ).^2 )./( D.Hansson2019.invitro.PKBS473.SEM(1,2) .^2 );    
cost_AS160 =   ( (  fold_diet_AS160 - D.Hansson2019.invitro.AS160.Mean(1,2)   ).^2 )./( D.Hansson2019.invitro.AS160.SEM(1,2)  .^2 );  

cost_in_vitro= cost_PKB + cost_AS160;

%%% 14D- GLUCOSE UPTAKE
id_GU          = ismember(SBvariables(utility.objModel),'CL_GIA'); 

AUC_chow4w_GU     = trapz(vitro_chow4w.variablevalues(:,id_GU));
AUC_hfd4w_GU      = trapz( vitro_hfd4w.variablevalues(:,id_GU));
fold_AUC_hfd4w_GU = AUC_hfd4w_GU/AUC_chow4w_GU; % fold on chow

cost_GU = (( fold_AUC_hfd4w_GU - D.Hansson2019.invitro.GU.Mean(1,2)).^2 )./( D.Hansson2019.invitro.GU.SEM(1,2)  .^2 );  



%% Total cost
cost= cost_bw +  cost_IPGTT + cost_in_vitro + bw_cost_new  + ss_cost  + adhoc_cost + adhoc_cost_insulin + cost_GU;






%% What to optimize
% out = pred_cost;
% out1 = IPGTT_hfd6w.variablevalues(end, end-1);       %% SIM - sista punkten,121
if  contains(e, 'PKB') 
    out = fold_diet_PKB;
elseif contains(e, 'AS160')
    out = fold_diet_AS160;
elseif contains(e, 'cost')

    id_tp_bw_8w_2 = [ 22 29 36 42 51 57];

    bw_cost8 = nansum (  ( (chow8w.variablevalues(id_tp_bw_8w_2, id_BW)    -  D.NEW.BW.Mean(4:end,1) ).^2 )./( D.NEW.BW.SEM(4:end,1).^2 ));
    bw_cost9 = nansum (  ( (hfd8w.variablevalues(id_tp_bw_8w_2,  id_BW)    -  D.NEW.BW.Mean(4:end,2) ).^2 )./( D.NEW.BW.SEM(4:end,2).^2 )); 
    BW_prediction_cost = bw_cost8 + bw_cost9;

    id_tp_IPGTT = [ 1 16 31 61 121 ];
    Glu_prediction_cost = nansum (  ( (IPGTT_hfd6w.variablevalues(id_tp_IPGTT, id_Gconc)    -  D.NEW.IPGTT.Mean(:,2) ).^2 )./( D.NEW.IPGTT.SEM(:,2).^2 ));

    pred_cost = BW_prediction_cost + Glu_prediction_cost; % <--- USE THIS
    out = pred_cost;
else
    sim=[];
    switch  e
        case 'chow2w_bw',              sim = chow2w;             ind = 1;
        case 'chow10dHFD4d_bw',        sim = chow10dHFD4d;       ind = 1;
        case 'chow8dHFD6d_bw',         sim = chow8dHFD6d;        ind = 1;
        case 'hfd2w_bw',               sim = hfd2w;              ind = 1;
        case 'chow8w_bw',              sim = chow8w;             ind = 1;
        case 'hfd8w_bw',               sim = hfd8w;              ind = 1;
        case 'IPGTT_chow2w_gcc',       sim = IPGTT_chow2w;       ind = size(sim.variablevalues,2)-1;
        case 'IPGTT_chow12dHFD2d_gcc', sim = IPGTT_chow12dHFD2d; ind = size(sim.variablevalues,2)-1;
        case 'IPGTT_chow8dHFD6d_gcc',  sim = IPGTT_chow8dHFD6d;  ind = size(sim.variablevalues,2)-1;
        case 'IPGTT_hfd2w_gcc',        sim = IPGTT_hfd2w;        ind = size(sim.variablevalues,2)-1;
        case 'IPGTT_chow2w_icc',       sim = IPGTT_chow2w;       ind = size(sim.variablevalues,2);
        case 'IPGTT_chow12dHFD2d_icc', sim = IPGTT_chow12dHFD2d; ind = size(sim.variablevalues,2);
        case 'IPGTT_chow8dHFD6d_icc',  sim = IPGTT_chow8dHFD6d;  ind = size(sim.variablevalues,2);
        case 'IPGTT_hfd2w_icc',        sim = IPGTT_hfd2w;        ind = size(sim.variablevalues,2);
        case 'IPGTT_hfd6w_gcc',        sim = IPGTT_hfd6w;        ind = size(sim.variablevalues,2)-1;
        otherwise
            disp(e)
            error('Not a valid experiment')
    end
    out = sim.variablevalues(ismember(sim.time,x),ind);
end

if isempty(out)
    disp(e)
    disp(x)
    error('Empty output! Time point missmatch?')
end
violation = max(cost - limit,0); %Got he limit from below. Intial parameter set was not valid under these conditions.
% p=-1; % -1 to get maximization, +1 to get minimization (max v = min -v)
tol=0;
if violation > tol
    out=pol*out+abs(out)+violation+limit;
else
    out=pol*out;
end



catch err
    %disp(getReport(err))
   cost =1e99;
   out=1e99;
end


end


function [penalty] = Constraint(lower, upper, offset)
viol =  lower - upper  ; 
if viol > 0
    penalty = viol + offset; 
else
    penalty = 0;
end
end



