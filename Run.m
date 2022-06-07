

clear all
clc
close all

format long
format compact

warning 'off'

modelName = 'MLMM_final_6'; 

addpath('./Models'); 
addpath('./Data'); 

addpath('./CODE/Simulations')
addpath('./CODE/Plot'); 
addpath('./CODE/Common Code'); 


D = Setup_data;

objModel = SBmodel(strcat(modelName,'.txt')); 
[pNames, p0] = SBparameters(objModel); ic0 = SBinitialconditions(objModel); func_mex_model = str2func(modelName);
utility.pNames = pNames;
utility.objModel = objModel;
utility.ic0=ic0;
utility.p0 = p0;
utility.modelName = modelName;

filedir = strcat('./Parameters/minmax');
files=dir(fullfile(filedir, '*.mat'));
nfiles = length(files);

msg = strcat("What figure to plot?;" );
tmp = ["Figure 3 "...
       "Figure 4 "...
       "Figure 5 "...
       "Figure 6 "...
       "Figure 7"];
 
choice = menu(msg,tmp);

if choice == 1
    Figure_3  
elseif choice == 2
    Figure_4     
elseif choice == 3
    Figure_5    
elseif choice == 4
    Figure_6    
elseif choice == 5    
    %% Setup
    
    modelName = 'MLMM_extended_v0_1';
    
    cd('./Models')
    SBPDmakeMEXmodel(SBmodel(strcat(modelName,'.txt')))
    cd ..
    
    objModel = SBmodel(strcat(modelName,'.txt'));
    [pNames, p0] = SBparameters(objModel); ic0 = SBinitialconditions(objModel); func_mex_model = str2func(modelName);
    utility.pNames = pNames;
    utility.objModel = objModel;
    utility.ic0=ic0;
    utility.p0 = p0;
    utility.modelName = modelName;
    
    
    ic0       = utility.ic0;
    objModel  = utility.objModel; % variables used for indexing
    modelName = utility.modelName;
    pNames    = utility.pNames ;
    p         = utility.p0; % Parmeter vector with values from model files. - Vector used for simulation
    
    
    Figure_7

end



