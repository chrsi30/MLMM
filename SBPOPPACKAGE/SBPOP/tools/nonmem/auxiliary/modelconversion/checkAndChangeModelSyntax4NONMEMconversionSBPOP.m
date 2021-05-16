function [ model ] = checkAndChangeModelSyntax4NONMEMconversionSBPOP( model,dosing,FLAG_CMT )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check reactions - for now we do not allow reactions in the SBmodel
% for NONMEM conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(SBreactions(model)) > 0,
    error('Model not allowed to contain reactions.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check dosing for other than BOLUS and INFUSION
% We will only allow BOLUS and INFUSION for NONMEM conversion ... not a
% limitation but simpler and cleaner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds = struct(dosing);
allTypes = {ds.inputs.type};
if ismember('ABSORPTION1',allTypes) || ismember('ABSORPTION0',allTypes),
    error('ABSORPTION0 and ABSORPTION1 dosings not allowed for NONMEM conversion.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if "T" is used in the model 
% and check other things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first check that "T" is not a state, variable, parameter or reaction name
sn = SBstates(model);
pn = SBparameters(model);
vn = SBvariables(model);
rn = SBreactions(model);
allelements = {sn{:} pn{:} vn{:} rn{:}};
if ~isempty(strmatchSB('T',allelements,'exact')),
    error('''T'' defined in the model, but NONMEM uses it as the time variable.');
end
if ~isempty(strmatchSB('F',allelements,'exact')),
    error('''F'' defined in the model, but NONMEM uses it as reserved word - please change it in the model.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPLACE some functions wth fortran equivalents accepted by NONMEM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nonmemfct = {'LOG','LOG10','EXP','SQRT','SIN','COS','ABS','TAN','ASIN','ACOS','ATAN','ABS'};
sbmodelfct = lower({'LOG','LOG10','EXP','SQRT','SIN','COS','ABS','TAN','ASIN','ACOS','ATAN','ABS'});
for k=1:length(nonmemfct),
    model = replaceelementSB(model,sbmodelfct{k},nonmemfct{k});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPLACE ^ by **
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelTextStructure = convertModelToTextSB(model);
modelText = setPartsToCompleteTextSB(modelTextStructure);
modelText = strrep(modelText,'^','**');
[SBstructure,errorMsg] = convertTextToModelSB(modelText);
model = SBmodel(SBstructure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle the FLAG_CMT thingy - but only if FLAG_CMT=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FLAG_CMT==1,
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if multiple inputs on same state - if yes => error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ms = struct(model);
for k=1:length(ms.inputs),
    if length(ms.inputs(k).stateindex)>1,
        error(sprintf('An INPUTn definition is used on more than one state. This can not be handled.'))
    end
end
if length(ms.inputs) ~= length(unique([ms.inputs.stateindex])),
    error(sprintf('Multiple INPUTn definitions on the same state.\n\tThe ADM/YTYPE version can not be used.\n\tPlease use the CMT version.'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now reorder the states according to the INPUTn numbers
% INPUT1 => state 1
% INPUT2 => state 2
% ...
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ms = struct(model);
input_numbers = [];
state_numbers = [];
for k=1:length(ms.inputs),
    input_numbers(end+1) = eval(strrep(ms.inputs(k).name,'INPUT',''));
    state_numbers(end+1) = ms.inputs(k).stateindex;
end
% Reorder states
ms.states(input_numbers) = ms.states(state_numbers);
% Need to reassign stateindex in ms.states.input
for k=1:length(ms.inputs),
    ms.inputs(k).stateindex = input_numbers(k);
end
% Create again a model and return
model = SBmodel(ms);
