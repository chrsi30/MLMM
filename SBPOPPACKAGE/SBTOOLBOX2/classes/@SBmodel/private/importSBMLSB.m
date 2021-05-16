function [SBstructure,errorMsg] = importSBMLSB(model)
% importSBMLSB: Just an interface to SBimportSBML, which is an ACADEMIC function
%               The dependency check (checkDependenciesSBPOP) then will ignore the 
%               importSBMLSB function.

% Run SBimportSBML
[SBstructure,errorMsg] = SBimportSBML(model);