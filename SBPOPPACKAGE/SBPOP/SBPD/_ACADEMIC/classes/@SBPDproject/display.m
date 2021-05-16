function [] = display(project)
% display: Displays information about SBPDproject. This function is 
% called by MATLAB whenever an object is the result of a statement that
% is not terminated by a semicolon. 

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLLECT INFORMATION ABOUT THE PROJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberModels = length(project.models);
numberExperiments = length(project.experiments);
numberMeasurements = length([project.experiments(1:end).measurements]);
numberEstimations = length(project.estimations);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text = sprintf('\tSBPDproject\n\t===========\n');
text = sprintf('%s\tName: %s\n',text,project.name);
text = sprintf('%s\tNumber Models:\t\t\t%d\n',text,numberModels);
text = sprintf('%s\tNumber Experiments:\t\t%d\n',text,numberExperiments);
text = sprintf('%s\tNumber Measurements:\t%d\n',text,numberMeasurements);
text = sprintf('%s\tNumber Estimations:\t\t%d\n',text,numberEstimations);
disp(text);