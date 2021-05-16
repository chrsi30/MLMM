function [] = SBPDplotmeasurements(project,varargin)
% SBPDplotmeasurements: Plots all measurements in the given SBPDproject
% using SBplot. Very useful to get a quick overview over measurement
% results. The function opens a simple GUI where you in the upper left
% corner can select the experiment and measurement to display. If error
% bound information is available in the measurement data this is displayed.
%
% USAGE:
% ======
% [] = SBPDplotmeasurements(project)        
% [] = SBPDplotmeasurements(project,experimentindices)        
%
% project:  SBPDproject object
% experimentindices: vector with indices of the experiments for which to
%   plot the measurement data. Per default the measurement data of all
%   experiments is plotted.

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org
academicWarningSBPD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isSBPDproject(project),
    error('Input argument is not an SBPDproject.');
end
project = SBPDstruct(project);
experiments = project.experiments;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle variable input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experimentindices = [1:length(experiments)];
if nargin == 2,
    experimentindices = varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get all measurements from all experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datastructures = {};
for k=1:length(experimentindices),
    es = SBstruct(experiments(experimentindices(k)).experiment);
    for k2=1:length(experiments(experimentindices(k)).measurements),
        % get SBplot datastructure for measurement
        ds = SBvisualizemeasurement(experiments(experimentindices(k)).measurements{k2});
        % update name of measurement in plotstructure
        ds.name = sprintf('E%d: %s, M%d: %s',experimentindices(k),experiments(experimentindices(k)).name,k2,ds.name);
        datastructures{end+1} = ds;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error if empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(datastructures),
    error('The project does not contain any measurements.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display using SBplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct call for SBplot
plotcall = 'SBplot(';
for k=1:length(datastructures),
    plotcall = sprintf('%sdatastructures{%d},',plotcall,k);
end
plotcall = [plotcall(1:end-1) ');'];
eval(plotcall);
