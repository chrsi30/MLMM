function SBplotselected(simdata,plotcomponents,varargin)
% SBplotselected: Takes a simulation result as input. Additionally, the
% names of states, variables, and reactions can be defined for which the
% values should be plotted. 
%
% USAGE:
% ======
% SBplotselected(simdata,plotcomponents)
% SBplotselected(simdata,plotcomponents,headers)
%
% simdata: Simulation results returned, e.g. from SBsimulate
% plotcomponents: cell-array with component names to plot (states,
%   variables and/or reactions). Alternatively, the cell-array can contain 
%   only cell-arrays of component names. This allows to define
%   plot-subgroups, which can be selected using the pulldown menu in the
%   upper left corner of the SBplot window.
% headers: if plotcomponents contains only cell-arrays, then this input
%   argument contains the names to be displayed in the pull-down menu.
%
% DEFAULT VALUES:
% ===============
% headers: {'plot 1', 'plot 2', ...}

% Information:
% ============
% Copyright (C) 2008 Henning Schmidt, henning@sbtoolbox2.org
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
% USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check plotcomponents argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(plotcomponents),
    error('SBT2:SBplotselected:componentsNotCellArray','Second input argument needs to be a cell-array.');
end
charfound = 0;
cellfound = 0;
for k=1:length(plotcomponents),
    if ischar(plotcomponents{k}),
        charfound = 1;
    elseif iscell(plotcomponents{k}),
        cellfound = 1;
    end
end
if charfound == 1 && cellfound == 1,
    error('SBT2:SBplotselected:charCellMix','Second input argument needs to contain either chars or cells. Not both!');
end
% convert char to cell model
if charfound == 1,
    plotcomponents = {plotcomponents};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
headers = {};
for k=1:length(plotcomponents),
    headers{k} = sprintf('plot %d',k);
end
if nargin == 3,
    headers = varargin{1};
end
if length(headers) ~= length(plotcomponents),
    error('SBT2:SBplotselected:wrongHeadersLength','Wrong number of elements in third input argument.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the different plot structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotstructures = {};
states = simdata.states; sind = 1:length(states);
variables = simdata.variables; vind = 1:length(variables);
reactions = simdata.reactions; rind = 1:length(reactions);
for k=1:length(plotcomponents),
    comp = plotcomponents{k};
    simdatak = simdata;
    sind_keep = [];
    vind_keep = [];
    rind_keep = [];
    for k2=1:length(comp),
        % reduce simulation results to only contain the components in comp
        sind_keep = [sind_keep strmatchSB(comp{k2},states,'exact')];
        vind_keep = [vind_keep strmatchSB(comp{k2},variables,'exact')];
        rind_keep = [rind_keep strmatchSB(comp{k2},reactions,'exact')];
    end
    simdatak.states = simdatak.states(sind_keep);
    simdatak.statevalues = simdatak.statevalues(:,sind_keep);
    simdatak.variables = simdatak.variables(vind_keep);
    simdatak.variablevalues = simdatak.variablevalues(:,vind_keep);
    simdatak.reactions = simdatak.reactions(rind_keep);
    simdatak.reactionvalues = simdatak.reactionvalues(:,rind_keep);
    % convert to plot structure 
    plotstructures{k} = createdatastruct2SBplotSB(simdatak,headers{k}); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct SBplot command and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotcommand = 'SBplot(';
for k=1:length(plotstructures),
    plotcommand = sprintf('%splotstructures{%d},',plotcommand,k);
end
plotcommand = [plotcommand(1:end-1) ')'];
eval(plotcommand)
