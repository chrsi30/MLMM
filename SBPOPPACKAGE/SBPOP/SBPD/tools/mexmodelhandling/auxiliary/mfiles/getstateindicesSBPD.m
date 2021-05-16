function [output] = getstateindicesSBPD(model,statenames)
% getstateindicesSBPD: Determines a vector of model state indices, the 
% order corresponding to the states argument.
% 
% USAGE:
% ======
% [output] = getstateindicesSBPD(model,states)
% 
% model: SBmodel, ODE file model, or MEX simulation model
% states: Cell-array of state names
% 
% OUTPUT:
% =======
% output: Vector of state indices

% Information:
% ============
% Copyright (C) 2005-2013 Henning Schmidt, henning@sbtoolbox2.org
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

% Convert to cell array if necessary
if ischar(statenames),
    statenames = {statenames};
end
% Get all state names in order as stored in the model
allnames = SBstates(model);
% Get the indices
output = zeros(1,length(statenames));
for k = 1:length(statenames),
    index = strmatchSB(statenames{k},allnames,'exact');
    if isempty(index),
        error(sprintf('State ''%s'' does not exist in the model.',statenames{k}));
    end
    output(k) = index;
end
return