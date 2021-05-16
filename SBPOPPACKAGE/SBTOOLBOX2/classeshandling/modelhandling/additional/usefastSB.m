function [output] = usefastSB(model)
% usefastSB: checks if an SBmodel contains "fast" reactions.
%
% Output Arguments:
% =================
% output: =1 if reactions with the "fast" flag set are present, =0 if not

% Information:
% ============
% Copyright (C) 2008  Henning Schmidt, henning@sbtoolbox2.org
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

if ~isSBmodel(model),
    error('Given model is not an SBmodel.');
end

ms = struct(model);
output = ~isempty(find([ms.reactions.fast]==1, 1));

