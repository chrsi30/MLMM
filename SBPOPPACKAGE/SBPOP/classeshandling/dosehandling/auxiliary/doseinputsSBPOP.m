function [names,types] = doseinputsSBPOP(dos)
% doseinputsSBPOP: Extracts names and types of dosing inputs from a
% SBPOPdosing object.
%
% USAGE:
% ======
% [names,types] = doseinputsSBPOP(dos) 
%
% dos: SBPOPdosing object
%
% Output Arguments:
% =================
% names: cell-array with input names
% types: cell-array with input types

% Information:
% ============
% Copyright © 2012 Novartis Pharma AG
% 
% This program is Free Open Source Software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 
% 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License 
% along with this program. If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isSBPOPdosing(dos),
    error('Input argument is not an SBPOPdosing object.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get names and types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds = struct(dos);
names = {ds.inputs.name};
types = {ds.inputs.type};
