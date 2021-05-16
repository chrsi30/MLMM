function [elements] = explodePCSB(text,varargin)
% explodePCSB: This function does not(!!!) lead to an explosion of your
% Personal Computer. It is an auxiliary function allowing to decompose a
% string expression into its comma separated elements. Commas within
% parentheses expressions are not considered. This is a useful function to
% obtain the arguments of a function call, where some arguments might have
% expressions involving commas but have parentheses around them.
% Alternatively, a separator character, other then a comma, can be
% specified by the user.
%
% USAGE:
% ======
% [elements] = explodePCSB(text)
% [elements] = explodePCSB(text,separatorCharacter)
% [elements] = explodePCSB(text,separatorCharacter,groupCharacterStart,groupCharacterEnd)
%
% text: text to decompose into its comma separated elements
% separatorCharacter: one character that should be used for the explosion
% of the string.
%
% DEFAULT VALUES:
% ===============
% separatorCharacter: ',' (comma)
% groupCharacterStart: '(' can also be a cell-array with several
%                      parenthesis types
% groupCharacterEnd: ')'  can also be a cell-array with several
%                      parenthesis types
%
% Output Arguments:
% =================
% elements: cell-array containing string elements, previously separated by
% the separator character.

% Information:
% ============
% Copyright (C) 2005-2007 Fraunhofer Chalmers Centre, Gothenburg, Sweden
% Main author: Henning Schmidt
% 
% Changes for the SBTOOLBOX2:
% 1/1/2008  Henning Schmidt, henning@sbtoolbox2.org
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
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
separatorCharacter = ',';
groupCharacterStart = {'('};
groupCharacterEnd = {')'};
if nargin == 1,
elseif nargin == 2,
    separatorCharacter = varargin{1};
elseif nargin == 4,
    separatorCharacter = varargin{1};
    groupCharacterStart = varargin{2};
    groupCharacterEnd = varargin{3};
else
    error('Incorrect number of input arguments.');
end

if ischar(groupCharacterStart),
    groupCharacterStart = {groupCharacterStart};
end
if ischar(groupCharacterEnd),
    groupCharacterEnd = {groupCharacterEnd};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO THE EXPLOSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elements = {};
openParenthesis = 0;
lastIndex = 1;
elementIndex = 1;
doubletext = double(text);
doublegroupCharacterStart = double(([groupCharacterStart{:}]));
doublegroupCharacterEnd = double(([groupCharacterEnd{:}]));
doubleseparatorCharacter = double(separatorCharacter);
for k2 = 1:length(text),
    if sum(doubletext(k2) == doublegroupCharacterStart),
        openParenthesis = openParenthesis + 1;
    elseif sum(doubletext(k2) == doublegroupCharacterEnd),
        openParenthesis = openParenthesis - 1;
    elseif (doubletext(k2) == doubleseparatorCharacter) && (openParenthesis == 0),
        elements{elementIndex} = strtrim(text(lastIndex:k2-1));
        elementIndex = elementIndex + 1;
        lastIndex = k2+1;
    end
end
elements{elementIndex} = strtrim(text(lastIndex:end));
return

