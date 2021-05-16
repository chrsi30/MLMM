function [output] = maplearraySB(input)
% maplearraySB: Creates a symbolic expression of the form: '{a,b,c,d}' 
% where a, b,c ,d are strings coming from a cellarray input argument.

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

academicWarningSB


for kindex=1:length(input)
    if kindex==1
        output = sprintf('{%s', input{kindex});
    else
        output = sprintf('%s,%s', output, input{kindex});
    end
end
output = sym(strcat(output,'}'));
return