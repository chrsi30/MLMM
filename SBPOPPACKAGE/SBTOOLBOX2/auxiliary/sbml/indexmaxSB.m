function [result] = indexmaxSB(varargin)
% matchmaxSB: This function searches for the maximum input argument and returns the
% index of the max.
% 
% USAGE:
% ======
% [result] = matchmaxSB(arg1,arg2,...,argn)
%
% arg1...argn: scalar input arguments of type double.
%
% Output Arguments:
% =================
% result: index of the largest element of X

% Information:
% ============
% Copyright (C) 2009  Henning Schmidt, henning@sbtoolbox2.org
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO IT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[v,result] = max([varargin{:}]);
return

