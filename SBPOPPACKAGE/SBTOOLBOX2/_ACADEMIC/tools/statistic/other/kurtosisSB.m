function retval = kurtosisSB(x, dim)
% If 'x' is a vector of length 'N', return the kurtosis
% of 'x'.  If 'x' is a matrix, return the kurtosis over the
% first non-singleton dimension. The optional argument 'dim'
% can be given to force the kurtosis to be given over that 
% dimension.
%
% USAGE:
% ======
% retval = kurtosisSB(x, dim)

% Information:
% ============
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>
% Adapted-By: jwe
% Adapted for the SBTOOLBOX2 by Henning Schmidt
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

  if (nargin ~= 1 && nargin ~= 2)
    error('Incorrect number of input arguments'); 
  end

  nd = ndims (x);
  sz = size (x);
  if (nargin ~= 2)
    % Find the first non-singleton dimension.
    dim  = 1;
    while (dim < nd + 1 && sz(dim) == 1)
      dim = dim + 1;
    end
    if (dim > nd)
      dim = 1;
    end
  else
    if (~(isscalar (dim) && dim == round (dim)) && dim > 0	&& dim < (nd + 1))
      error ('dim must be an integer and valid dimension');
    end
  end
  
  if (~ismatrixSB(x) && ~isvector(x))
    error ('x has to be a matrix or a vector');
  end

  c = sz(dim);
  sz(dim) = 1;
  idx = ones (1, nd);
  idx(dim) = c;
  x = x - repmat (mean (x, dim), idx);
  retval = zeros (sz);
  s = std (x, [], dim);
  x = sum(x.^4, dim);
  ind = find (s > 0);
  retval(ind) = x(ind) / (c*s(ind) .^ 4);

return
