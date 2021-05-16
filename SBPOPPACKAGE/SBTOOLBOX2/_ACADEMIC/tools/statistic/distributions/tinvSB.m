function inv = tinvSB(x, n)
% Quantile function of the t distribution
%
% USAGE:
% ======
% inv = tinvSB(x, n)
%
% For each component of 'x', compute the quantile (the inverse of
% the CDF) at 'x' of the t (Student) distribution with parameter
% 'n'.

% Information:
% ============
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>
% This function has been taken from Octave and adapted for the SBTOOLBOX2
% by Henning Schmidt
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

  if (nargin ~= 2)
    error('Incorrect number of input arguments'); 
  end

  if (~isscalar (n))
    [retval, x, n] = common_size (x, n);
    if (retval > 0)
      error ('x and n must be of common size or scalar');
    end
  end

  inv = zeros (size (x));

  k = find ((x < 0) | (x > 1) | isnan (x) | ~(n > 0));
  if (any (k))
    inv(k) = NaN;
  end

  k = find ((x == 0) & (n > 0));
  if (any (k))
    inv(k) = -Inf;
  end

  k = find ((x == 1) & (n > 0));
  if (any (k))
    inv(k) = Inf;
  end

  k = find ((x > 0) & (x < 1) & (n > 0) & (n < 10000));
  if (any (k))
    if (isscalar (n))
      inv(k) = (sign (x(k) - 1/2) .* sqrt (n .* (1 ./ betainvSB(2*min (x(k), 1 - x(k)), n/2, 1/2) - 1)));
    else
      inv(k) = (sign (x(k) - 1/2) .* sqrt (n(k) .* (1 ./ betainvSB(2*min (x(k), 1 - x(k)), n(k)/2, 1/2) - 1)));
    end
  end

  % For large n, use the quantiles of the standard normal
  k = find ((x > 0) & (x < 1) & (n >= 10000));
  if (any (k))
    inv(k) = stdnormalinvSB(x(k));
  end

return
