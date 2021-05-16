function pdf = tpdfSB(x, n)
% Probability density function of the t distribution
%
% USAGE:
% ======
% pdf = tpdfSB(x, n)
%
% For each element of 'x', compute the probability density function
% (PDF) at 'x' of the t (Student) distribution with 'n'
% degrees of freedom. 

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

  pdf = zeros (size (x));

  k = find (isnan (x) | ~(n > 0) | ~(n < Inf));
  if (any (k))
    pdf(k) = NaN;
  end

  k = find (~isinf (x) & ~isnan (x) & (n > 0) & (n < Inf));
  if (any (k))
    if (isscalar (n))
      pdf(k) = (exp (- (n + 1) .* log (1 + x(k) .^ 2 ./ n)/2) / (sqrt (n) * beta(n/2, 1/2)));
    else
      pdf(k) = (exp (- (n(k) + 1) .* log (1 + x(k) .^ 2 ./ n(k))/2) ./ (sqrt (n(k)) .* beta(n(k)/2, 1/2)));
    end
  end

return
