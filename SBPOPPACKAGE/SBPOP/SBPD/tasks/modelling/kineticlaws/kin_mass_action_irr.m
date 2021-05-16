function [R] = kin_mass_action_irr(k,substrate)
% kin_mass_action_irr: Mass action (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate (several substrates can be realized by
%   substrate=substrate1*substrate2*...
%
% USAGE:
% ======
% R = kin_mass_action_irr(k,substrate)
%
% Output Arguments:
% =================
% R = k*substrate

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = k*substrate;

