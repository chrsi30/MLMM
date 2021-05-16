function [R] = kin_hill_cooperativity_irr(V,substrate,h,Shalve)
% kin_hill_cooperativity_irr: Hill type (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_hill_cooperativity_irr(V,substrate,h,Shalve)
%
% Output Arguments:
% =================
% R = V*substrate^h / ( Shalve^h + substrate^h )

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = V*substrate^h / ( Shalve^h + substrate^h );

