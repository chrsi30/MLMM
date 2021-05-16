function [R] = kin_catalytic_activation_irr(V,substrate,activator,Kms,Ka)
% kin_catalytic_activation_irr: Catalytic activation (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_catalytic_activation_irr(V,substrate,activator,Kms,Ka)
%
% Output Arguments:
% =================
% R = V*substrate*activator / ( (Kms + substrate)*(Ka+activator) )

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = V*substrate*activator / ( (Kms + substrate)*(Ka+activator) )

