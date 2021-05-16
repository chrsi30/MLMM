function [R] = kin_michaelis_menten_irr(V,substrate,Km)
% kin_michaelis_menten_irr: Michaelis Menten (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_michaelis_menten_irr(V,substrate,Km)
%
% Output Arguments:
% =================
% R = V*substrate / ( Km + substrate )

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = V*substrate / ( Km + substrate );

