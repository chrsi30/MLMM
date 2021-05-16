function [R] = kin_substrate_inihib_irr(V,substrate,Km,Ki)
% kin_substrate_inihib_irr: Substrate inhibition (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_substrate_inihib_irr(V,substrate,Km,Ki)
%
% Output Arguments:
% =================
% R = V*substrate / ( Km + substrate + Km*(substrate/Ki)^2 )

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = V*substrate / ( Km + substrate + Km*(substrate/Ki)^2 );

