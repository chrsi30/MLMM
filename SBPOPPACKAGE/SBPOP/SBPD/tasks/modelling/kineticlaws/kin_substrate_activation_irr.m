function [R] = kin_substrate_activation_irr(V,substrate,Ksa,Ksc)
% kin_substrate_activation_irr: Substrate activation (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_substrate_activation_irr(V,substrate,Ksa,Ksc)
%
% Output Arguments:
% =================
% R = V*(substrate/Ksa)^2 / ( 1+substrate/Ksc+substrate/Ksa+(substrate/Ksa)^2 ) 

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = V*(substrate/Ksa)^2 / ( 1+substrate/Ksc+substrate/Ksa+(substrate/Ksa)^2 );

