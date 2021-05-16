function [R] = kin_substrate_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,Ki)
% kin_substrate_inihib_rev: Substrate inhibition (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_substrate_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,Ki)
%
% Output Arguments:
% =================
% R = (Vf*substrate/Kms-Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp+(substrate/Ki)^2 ) 

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = (Vf*substrate/Kms-Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp+(substrate/Ki)^2 );

