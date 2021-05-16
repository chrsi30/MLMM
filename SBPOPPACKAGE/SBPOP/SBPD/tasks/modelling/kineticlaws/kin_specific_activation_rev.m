function [R] = kin_specific_activation_rev(Vf,substrate,Kms,Vr,product,Kmp,activator,Ka)
% kin_specific_activation_rev: Specific activation (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_specific_activation_rev(Vf,substrate,Kms,Vr,product,Kmp,activator,Ka)
%
% Output Arguments:
% =================
% R = ( Vf*substrate/Kms-Vr*product/Kmp )*activator / ( Ka+(1+substrate/Kms+product/Kmp)*activator )  

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = ( Vf*substrate/Kms-Vr*product/Kmp )*activator / ( Ka+(1+substrate/Kms+product/Kmp)*activator );

