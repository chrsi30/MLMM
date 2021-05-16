function [R] = kin_michaelis_menten_rev(Vf,substrate,Kms,Vr,product,Kmp)
% kin_michaelis_menten_rev: Michaelis Menten (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_michaelis_menten_rev(Vf,substrate,Kms,Vr,product,Kmp)
%
% Output Arguments:
% =================
% R = (Vf*substrate/Kms+Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp ) 

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = (Vf*substrate/Kms+Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp ) ;

