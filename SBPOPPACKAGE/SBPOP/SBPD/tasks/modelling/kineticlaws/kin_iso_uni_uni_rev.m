function [R] = kin_iso_uni_uni_rev(Vf,substrate,product,Keq,Kii,Kms,Kmp)
% kin_iso_uni_uni_rev: enzyme isomerization product inhibition
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_iso_uni_uni_rev(Vf,substrate,product,Keq,Kii,Kms,Kmp)
%
% Output Arguments:
% =================
% R = Vf*(substrate-product/Keq) / ( substrate*(1+product/Kii) + Kms*(1+product/Kmp) )

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = Vf*(substrate-product/Keq) / ( substrate*(1+product/Kii) + Kms*(1+product/Kmp) );

