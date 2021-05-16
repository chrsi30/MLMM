function [R] = kin_allosteric_inihib_empirical_rev(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki,n)
% kin_allosteric_inihib_empirical_rev: Allosteric inhibition (reversible) kinetics
% 
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_allosteric_inihib_empirical_rev(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki,n)
%
% Output Arguments:
% =================
% R = (Vf*substrate/Kms - Vr*product/Kmp) / (1 + substrate/Kms + product/Kmp + (inhibitor/Ki)^n)

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = (Vf*substrate/Kms - Vr*product/Kmp) / (1 + substrate/Kms + product/Kmp + (inhibitor/Ki)^n);

