function [R] = kin_allosteric_inihib_mwc_irr(V,substrate,Ks,n,L,inhibitor,Ki)
% kin_allosteric_inihib_mwc_irr: Allosteric inhibition (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_allosteric_inihib_mwc_irr(V,substrate,Ks,n,L,inhibitor,Ki)
%
% Output Arguments:
% =================
% R = V*substrate*(Ks+substrate)^(n-1) / ( L*(Ks*(1+inhibitor/Ki))^n + (Ks+substrate)^n )

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = V*substrate*(Ks+substrate)^(n-1) / ( L*(Ks*(1+inhibitor/Ki))^n + (Ks+substrate)^n );

