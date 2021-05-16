function [R] = kin_degradation(kdeg,substrate)
% kin_degradation: Linear degradation kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_degradation(kdeg,substrate)
%
% Output Arguments:
% =================
% R = kdeg*substrate

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

R = kdeg*substrate;

