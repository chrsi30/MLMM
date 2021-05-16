function [SBstructure,errorMsg] = SBimportSBML(SBMLmodelFile)
% SBimportSBML 
% imports a SBML model using the TranslateSBML function from libSBML
% Supported SBML levels: 1 and 2
%
% SBMLmodelFile: SBMLmodelfilename.xml (string)
%
% SBstructure: empty if error occurred

% Information:
% ============
% Copyright (C) 2005-2007 Fraunhofer Chalmers Centre, Gothenburg, Sweden
% Main author: Henning Schmidt
% 
% Changes for the SBTOOLBOX2:
% 1/1/2008  Henning Schmidt, henning@sbtoolbox2.org
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
% USA.

academicWarningSB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorMsg = '';
SBstructure = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBML -> MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the TranslateSBML function from libSBML
try
    SBMLmodel = TranslateSBML(SBMLmodelFile);
catch
    errTxt = sprintf('An error during the SBML import using the SBML Toolbox occurred.\nThis might be due to the fact that you have an outdated installation\nof the SBML Toolbox. Please have a look at www.sbtoolbox2.org for which\nversion of the SBML Toolbox is required.\n\nLast error message:\n\n%s',lasterr);
    error(errTxt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERT INTO OWN SB STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(SBMLmodel.typecode,'SBML_MODEL'),
    errorMsg = 'Model is not an SBML model';
    return;
end
if SBMLmodel.SBML_level == 1,
    % Convert level 1
    [SBstructure,errorMsg] = SBconvertSBML1(SBMLmodel);
elseif SBMLmodel.SBML_level == 2 && SBMLmodel.SBML_version == 1,
    % Convert level 2
    [SBstructure,errorMsg] = SBconvertSBML2(SBMLmodel);
else
    % Not supported SBML level
    errorMsg = sprintf('Level %d Version %d of SBML not yet implemented',SBMLmodel.SBML_level,SBMLmodel.SBML_version);
end
