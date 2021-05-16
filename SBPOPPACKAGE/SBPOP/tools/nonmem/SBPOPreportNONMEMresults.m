function SBPOPreportNONMEMresults(projectPath)
% SBPOPreportNONMEMresults: parses the results of a NONMEM run and reports
% them in a similar manner as in the MONOLIX pop_parameters.txt file.
%
% The function saves a text file version in the projectPath/RESULTS
% folder. Additionally, the result is shown in the command window.
%
% USAGE:
% ======
% SBPOPreportNONMEMresults(projectPath)
%
% projectPath:   path to the .nmctl NONMEM project file

% Information:
% ============
% Copyright (C) 2012 Novartis Pharma AG
%
% This program is Free Open Source Software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = parseNONMEMresultsSBPOP(projectPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start output text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUTPUT = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUTPUT  = sprintf('%s===========================================================================\n',OUTPUT);
OUTPUT  = sprintf('%s    Summary results\n',OUTPUT);
[~,project] = fileparts(x.path);
OUTPUT  = sprintf('%s    Project: %s\n',OUTPUT,project);
OUTPUT  = sprintf('%s===========================================================================\n',OUTPUT);
OUTPUT  = sprintf('%s\n',OUTPUT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print termination information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(x.termination_info),
    OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
    method = sprintf('%s,',x.PROJECTINFO.METHOD{:});
    OUTPUT = sprintf('%sTermination information (Method(s): %s)\n',OUTPUT,method(1:end-1));
    OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
    for k=1:length(x.termination_info),
        OUTPUT = sprintf('%s%s\n',OUTPUT,x.termination_info{k});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if major problems with the fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isnan(x.objectivefunction.OBJ),
    OUTPUT = sprintf('%sMajor problems with the project results. Please check!\n',OUTPUT);
    % Save the text
    filename = sprintf('%s/RESULTS/project_results.txt',projectPath);
    fid = fopen(filename,'w');
    fprintf(fid,'%s',OUTPUT);   
    fclose(fid);
    % Print out in command window
    disp(OUTPUT)    
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the fixed effects (still transformed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
OUTPUT  = sprintf('%sName                     Value      stderr     RSE (%%)          95%%CI\n',OUTPUT);
OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
fe = x.rawParameterInfo.fixedEffects;
names = {};
values = {};
stderrs = {};
RSEs = {};
for k=1:length(fe.names)
    % Names
    if strcmp(fe.trans{k},'(phi)'),
        names{k} = fe.names{k};
    elseif strcmp(fe.trans{k},'exp(phi)'),
        names{k} = ['log(' fe.names{k} ')'];
    elseif strcmp(fe.trans{k},'exp(phi)./(1+exp(phi))'),
        names{k} = ['logit(' fe.names{k} ')'];
    else
        error('Unknown transformation');
    end
    % Values
    values{k} = sprintf('%1.4g',fe.values(k));
    if fe.estimated(k) && ~isnan(fe.stderr(k)),
        % Stderrs
        stderrs{k} = sprintf('%1.4g',fe.stderr(k));
        % RSEs
        RSEs{k} = sprintf('%1.4g',fe.rse(k));
    else
        stderrs{k} = '-';
        RSEs{k} = '-';
    end
end
for k=1:length(names),
    text = sprintf('%s%s%s%s\n',postFillCharSB(names{k},18,' '),...
        preFillCharSB(values{k},12,' '), ...
        preFillCharSB(stderrs{k},12,' '), ...
        preFillCharSB(RSEs{k},12,' '));
    OUTPUT = sprintf('%s%s',OUTPUT,text);
end
OUTPUT  = sprintf('%s\n',OUTPUT);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the fixed effects (back transformed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe = x.rawParameterInfo.fixedEffects;
names = {};
values = {};
stderrs = {};
RSEs = {};
for k=1:length(fe.names)
    % Names
    names{k} = fe.names{k};
    % Value
    phi = fe.values(k);
    value = eval(fe.trans{k});
    % sample standard error
    if fe.stderr(k)==0,
        stderr = 0;
    else
        phi = fe.values(k)+fe.stderr(k)*randn(1,100000);
        stderr = std(eval(fe.trans{k}));
    end
    
    if fe.estimated(k) && ~isnan(stderr),
        % RSEs
        rse = sprintf('%1.4g*',100*stderr/value);
        % Stderrs
        stderr = sprintf('%1.4g*',stderr);
    else
        stderr = '-';
        rse = '-';
    end
    
    % Values
    values{k} = sprintf('%1.4g',value);
    % Stderrs
    stderrs{k} = stderr;
    % RSEs
    RSEs{k} = rse;
end
for k=1:length(names),
    text = sprintf('%s%s%s%s\n',postFillCharSB(names{k},18,' '),...
        preFillCharSB(values{k},12,' '), ...
        preFillCharSB(stderrs{k},12,' '), ...
        preFillCharSB(RSEs{k},12,' '));
    OUTPUT = sprintf('%s%s',OUTPUT,text);
end
OUTPUT  = sprintf('%s                                    (*approximation by sampling)\n',OUTPUT);
OUTPUT  = sprintf('%s\n',OUTPUT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe = x.rawParameterInfo.covariate;
names = {};
values = {};
stderrs = {};
RSEs = {};
CI = {};
for k=1:length(fe.names)
    names{k} = fe.names{k};
    values{k} = sprintf('%1.4g',fe.values(k));
    if fe.estimated(k) && ~isnan(fe.stderr(k)),
        stderrs{k} = sprintf('%1.4g',fe.stderr(k));
        RSEs{k} = sprintf('%1.4g',fe.rse(k));
        CI{k} = sprintf('[%1.2f %1.2f]',fe.values(k)-1.96*fe.stderr(k),fe.values(k)+1.96*fe.stderr(k));
    else
        stderrs{k} = '-';
        RSEs{k} = '-';
        CI{k} = '-';
    end
end
for k=1:length(names),
    text = sprintf('%s%s%s%s%s\n',postFillCharSB(names{k},18,' '),...
        preFillCharSB(values{k},12,' '), ...
        preFillCharSB(stderrs{k},12,' '), ...
        preFillCharSB(RSEs{k},12,' '), ...
        preFillCharSB(CI{k},20,' '));
    OUTPUT = sprintf('%s%s',OUTPUT,text);
end
if length(fe.names) > 0,
    OUTPUT  = sprintf('%s\n',OUTPUT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the random effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe = x.rawParameterInfo.randomEffects;
names = {};
values = {};
stderrs = {};
RSEs = {};
for k=1:length(fe.names)
    names{k} = fe.names{k};
    values{k} = sprintf('%1.4g',fe.values(k));
    if fe.estimated(k) && ~isnan(fe.stderr(k)),
        stderrs{k} = sprintf('%1.4g',fe.stderr(k));
        RSEs{k} = sprintf('%1.4g',fe.rse(k));
    else
        stderrs{k} = '-';
        RSEs{k} = '-';
    end
end
for k=1:length(names),
    text = sprintf('%s%s%s%s\n',postFillCharSB(names{k},18,' '),...
        preFillCharSB(values{k},12,' '), ...
        preFillCharSB(stderrs{k},12,' '), ...
        preFillCharSB(RSEs{k},12,' '));
    OUTPUT = sprintf('%s%s',OUTPUT,text);
end
if length(fe.names) > 0,
    OUTPUT  = sprintf('%s\n',OUTPUT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(x.rawParameterInfo.correlation.names),
    fe = x.rawParameterInfo.correlation;
    names = {};
    values = {};
    stderrs = {};
    RSEs = {};
    CI = {};
    for k=1:length(fe.names)
        names{k} = fe.names{k};
        values{k} = sprintf('%1.4g',fe.values(k));
        if ~isnan(fe.stderr(k)),
            stderrs{k} = sprintf('%1.4g',fe.stderr(k));
            RSEs{k} = sprintf('%1.4g',fe.rse(k));
            CI{k} = sprintf('[%1.2f %1.2f]',fe.values(k)-1.96*fe.stderr(k),fe.values(k)+1.96*fe.stderr(k));
        else
            stderrs{k} = '-';
            RSEs{k} = '-';
            CI{k} = '-';
        end
    end
    for k=1:length(names),
        text = sprintf('%s%s%s%s%s\n',postFillCharSB(names{k},18,' '),...
            preFillCharSB(values{k},12,' '), ...
            preFillCharSB(stderrs{k},12,' '), ...
            preFillCharSB(RSEs{k},12,' '), ...
            preFillCharSB(CI{k},20,' '));
        OUTPUT = sprintf('%s%s',OUTPUT,text);
    end
    OUTPUT  = sprintf('%s\n',OUTPUT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe = x.rawParameterInfo.errorParameter;
names = {};
values = {};
stderrs = {};
RSEs = {};
for k=1:length(fe.names)
    names{k} = fe.names{k};
    values{k} = sprintf('%1.4g',fe.values(k));
    if ~isnan(fe.stderr(k)),
        stderrs{k} = sprintf('%1.4g',fe.stderr(k));
        RSEs{k} = sprintf('%1.4g',fe.rse(k));
    else
        stderrs{k} = '-';
        RSEs{k} = '-';
    end
end
for k=1:length(names),
    text = sprintf('%s%s%s%s\n',postFillCharSB(names{k},18,' '),...
        preFillCharSB(values{k},12,' '), ...
        preFillCharSB(stderrs{k},12,' '), ...
        preFillCharSB(RSEs{k},12,' '));
    OUTPUT = sprintf('%s%s',OUTPUT,text);
end
OUTPUT  = sprintf('%s\n',OUTPUT);


if ~isempty(x.parameters.correlationmatrix),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correlation fE and beta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
    OUTPUT = sprintf('%sCorrelation of fixed effects and covariate coefficients\n',OUTPUT);
    OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
    names = [x.rawParameterInfo.fixedEffects.names x.rawParameterInfo.covariate.names];
    ix_all = [];
    for k=1:length(names),
        ix_all(k) = strmatchSB(names{k},x.parameters.names);
    end
    cor = x.parameters.correlationmatrix(ix_all,ix_all);
    for row=1:length(cor),
        rowtext = postFillCharSB(names{row},18,' ');
        for col=1:row,
            rowtext = sprintf('%s%s',rowtext,preFillCharSB(sprintf('%1.2g',0.01*round(100*cor(row,col))),7,' '));
        end
        OUTPUT = sprintf('%s%s\n',OUTPUT,rowtext);
    end
    OUTPUT  = sprintf('%s\n',OUTPUT);
    eigM    = eig(cor);
    eigMmin = min(eigM);
    eigMmax = max(eigM);
    OUTPUT  = sprintf('%sEigenvalues (min, max, max/min): %1.2f  %1.2f  %1.2f\n',OUTPUT,eigMmin,eigMmax,eigMmax/eigMmin);
    OUTPUT  = sprintf('%s\n',OUTPUT);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correlation of omegas and error parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
    OUTPUT = sprintf('%sCorrelation of random effects (variances) and error parameters\n',OUTPUT);
    OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
    names = [strrep(x.rawParameterInfo.randomEffects.names,'omega','omega2') x.rawParameterInfo.errorParameter.names];
    ix_all = [];
    for k=1:length(names),
        ix_all(k) = strmatchSB(names{k},x.parameters.names,'exact');
    end
    cor = x.parameters.correlationmatrix(ix_all,ix_all);
    for row=1:length(cor),
        rowtext = postFillCharSB(names{row},18,' ');
        for col=1:row,
            rowtext = sprintf('%s%s',rowtext,preFillCharSB(sprintf('%1.2g',0.01*round(100*cor(row,col))),7,' '));
        end
        OUTPUT = sprintf('%s%s\n',OUTPUT,rowtext);
    end
    OUTPUT  = sprintf('%s\n',OUTPUT);
    eigM    = eig(cor);
    eigMmin = min(eigM);
    eigMmax = max(eigM);
    OUTPUT  = sprintf('%sEigenvalues (min, max, max/min): %1.2f  %1.2f  %1.2f\n',OUTPUT,eigMmin,eigMmax,eigMmax/eigMmin);
    OUTPUT  = sprintf('%s\n',OUTPUT);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correlation of correlations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    names = [strrep(x.rawParameterInfo.correlation.names,'corr(','omega2(')];
    if length(names) > 1,
        OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
        OUTPUT = sprintf('%sCorrelation of random effect covariances\n',OUTPUT);
        OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
        ix_all = [];
        for k=1:length(names),
            ix_all(k) = strmatchSB(names{k},x.parameters.names);
        end
        cor = x.parameters.correlationmatrix(ix_all,ix_all);
        for row=1:length(cor),
            rowtext = postFillCharSB(names{row},18,' ');
            for col=1:row,
                rowtext = sprintf('%s%s',rowtext,preFillCharSB(sprintf('%1.2g',0.01*round(100*cor(row,col))),7,' '));
            end
            OUTPUT = sprintf('%s%s\n',OUTPUT,rowtext);
        end
        OUTPUT  = sprintf('%s\n',OUTPUT);
        eigM    = eig(cor);
        eigMmin = min(eigM);
        eigMmax = max(eigM);
        OUTPUT  = sprintf('%sEigenvalues (min, max, max/min): %1.2f  %1.2f  %1.2f\n',OUTPUT,eigMmin,eigMmax,eigMmax/eigMmin);
        OUTPUT  = sprintf('%s\n',OUTPUT);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the OFV / AIC / BIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
METHOD = x.PROJECTINFO.METHOD{end};
OUTPUT = sprintf('%sObjective function (%s)\n',OUTPUT,METHOD);
if strcmp(METHOD,'SAEM'),
    OUTPUT = sprintf('%sThe SAEM objective function should not be used for statistical testing.\n',OUTPUT);
    OUTPUT = sprintf('%sPlease consider the use of the IMPORTANCESAMPLING option!\n',OUTPUT);
end
OUTPUT = sprintf('%s---------------------------------------------------------------------------\n',OUTPUT);
OUTPUT = sprintf('%sOFV:    %g\n',OUTPUT,x.objectivefunction.OBJ);
OUTPUT = sprintf('%sAIC:    %g\n',OUTPUT,x.objectivefunction.AIC);
OUTPUT = sprintf('%sBIC:    %g\n',OUTPUT,x.objectivefunction.BIC);
OUTPUT  = sprintf('%s\n',OUTPUT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = sprintf('%s/RESULTS/project_results.txt',projectPath);
fid = fopen(filename,'w');
fprintf(fid,'%s',OUTPUT);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print out in command window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(OUTPUT)
