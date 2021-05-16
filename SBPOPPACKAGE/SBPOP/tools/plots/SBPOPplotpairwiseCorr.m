function [] = SBPOPplotpairwiseCorr(data,OPTIONS)
% [DESCRIPTION]
% This function plots the pairwise correlation between variables passed in
% columns of a matrix or passed as a dataset.
%
% [SYNTAX]
% [] = SBPOPplotpairwiseCorr(data)
% [] = SBPOPplotpairwiseCorr(data,OPTIONS)
%
% [INPUT]
% data:         Matrix or dataset. Each column corresponds to a variable
%               and each row to a sample
% OPTIONS:      MATLAB structure with optional arguments
%
%                   options.names:     cell-array with variable names. In
%                       case of data as dataset names will be taken from the
%                       header but can be overwritten with this option. If
%                       variable values are provided in a matrix, it is better
%                       to provide the names using this option
%                   options.LogFlag:   =1 do log transform the variables,
%                                      =0 do not transform (default: 0)
%                   options.CorrThres:   Value between 0 and 1 indicating the
%                       threshold for the Pearson correlation coefficient from which on a
%                       different color as background should be used (default:
%                       0.3)
%                   options.AxisColor:  [r g b] numeric values to use as color 
%                       for background if Corr>CorrThres =1 do log transform
%                       the variables (default: [1 0.6 0.6])
%
% [OUTPUT]
% Pairwise correlation plots.
%
% [AUTHOR]
% Original author: Andy Stein
% Adaptation to SBPOP: Henning Schmidt
%
% [DATE]
% 8th February 2013
%
% [PLATFORM]
% Windows XP Engine, MATLAB R2009a, MATLAB
%
% [KEYWORDS]
% MATLAB, SBPOP
% 
% [TOOLBOXES USED]
% Statistics Toolbox

% Information:
% ============
% Copyright � 2012 Novartis Pharma AG
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

% Check if dataset or if matrix
% If dataset ten use columnnames as default values for names
if strcmp(class(data),'dataset'),
    names = get(data,'VarNames');
    % Make a matrix out of it
    data = double(data);
elseif strcmp(class(data),'double'),
    names = {};
    for k=1:size(data,2),
        names{k} = sprintf('Data#%d',k);
    end
end

% Default settings
% Flag to take log of values
LogFlag = 0;
% Coloring
CorrThres = 0.3; % Pearson correlation coefficient threshold for different background color
AxisColor = [1 .6 .6]; %color of axis when 

% Get optional values if defined
try names = OPTIONS.names; catch, end;
try LogFlag = OPTIONS.LogFlag; catch, end;
try CorrThres = OPTIONS.CorrThres; catch, end;
try AxisColor = OPTIONS.AxisColor; catch, end;

% Subindex properties
Spacing = 0;
Padding = 0;
Margin  = .1;
    
clf;
n = length(names);
for ir=1:n
    ystr = names{ir};
    y = data(:,ir);
    if LogFlag==1
        y = log(y);
        ystr = {'log',ystr}; %#ok<*AGROW>
    end
    for ic=1:ir
        xstr = names{ic};
        x = data(:,ic);
        if LogFlag==1
            x = log(x);
            xstr = {'log',xstr};
        end
        ip = (ir-1)*n+ic;
        subaxis(n,n,ip,'Spacing',Spacing,'Padding',Padding,'Margin',Margin);
        if ir==ic %plot histogram 
            % Only do this if not only NaN values are present
            if ~isempty(find(isnan(x)==0)),
                [b xbin] = hist(x,20);
                h = bar(xbin,b,1);
                set(gca,'XLim',[min(x) max(x)]);
            end
            set(h,'FaceColor',.6*[1 1 1],'LineStyle','none')
            title(xstr,'Interpreter','none')
        else %plot correlation
            % Only do this is non NaN pairs do exist, which might not
            % always be the case and then will lead to an error
            if ~isempty(find(double(isnan(x))+double(isnan(y))==0)),
                % Determine pearsons coefficient of correlation
                [rho,pval] = corr(x,y,'type','Pearson','rows','pairwise','tail','both');

                % Plot
                if abs(rho)>CorrThres
                    optcorr.Color     = AxisColor;
                else
                    optcorr.Color     = 0.6*[1 1 1];
                end
                optcorr.LineColor = [0 0 0];
                optcorr.TitleType = 'none';
                optcorr.LineStyle = '-';
                optcorr.LineWidth = 2;
                plotcorrSBPOP(x,y,optcorr);
                
                xt = (min(x)+max(x))/2;
                yt = min(y)+.6*(max(y)-min(y));
                if abs(rho)>=0.01
                    if pval>=0.01,
                        str = sprintf('corr=%1.2f\np=%1.2f',rho,pval);
                    else
                        str = sprintf('corr=%1.2f\np<0.01',rho);
                    end
                else
                    str = '|corr|<0.01';
                end
                text(xt,yt,str,'Color',[0 0 0],'Hor','Center','Ver','Middle','FontWeight','Bold','Interpreter','none');
                
                set(gca,'XLim',[min(x) max(x)]);
                set(gca,'YLim',[min(y) max(y)]);
            end
        end
        
        if ic==1
            ylabel(ystr,'Interpreter','none');
        end
        set(gca,'YTick',[]);
        
        if ir==n
            xlabel(xstr,'Interpreter','none');
        end
        set(gca,'XTick',[]);
        
    end
end