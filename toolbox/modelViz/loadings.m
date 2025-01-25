
function figH = loadings(model,varargin)

% Compute and plot loadings.
%
% figH =loadings(model) % minimum call
%
%
% INPUTS:
%
% model (structure): structure with model parameters. 
%   lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%       first two LVs).
%   loads: [MxA] model parameters.
%   scores: [NxA] data scores. 
%
%
% Optional Inputs (parameter):
%
% 'PlotType': str
%      'Scatter': scatterplot (by default)
%      'Bars': bar plot (by default)
%
% 'Title': (str) title for the plots. Empty by default;
%
% 'VarsLabel': [Mx1] name of the variables (numbers are used by default)
%
% 'VarsClass': [Mx1] groups for different visualization (a single group 
%   by default)
%
% 'BlurIndex': [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default)
%
% 'Color': Choose a color for your data.  
%   - 'hsv' for hsv palette 
%   - 'parula' for parula palette
%   - 'okabeIto' for color blindness (by default for multiple classes) 
%
%
% OUTPUTS:
%
% figH: set of figure handles
%
%
% EXAMPLE OF USE: Random data
% 
% A = cell(1, 10);
% 
% for i = 1:10
%     A{i} = ['A_{', num2str(i), '}'];
% end
% 
% X = simuleMV(20,10,'LevelCorr',8);
% model = pcaEig(X,'PCs',1:2);
% 
% P = loadings(model,'VarsLabel',A);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 15/Jan/2025
%
% Copyright (C) 2025  University of Granada, Granada
%  
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(model.scores, 1);
M = size(model.loads, 1);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'PlotType','Scatter');
addParameter(p,'Title',' ');
addParameter(p,'VarsLabel',1:M);
addParameter(p,'VarsClass',ones(M,1));   
addParameter(p,'BlurIndex',1);   
addParameter(p,'Color',[]);  
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
plottype = p.Results.PlotType;
tit = p.Results.Title;
label = p.Results.VarsLabel;
classes = p.Results.VarsClass;
blur = p.Results.BlurIndex;
color = p.Results.Color;

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
assert (isequal(size(label), [M 1]), 'Dimension Error: parameter''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: parameter ''VarsClass'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  

%% Main code

P = model.loads;


%% Show results

if ~isfield(model,'type') || strcmp(model.type,'PCA')
    dim = 'PC';
elseif strcmp(model.type,'PLS')
    dim = 'LV';
else
    dim = 'PC';
end

figH = [];
if length(model.lvs) == 1 || strcmp(plottype,'Bars')
    for i=1:length(model.lvs)
        figH = [figH plotVec(P(:,i), 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('Loadings %s %d',dim,model.lvs(i))},'Color',color)];
        title(tit);
    end
elseif strcmp(plottype,'Scatter')
    for i=1:length(model.lvs)-1
        for j=i+1:length(model.lvs)
            figH = [figH plotScatter([P(:,i),P(:,j)], 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{sprintf('Loadings %s %d',dim,model.lvs(i)),sprintf('Loadings %s %d',dim,model.lvs(j))}','BlurIndex',blur,'Color',color)];
            title(tit);
        end
    end
end
        