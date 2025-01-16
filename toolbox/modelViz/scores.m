
function figH = scores(model,varargin)

% Compute and plot scores.
%
% figH = scores(model) % minimum call
%
%
% INPUTS:
%
% model (structure): structure with model parameters. 
%   var: [1x1] Total variance
%   lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%       first two LVs).
%   av: [1xM] centering parameters. 
%   sc: [1xM] scaling parameters. 
%   loads: [MxA] model parameters.
%   scores: [NxA] data scores.
%
%
% Optional INPUTS(parameters):
%
% 'ObsTest': [LxM] data set with the observations to be visualized in the model
%   space. By default, model.scores are plotted.
%
% 'PlotType': str
%      'Scatter': scatterplot (by default)
%      'Bars': bar plot (by default)
%
% 'PlotCal': bool
%      false: plot only test data
%      true: plot both calibration and test (by default)
%
% 'Title': (str) title for the plots. Empty by default;
%
% 'ObsLabel': [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% 'ObsClass': [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
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
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(20,10,'LevelCorr',8);
% 
% model.lvs = 1:3;
% [Xcs,model.av,model.sc] = preprocess2D(X);
% model.var = trace(Xcs'*Xcs);
% model = pcaEig(Xcs,'PCs',model.lvs);
% 
% scores(model);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% nObs = 100;
% nVars = 10;
% X = simuleMV(nObs,nVars,'LevelCorr',8);
% 
% model.lvs = 1:2;
% [Xcs,av,sc] = preprocess2D(X);
% model = pcaEig(Xcs,'PCs',model.lvs);
% model.var = trace(Xcs'*Xcs);
% model.av = av;
% model.sc = sc;
% 
% nObst = 10;
% test = simuleMV(nObst,nVars,'LevelCorr',8,'Covar',corr(X)*(nObst-1)/(nObs-1));
% 
% scores(model,'ObsTest',test);
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 12/Jan/2025
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
M = size(model.loads,1);
A = length(model.lvs);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'ObsTest',[]);
addParameter(p,'PlotType','Scatter');
addParameter(p,'PlotCal',true);
addParameter(p,'Title',' ');  
addParameter(p,'BlurIndex',1);
addParameter(p,'ObsLabel',[]);
addParameter(p,'ObsClass',[]);
addParameter(p,'Color',[]);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility

test = p.Results.ObsTest;
plottype = p.Results.PlotType;
plotcal = p.Results.PlotCal;
tit = p.Results.Title;
blur = p.Results.BlurIndex;
label = p.Results.ObsLabel;
classes = p.Results.ObsClass;
color = p.Results.Color;

L = size(test, 1);

if plotcal
    K = N+L;
else
    K = L;
end

if isempty(label) 
    if plotcal
        label = [1:N 1:L]; 
    else
        label = 1:L;
    end
end
if isempty(classes)
    if plotcal
        classes = [ones(N,1);2*ones(L,1)];  
    else 
        classes = ones(L,1); 
    end
end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: parameter ''ObsTest'' must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (isequal(size(label), [K 1]), 'Dimension Error: parameter ''ObsLabel'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: parameter ''ObsClass'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  

%% Main code

T = model.scores;
d = diag(T'*T);

if isfield(model,'scoresV')
    T = model.scoresV;
end

if ~isempty(test)
    if isfield(model,'av')
        testcs = preprocess2Dapp(test,model.av,'Scale',model.sc);
    else
        testcs = test;
    end
    TT = testcs*model.loads;
else
    TT = [];
end


%% Show results

if plotcal
    Tt = [T;TT];
else
    Tt = TT;
end

if ~isfield(model,'type') || strcmp(model.type,'PCA')
    dim = 'PC';
elseif strcmp(model.type,'PLS')
    dim = 'LV';
else
    dim = 'PC';
end

figH = [];
if strcmp(plottype,'Bars') || A == 1
    for i=1:length(model.lvs)
        figH = [figH plotVec(Tt(:,i), 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('Scores %s %d (%.0f%%)',dim,model.lvs(i),100*d(i)/model.var)},'Color',color);];
        title(tit);
    end
elseif strcmp(plottype,'Scatter')
    for i=1:length(model.lvs)-1
        for j=i+1:length(model.lvs)
            figH = [figH plotScatter([Tt(:,i),Tt(:,j)],'EleLabel',label,'ObsClass',classes,'XYLabel',{sprintf('Scores %s %d (%.0f%%)',dim,model.lvs(i),100*d(i)/model.var),sprintf('Scores %s %d (%.0f%%)',dim,model.lvs(j),100*d(j)/model.var)}','BlurIndex',blur,'Color',color)];
            title(tit);
        end
    end
end
        