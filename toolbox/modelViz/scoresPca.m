
function [T,TT] = scoresPca(x,varargin)


% Compute and plot scores in PCA. This routine is deprecated and superseded 
% by scores.m (please, use the latter)
%
% T = scoresPca(x) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
%
% Optional INPUTS (parameters):
%
% 'PCs': [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 1:rank(xcs)
%
% 'ObsTest': [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% 'Preprocessing': [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
%
% 'PlotType': str
%      'Scatter': scatterplot (by default)
%      'Bars': bar plot (by default)
%
% 'PlotCal': bool
%      false: plot only test data
%      true: plot both calibration and test (by default)
%
% 'ObsLabel': [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% 'ObsClass': [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
%
% 'BlurIndex': [1x1] to avoid blur when adding labels. It reflects the
%   minimum distance (normalized to [0,1]) where a cluttered label is 
%   allowed to be visualized. For a value of 0, no cluttered labels are 
%   printed, while for a value of 1 all labels are printed, and thus the 
%   highest blur. By default 0.5 is chosen.
%
%
% OUTPUTS:
%
% T: [NxA] calibration scores
%
% TT: [NxA] test scores
%
%
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(20,10,'LevelCorr',8);
% T = scoresPca(X,'PCs',1:3);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% nObs = 100;
% nVars = 10;
% X = simuleMV(nObs,nVars,'LevelCorr',8);
% 
% nObst = 10;
% test = simuleMV(nObst,nVars,'LevelCorr',8,'Covar',corr(X)*(nObst-1)/(nObs-1));
% 
% scoresPca(X,'PCs',1,'ObsTest',test);
% scoresPca(X,'PCs',1:2,'ObsTest',test);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 13/Jan/2025
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
N = size(x, 1);
M = size(x, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'PCs',1:rank(x));
addParameter(p,'PlotType','Scatter');
addParameter(p,'PlotCal',true);   
addParameter(p,'ObsTest',[]);  
addParameter(p,'Preprocessing',2);
addParameter(p,'ObsLabel',[]);
addParameter(p,'ObsClass',[]);
addParameter(p,'BlurIndex',0.5);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
test = p.Results.ObsTest;
plottype = p.Results.PlotType;
plotcal = p.Results.PlotCal;
prep = p.Results.Preprocessing;
pcs = p.Results.PCs;
label = p.Results.ObsLabel;
classes = p.Results.ObsClass;
blur = p.Results.BlurIndex;

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

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
A = length(pcs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''PCs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''PCs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: parameter ''ObsTest'' must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: parameter ''ObsLabel'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: parameter ''ObsClass'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''PCs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[xcs,m,sd] = preprocess2D(x,'Preprocessing',prep);
model = pcaEig(xcs,'PCs',pcs);
T = model.scores;

if ~isempty(test)
    testcs = preprocess2Dapp(test,m,'Scale',sd);
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

if strcmp(plottype,'Bars') || A == 1
    for i=1:length(pcs)
        plotVec(Tt(:,i), 'EleLabel', label, 'ObsClass', classes, 'XYLabel', {'',sprintf('Scores PC %d (%.0f%%)',pcs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2)))});
    end
elseif strcmp(plottype,'Scatter')
    for i=1:length(pcs)-1
        for j=i+1:length(pcs)
            plotScatter([Tt(:,i),Tt(:,j)], 'EleLabel', label, 'ObsClass', classes, 'XYLabel', {sprintf('Scores PC %d (%.0f%%)',pcs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2))), sprintf('Scores PC %d (%.0f%%)',pcs(j),100*sum(T(:,j).^2)/sum(sum(xcs.^2)))}', 'BlurIndex', blur);
        end
    end
end
        