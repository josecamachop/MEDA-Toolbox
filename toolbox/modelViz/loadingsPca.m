
function P = loadingsPca(x,varargin)


% Compute and plot loadings in PCA. This routine is deprecated and superseded 
% by loadings.m (please, use the latter)
%
% P = loadingsPca(x) % minimum call
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
% 'Preprocessing': [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
%
% 'PlotType': str
%      'Scatter': scatterplot (by default)
%      'Bars': bar plot (by default)
%
% 'VarsLabel': [Mx1] name of the variables (numbers are used by default)
%
% 'VarsClass': [Mx1] groups for different visualization (a single group 
%   by default)
%
% 'BlurIndex': [1x1] to avoid blur when adding labels. It reflects the
%   minimum distance with other points where a label is allowed to be 
%   visualized. For a value of 0, all labels are printed, while for a 
%   large value only uncluttered labels are printed. When Inf is chosen, 
%   only indices as visualized (by default 1).
%
%
% OUTPUTS:
%
% P: [MxA] loadings
%
%
% EXAMPLE OF USE: Scatter plot of random scores
%
% A = cell(1, 10);
% 
% for i = 1:10
%     A{i} = ['A_{', num2str(i), '}'];
% end
% 
% X = simuleMV(20,10,'LevelCorr',8);
% P = loadingsPca(X,'PCs',1:3,'VarsLabel',A);
%
%
% EXAMPLE OF USE: Line plot of random scores
%
% X = real(adicov(randn(10,10).^19,randn(100,10),10));
% P = loadingsPca(X,'PCs',1:3,'PlotType','Bars');
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 14/Aug/2025
% Dependencies: Matlab R2017b, MEDA v1.9
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
PCS = 1:rank(x);
addParameter(p,'PCs',PCS);  
addParameter(p,'Preprocessing',2);
addParameter(p,'PlotType','Scatter');
addParameter(p,'VarsLabel',1:M);
addParameter(p,'VarsClass',ones(M,1));   
addParameter(p,'BlurIndex',1);     
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
pcs = p.Results.PCs;
prep = p.Results.Preprocessing;
plottype = p.Results.PlotType;
label = p.Results.VarsLabel;
classes = p.Results.VarsClass;
blur = p.Results.BlurIndex;

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
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: parameter ''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: parameter ''VarsClass'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''PCs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,'Preprocessing',prep);
model = pcaEig(xcs,'PCs',pcs);
P = model.loads;
T = model.scores;

%% Show results

if length(pcs) == 1 || strcmp(plottype,'Bars')
    for i=1:length(pcs)
        plotVec(P(:,i), 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('Loadings PC %d (%.0f%%)',pcs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2)))});
    end
elseif strcmp(plottype,'Scatter')
    for i=1:length(pcs)-1
        for j=i+1:length(pcs)
            plotScatter([P(:,i),P(:,j)], 'EleLabel',label,'ObsClass' ,classes, 'XYLabel',{sprintf('Loadings PC %d (%.0f%%)',pcs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2))),sprintf('Loadings PC %d (%.0f%%)',pcs(j),100*sum(T(:,j).^2)/sum(sum(xcs.^2)))}', 'BlurIndex',blur);
        end
    end
end
        