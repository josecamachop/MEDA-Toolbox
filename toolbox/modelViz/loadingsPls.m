
function [P,W,Q] =loadingsPls(x,y,varargin)

% Compute and plot loadings in PLS. This routine is deprecated and superseded 
% by loadings.m (please, use the latter)
%
% P =loadingsPls(x,y) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
%
% Optional INPUTS (parameters):
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 1:rank(x)
%
% 'PreprocessingX': [1x1] preprocesing of the x-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)  
%
% 'PreprocessingY': [1x1] preprocesing of the y-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)   
%
% 'PlotType': str
%      'Scatter': scatterplot (by default)
%      'Bars': bar plot (by default)
%
% 'LoadingsType': str
%      'Weights': PLS weights
%      'Loadings': PLS loadings (by default)
%
% 'VarsLabel': [Mx1] name of the variables (numbers are used by default)
%
% 'VarsClass': [Mx1] groups for different visualization (a single group 
%   by default)
%
% 'BlurIndex': [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default).
%
%
% OUTPUTS:
%
% P: [MxA] X-block loadings
%
% W: [MxA] X-block weights
%
% Q: [OxA] Y-block loadings
%
%
% EXAMPLE OF USE: Random loadings: bar and scatter plot of loadings
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% A = cell(1, 10);
% 
% for i = 1:10
%     A{i} = ['A_{', num2str(i), '}'];
% end
% 
% loadingsPls(X,Y,'LVs',1);
% [P,W,Q] =loadingsPls(X,Y,'LVs',1:3,'VarsLabel',A);
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
O = size(y, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
LVS = 1:rank(x);
addParameter(p,'LVs',LVS);  
addParameter(p,'PreprocessingX',2);
addParameter(p,'PreprocessingY',2);
addParameter(p,'PlotType','Scatter');
addParameter(p,'LoadingsType','Loadings'); 
addParameter(p,'VarsLabel',1:M);
addParameter(p,'VarsClass',ones(M,1));   
addParameter(p,'BlurIndex',1);     
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
loadingstype = p.Results.LoadingsType;
plottype = p.Results.PlotType;
label = p.Results.VarsLabel;
classes = p.Results.VarsClass;
blur = p.Results.BlurIndex;

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''LVs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: parameter ''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: parameter ''VarsClass'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,'Preprocessing',prepx);
ycs = preprocess2D(y,'Preprocessing',prepy);

model = simpls(xcs,ycs,'LVs',lvs); 
W = model.weights;
P = model.loads;
Q = model.yloads;

%% Show results


if strcmp(loadingstype,'Weights')
    Pt = W;
    text = 'Weights';
elseif strcmp(loadingstype,'Loadings')
    Pt = P;
    text = 'X-block loadings';
end

if length(lvs) == 1 || strcmp(plottype,'Bars')
    for i=1:length(lvs)
        plotVec(Pt(:,i),  'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('%s LV %d',text,lvs(i))});
    end
elseif strcmp(plottype,'Scatter')
    for i=1:length(lvs)-1
        for j=i+1:length(lvs)
            plotScatter([Pt(:,i),Pt(:,j)], 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{sprintf('%s LV %d',text,lvs(i)),sprintf('%s LV %d',text,lvs(j))}', 'BlurIndex',blur);
        end
    end
end
        