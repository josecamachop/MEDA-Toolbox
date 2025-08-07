
function [yvar,tvar] = varPls(x,y,varargin)

% Variability captured in terms of the number of LVs.
%
% varPls(x,y) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
%
% Optional INPUTS (parameter):
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x). The value for 0 PCs is
%   added at the begining if not specified.
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
% 'PlotScores': bool
%       false: Residual Variance in Y (by default)
%       true: Residual Variance in Y and Scores  
%
%
% OUTPUTS:
%
% yvar: [Ax1] Percentage of captured variance of Y.
%
% tvar: [Ax1] Percentage of captured variance of the scores.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% lvs = 0:10;
% xvar = varPls(X,Y,'LVs',lvs);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 25/Jan/2025
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
O = size(y, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'LVs',0:rank(x));   
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);   
addParameter(p,'PlotScores',false);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;
scoresplot = p.Results.PlotScores;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
A = length(lvs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''LVs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique([0 lvs]);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LVs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,'Preprocessing',prepx); 
ycs = preprocess2D(y,'Preprocessing',prepy); 

%model = kernelpls(xcs'*xcs,xcs'*ycs,'LVs',1:max(lvs));
model = simpls(xcs,ycs,'LVs',1:max(lvs));
Q = model.yloads;
R = model.altweights;

lvs(find(lvs>size(R,2))) = [];

totalVt = sum(sum(xcs.^2));
tvar = ones(length(lvs),1);
totalVy = sum(sum(ycs.^2));
yvar = ones(length(lvs),1);
for i=1:length(lvs)
    tvar(i) = tvar(i) - sum(eig(R(:,1:lvs(i))'*xcs'*xcs*R(:,1:lvs(i))))/totalVt;
    yvar(i) = yvar(i) - sum(eig(Q(:,1:lvs(i))*R(:,1:lvs(i))'*xcs'*xcs*R(:,1:lvs(i))*Q(:,1:lvs(i))'))/totalVy;
end
    
%% Show results

if scoresplot
    plotVec([yvar tvar],'EleLabel',lvs,'XYLabel',{'#LVs','% Residual Variance'},'PlotType','Lines','VecLabel',{'Y','Scores'});
else
    plotVec(yvar,'EleLabel',lvs,'XYLabel',{'#LVs','% Residual Variance in Y'},'PlotType','Lines');
end

