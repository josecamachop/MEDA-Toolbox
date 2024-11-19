function [ypred,testypred] = predPls(x,y,varargin)

% Compute and plot prediction in PLS.
%
% ypred = predPls(x,y) % minimum call
%
% See also: simpls, kernelpls
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
% 'ObsTest': [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
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
%
% OUTPUTS:
%
% ypred: [NxO] calibration data prediction of y-block
%
% testypred: [LxO] test data prediction of y-block
%
%
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% ypred = predPls(X,Y,'LVs',1:3);
%
%
% EXAMPLE OF USE: Calibration and Test for the prediction of two variables
%
% nObs = 100;
% nVars = 10;
% X = simuleMV(nObs,nVars,'LevelCorr',6);
% Y = 0.1*randn(nObs,2) + X(:,1:2);
% 
% nObst = 10;
% test = simuleMV(nObst,nVars,'LevelCorr',6,'Covar',corr(X)*(nObst-1)/(nObs-1));
% Ytest = 0.1*randn(nObst,2) + test(:,1:2);
% 
% [ypred,testypred] = predPls(X,Y,'LVs',1:2,'ObsTest',test);
% plotScatter([Ytest(:,1),testypred(:,1)],'XYLabel',{'Real var1' 'Predicted var1'});
% plotScatter([Ytest(:,2),testypred(:,2)],'XYLabel',{'Real var2' 'Predicted var2'});
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 18/Nov/2024
%
% Copyright (C) 2024  University of Granada, Granada
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
addParameter(p,'LVs',1:rank(x)); 
addParameter(p,'ObsTest',[]);
L = size('ObsTest', 1);
addParameter(p,'PreprocessingX',2);  
addParameter(p,'PreprocessingY',2); 
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;
test = p.Results.ObsTest;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
L = size(test, 1);
K = N+L;


% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''LVs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: parameter ''ObsTest'' must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
  
% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[xcs,m,sd] = preprocess2D(x,'Preprocessing',prepx);
[ycs,my,sdy] = preprocess2D(y,'Preprocessing',prepy);

model = simpls(xcs,ycs,'LVs',lvs);
ypred = (xcs*model.beta).*(ones(N,1)*sdy) + (ones(N,1)*my);

if ~isempty(test)
    testcs = preprocess2Dapp(test,m,'Scale',sd);
    testypred = (testcs*model.beta).*(ones(L,1)*sdy) + (ones(L,1)*my);
else
    testypred = [];
end

        