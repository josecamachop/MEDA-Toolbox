
function [omedaVec,lim] = omedaPca(x,pcs,test,dummy,varargin)

% Observation-based Missing data methods for Exploratory Data Analysis 
% (oMEDA) for PCA. The original paper is Journal of Chemometrics, 2011, 25 
% (11): 592-600. This algorithm follows the direct computation for
% Known Data Regression (KDR) missing data imputation.
%
% omedaVec = omedaPca(x,pcs,test,dummy) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 1:rank(xcs)
%
% test: [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% dummy: [Lx1] dummy variable containing weights for the observations to 
%   compare, and 0 for the rest of observations. By default all test
%   observations are set to 1.
%
%
% Optional INPUTS (parameters):
%
% 'Preprocessing': [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
%
% 'ControlLim': bool
%       false: no control limits (by default)
%       true: plot control limits 
%
% 'Normalize': bool
%       false: no normalization (by default)
%       true: normalize by control limits
%
% 'VarsLabel': [Mx1] name of the variables (numbers are used by default)
%
% 'VarsClass': [Mx1] groups of variables (one group by default)
%
%
% OUTPUTS:
%
% omedaVec: [Mx1] oMEDA vector.
%
% lim: [Mx1] oMEDA limits.
%
%
% EXAMPLE OF USE: Anomaly on first observation and first 2 variables.
% 
% nobs = 100;
% nvars = 10;
% nPCs = 10;
% X = simuleMV(nobs,nvars,'LevelCorr',6);
% 
% nobst = 10;
% test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',cov(X)*(nobst-1));
% test(1,1:2) = 10*max(abs(X(:,1:2))); 
% dummy = zeros(10,1);
% dummy(1) = 1;
% 
% pcs = 1:nPCs;
% 
% omedaVec = omedaPca(X,pcs,test,dummy,'Preprocessing',1);
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
if isempty(pcs), pcs = 1:rank(x); end;
if isempty(test), test = x; end;
L = size(test, 1);
if nargin < 4 || isempty(dummy), dummy = ones(L,1); end;


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Preprocessing',2);     
addParameter(p,'ControlLim',false);    
addParameter(p,'Normalize',false);  
addParameter(p,'VarsLabel',1:M);  
addParameter(p,'VarsClass',ones(M,1));  
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
prep = p.Results.Preprocessing; 
ctrl = p.Results.ControlLim;
norm = p.Results.Normalize;
label = p.Results.VarsLabel;
classes = p.Results.VarsClass;

% Convert row arrays to column arrays
if size(label,1) == 1, label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
A = length(pcs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''pcs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''pcs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(test), [L M]), 'Dimension Error: parameter ''test'' must be L-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(dummy), [L 1]), 'Dimension Error: parameter ''dummy''must be L-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing''must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: parameter ''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(classes), [M 1]), 'Dimension Error: parameter ''VarsClass'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 

% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''pcs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[xcs,m,sd] = preprocess2D(x,'Preprocessing',prep);

model = pcaEig(xcs,'PCs',pcs);
P = model.loads;
    
testcs = preprocess2Dapp(test,m,'Scale',sd);
omedaVec = omeda(testcs,dummy,P);

% heuristic: 95% limit for one-observation-dummy
xr = xcs*P*P';
omedax = abs((2*xcs-xr).*(xr));
lim = prctile(omedax,95)';
    

%% Show results
    
vec = omedaVec;

if ctrl
    limp = lim;
else
    limp = [];
end

if norm
    ind = find(lim>1e-10);
    vec(ind) = vec(ind)./lim(ind);
    if ~isempty(limp)
        limp(ind) = limp(ind)./lim(ind);
    end
end

plotVec(vec,'EleLabel',label,'ObsClass',classes,'XYLabel',{[],'d^2_A'},'LimCont',[limp -limp]);
    

        