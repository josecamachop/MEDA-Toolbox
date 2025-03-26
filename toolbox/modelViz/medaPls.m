
function [map,ind,ord] = medaPls(x,y,varargin)

% Missing data methods for exploratory data analysis in PCA. The original
% paper is Chemometrics and Intelligent Laboratory Systems 103(1), 2010, pp.
% 8-18. 
%
% medaMap = medaPls(x,y) % minimum call
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
% LVs: [1xA] Latent Variables considered (e.g. pcs = 1:2 selects the
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
% 'Threshold': [1x1] threshold (0,1] for discretization and discarding (0.1 by default)
%
% 'Seriated': bool
%      false: no seriated (by default)
%      true: seriated
%
% 'Discard': bool
%      false: not discard
%      true: discard variables with diagonal below threshold      
%
% 'VarsLabel': [Mx1] name of the variables (numbers are used by default)
%
% 'Vars': [Sx1] Subset of variables to plot (1:M by default)
%
%
% OUTPUTS:
%
% map: [MxM] non-seriated MEDA matrix.
%
% ind: [M2x1] indices of seriated variables over the input threshold (M2 < M).
%
% ord: [M2x1] order of seriated variables.
%
%
%
% EXAMPLE OF USE: Seriation and discarding uninformative variables
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% lvs = 1;
% map = medaPls(X,Y,'LVs',lvs,'Threshold',0.3);
% map = medaPls(X,Y,'LVs',lvs,'Threshold',0.3,'Seriated',true);
% map = medaPls(X,Y,'LVs',lvs,'Threshold',0.3,'Seriated',true,'Discard',true);
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
addParameter(p,'Threshold',0.1);
addParameter(p,'Seriated',false);
addParameter(p,'Discard',false);  
addParameter(p,'VarsLabel',1:M);
addParameter(p,'Vars',1:M);      
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
thres = p.Results.Threshold;
seriated = p.Results.Seriated;
discard = p.Results.Discard;
label = p.Results.VarsLabel;
vars = p.Results.Vars;

% Convert row arrays to column arrays
if size(label,1) == 1, label = label'; end;
if size(vars,1)  == 1, vars = vars'; end;

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
assert (isequal(size(thres), [1 1]), 'Dimension Error: parameter ''Threshold'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: parameter ''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(size(vars) > [M 1])), 'Dimension Error: parameter ''Vars'' must be at most M-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (thres>0 && thres<=1, 'Value Error: parameter ''Threshold'' must be in (0,1]. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(vars<=0)) && isequal(fix(vars), vars) && isempty(find(vars>M)), 'Value Error: parameter ''Vars'' must contain positive integers below or equal to M. Type ''help %s'' for more info.', routine(1).name);


%% Main code

x = x(:,vars);
label = label(vars);

x2 = preprocess2D(x,'Preprocessing',prepx);
y2 = preprocess2D(y,'Preprocessing',prepy);

model = simpls(x2,y2,'LVs',lvs);
R = model.altweights;
P = model.loads;

map = meda(x2'*x2,R,'OutSubspace',P);

if discard
    Dmap = diag(map);
    ind = find(Dmap > thres);
else
    ind = 1:length(vars);
end
label = label(ind);
map = map(ind,ind);

if seriated
    [map, ord] = seriation(map);
else
    ord = 1:length(ind);
end
label = label(ord);
    

%% Show results

plotMap(map,'VarsLabel',label);
 