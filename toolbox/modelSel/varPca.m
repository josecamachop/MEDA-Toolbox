
function [xvar,cumpress] = varPca(x,varargin)

% Variability captured in terms of the number of PCs. It includes the ckf
% algorithm.
%
% varPca(x,pcs) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
%
% Optional INPUTS (parameters):
%
% 'Pcs': [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:rank(x). The value for 0 PCs is
%   added at the begining if not specified.
%
% 'Preprocessing': [1x1] preprocesing
%       0: no preprocessing 
%       1: mean-centering 
%       2: auto-scaling (default)  
%
% 'PlotCkf': bool
%       false: Residual Variance in X  
%       true: Residual Variance in X and ckf
%
%
% OUTPUTS:
%
% xvar: [Ax1] Percentage of captured variance of X.
%
% cumpress: [Ax1] ckf curve.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,'LevelCorr',8);
% pcs = 0:10;
% xvar = varPca(X,'PCs',pcs);
%
%
% codified by: Jose Camacho (josecamacho@ugr.es)
% last modification: 16/Jan/2025
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'PCs',0:rank(x));   
addParameter(p,'Preprocessing',2);   
addParameter(p,'PlotCkf',false);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
pcs = p.Results.PCs;
ckfplot = p.Results.PlotCkf;
prep = p.Results.Preprocessing;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
A = length(pcs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''Pcs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''Pcs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
pcs = unique([0 pcs]);

% Validate values of input data
assert (isempty(find(pcs<0)), 'Value Error: parameter ''Pcs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(pcs), pcs), 'Value Error: parameter ''Pcs'' contain integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,'Preprocessing',prep); 

model = pcaEig(xcs,'PCs',1:max(pcs));
P = model.loads;
T = model.scores;
pcs(find(pcs>size(P,2))) = [];

totalVx = sum(sum(xcs.^2));
xvar = ones(length(pcs),1);

for i = 1:length(pcs)
    xvar(i) = xvar(i) - sum(eig(T(:,1:pcs(i))'*T(:,1:pcs(i))))/totalVx;
end

cumpress = zeros(length(pcs),1);
if nargout>1 || ckfplot
    for i = 1:length(pcs)
         c = ckf(xcs,T(:,1:pcs(i)),P(:,1:pcs(i)),'Plot',false);
         cumpress(i) = c(end);
    end
end

    
%% Show results

if ckfplot
    plotVec([xvar cumpress/cumpress(1)],'PlotType','Lines','EleLabel',pcs,'XYLabel',{'#PCs','% Residual Variance'},'VecLabel',{'X','ckf'});
else
    plotVec(xvar,'PlotType','Lines','EleLabel',pcs,'XYLabel',{'#PCs','% Residual Variance'});
end

        