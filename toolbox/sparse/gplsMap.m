function  [beta,W,P,Q,R,bel,stree] = gplsMap(xcs,ycs,varargin)

% Group-wise Partial Least Squares based on correlation. This routine includes the 
% map estimation with MEDA, groups identification with GIA and model
% calibration with GPLS. The original paper is Camacho, J., 
% Saccenti, E. Group-wise Partial Least Squares Regression. Journal of 
% Chemometrics, 2017.
%
% beta = gplsMeda(xcs,ycs)     % minimum call
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% ycs: [NxO] preprocessed billinear data set of predicted variables
%
%
% Optional INPUTS (parameters):
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(xcs)
%
% 'Gamma': [1x1] correlation threshold to identify groups (0.7 by default)
%
% 'Stree': [struct] tree with GIA division, if previously executed. This
%   structure makes reiterative GIA computations faster (empty by default)
%   - tree: cell with division tree 
%   - indm: number of variables above a given threshold
%   - index: value for each tree division
%
% 'Type': 'Pearson' (by default) | 'Kendall' | 'Spearman' | 'MEDA'
%
% 'Map': 'XX' (by default) | 'XY' (only for MEDA)
%
%
% OUTPUTS:
%
% beta: [MxO] matrix of regression coefficients: W*inv(P'*W)*Q'
%
% W: [MxA] matrix of weights
%
% P: [MxA] matrix of x-loadings
%
% Q: [OxA] matrix of y-loadings
%
% R: [MxA] matrix of modified weights: W*inv(P'*W)
%
% bel: [Ax1] correspondence between LVs and States.
%
% stree: [struct] tree with GIA division, if previously executed. This
%   structure makes reiterative GIA computations faster.
%   - tree: cell with division tree 
%   - indm: number of variables above a given threshold
%   - index: value for each tree division
%
%
% EXAMPLE OF USE: Random data:
%
% obs = 20;
% vars = 100;
% X = simuleMV(obs,vars,'LevelCorr',5);
% X = [0.1*randn(obs,5)+X(:,1)*ones(1,5) X(:,6:end)];
% Y = sum((X(:,1:5)),2);
% Y = 0.1*randn(obs,1)*std(Y) + Y;
% lvs = 1;
% [beta,W,P,Q,R,bel] = gplsMap(X,Y,'LVs',lvs);
% 
% plotVec(beta,'XYLabel',{'','Regression coefficients'});
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 26/Mar/2025
% Dependencies: Matlab R2017b, MEDA v1.8
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine.name);
N = size(xcs, 1);
M = size(xcs, 2);
O = size(ycs, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
LVS = 0:rank(xcs);
addParameter(p,'LVs',LVS);  
addParameter(p,'Gamma',0.7);     
addParameter(p,'Stree',{});  
addParameter(p,'Type','Pearson');
addParameter(p,'Map','XX');
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;
gamma = p.Results.Gamma;
stree = p.Results.Stree;
type = p.Results.Type;
mapdata = p.Results.Map;


% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
lvs(find(lvs>M)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (isequal(size(ycs), [N O]), 'Dimension Error: parameter ''ycs'' must be N-by-O. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(gamma), [1 1]), 'Dimension Error: parameter ''Gamma'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain positive integers. Type ''help %s'' for more info.', routine.name);
assert (gamma>=0 && gamma<=1, 'Value Error: parameter ''Gamma'' must be between 0 and 1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

M=size(xcs,2);
O=size(ycs,2);

if isempty(stree) 

    if type=="MEDA"
        if mapdata=="XX"
            model = pcaEig(xcs,'PCs',lvs);
            map = meda(xcs'*xcs,model.loads);
        elseif mapdata=="XY"    
            model = simpls(xcs,ycs,'LVs',lvs);
            map = meda(xcs'*xcs,model.altweights,'OutSubspace',model.loads);
        end
    else
        map = corr(xcs,xcs,"Type",type);
    end

else
    map = zeros(M);
end

[bel,states,stree] = gia(map,'Gamma',gamma,'MinSize',1,'Stree',stree);

[beta,W,P,Q,R,bel] = gpls(xcs,ycs,states,'LVs',lvs);
