function [Qm,Q,lvso,gammaso] = dcrossvalGpls(x,y,varargin)

% Row-wise k-fold (rkf) double cross-validation for square-prediction-errors computing in GPLS.
%
% Qm = dcrossvalGpls(x,y) % minimum call
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
%   first two LVs). By default, lvs = 0:rank(x)
%
% 'Gamma': [1xJ] gamma values considered. By default, gammas = 0:0.1:1
%
% 'Alpha': [1x1] Trade-off controlling parameter that goes from -1 (maximum 
%   completeness), through 0 (pure prediction, by default) to 1 (maximum 
%   parsimony) 
%
% 'MaxBlock': [1x1] maximum number of blocks of samples (N by default)
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
% 'Plot': (bool) plot results
%       false: no plots.
%       true: plot (default)
%
%
% OUTPUTS:
%
% Qm: [1x1] Mean Goodness of Prediction
%
% Q: [blocksrx1] Goodness of Prediction
%
% lvso: [blocksrx1] optimum number of LVs in the inner loop
%
% gammaso: [blocksrx1] optimum gamma in the inner loop
%
%
% EXAMPLE OF USE: Random data with structural relationship
% 
% obs = 20;
% vars = 100;
% X = simuleMV(obs,vars,'LevelCorr',5);
% X = [0.1*randn(obs,5)+X(:,1)*ones(1,5) X(:,6:end)];
% Y = sum((X(:,1:5)),2);
% Y = 0.1*randn(obs,1)*std(Y) + Y;
% 
% lvs = 0:10;
% gammas = [0 0.5:0.1:1];
% [Qm,Q,lvso,gammaso] = dcrossvalGpls(X,Y,'LVs',lvs,'Gamma',gammas,'MaxBlock',5)
% [Qmsimple,Qsimple,lvsosimple,gammasosimple] = dcrossvalGpls(X,Y,'LVs',lvs,'Gamma',gammas,'Alpha',0.5,'MaxBlock',5)
% 

% coded by: Jose Camacho (josecamacho@ugr.es)
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
O = size(y, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
lat=0:rank(x);
addParameter(p,'LVs',lat'); 
gam = 0:0.1:1;
addParameter(p,'Gamma',gam);
addParameter(p,'Alpha',0);
addParameter(p,'MaxBlock',N);
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);
addParameter(p,'Plot',true);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility

lvs = p.Results.LVs;
gammas = p.Results.Gamma;
alpha = p.Results.Alpha;
blocksr = p.Results.MaxBlock;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
opt = p.Results.Plot;

% Extract LVs and Gamma length
A = length(lvs);
J =  length(gammas);

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;
if size(gammas,2) == 1, gammas = gammas'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N O]), 'Dimension Error: parameter ''y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(gammas), [1 J]), 'Dimension Error: parameter ''Gamma'' must be 1-by-J. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(alpha), [1 1]), 'Dimension Error: parameter ''Alpha'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksr), [1 1]), 'Dimension Error: parameter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);
gammas = unique(gammas);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LVs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(gammas<0 | gammas>1)), 'Value Error: parameter ''Gamma'' must not contain values out of [0,1]. Type ''help %s'' for more info.', routine(1).name);
assert (alpha>=-1 & alpha<=1, 'Value Error: parameter ''Alpha'' must not be out of [1,1]. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocksr), blocksr), 'Value Error: parameter ''MaxBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr>3, 'Value Error: parameter ''MaxBlock'' must be above 3. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr<=N, 'Value Error: parameter ''MaxBlock'' must be at most N. Type ''help %s'' for more info.', routine(1).name);


%% Main code

% Cross-validation

rows = rand(1,N);
[a,rind]=sort(rows);
elemr=N/blocksr;
        
for i=1:blocksr
    disp(sprintf('Crossvalidation block %i of %i',i,blocksr))
    indi = rind(round((i-1)*elemr+1):round(i*elemr)); % Sample selection
    i2 = ones(N,1);
    i2(indi)=0;
    val = x(indi,:);
    rest = x(find(i2),:); 
    valy = y(indi,:);
    resty = y(find(i2),:);
        
    [cumpress,kk,nze] =  crossvalGpls(rest,resty,'LVs',lvs,'Gamma',gammas,'Maxblock',blocksr-1,'PreprocessingX',prepx,'PreprocessingY',prepy,'Plot',false);
       
    cumpressb = (1-abs(alpha))*cumpress/max(max(cumpress)) + alpha*nze/max(max(nze));
    
    [l,g]=find(cumpressb==min(min(cumpressb)));
    lvso(i) = lvs(l(1));
    gammaso(i) = gammas(g(1));
    
    [ccs,av,st] = preprocess2D(rest,'Preprocessing',prepx);
    [ccsy,avy,sty] = preprocess2D(resty,'Preprocessing',prepy);
    
    vcs = preprocess2Dapp(val,av,'Scale',st);
    vcsy = preprocess2Dapp(valy,avy,'Scale',sty);
    
    beta = gplsMeda(ccs,ccsy,'LVs',1:lvso(i),'Gamma',gammaso(i));
    srec = vcs*beta;
    
    Q(i) = 1 - sum(sum((vcsy-srec).^2))/sum(sum(vcsy.^2));
    
end

Qm = mean(Q);

%% Show results

if opt
    plotVec(Q,'XYLabel',{'#Split','Goodness of Prediction'},'Plot','Lines'); 
end

