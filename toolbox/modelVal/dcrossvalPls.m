function [Qm,Q,lvso,keepXso] = dcrossvalPls(x,y,varargin)

% Row-wise k-fold (rkf) double cross-validation in PLS. The algorithm uses 
% repetitions of the dCV loop to estimate the stability: see Szymanska, E.,
% Saccenti, E., Smilde, A.K., Weterhuis, J. Metabolomics (2012) 8: 3. The 
% algorithm also considers the best trade-off when non-significant 
% differences among solutions in J. Camacho, J. González-Martínez and E. 
% Saccenti. Rethinking cross-validation in PLS with Variable Selection. 
% Submitted to Journal of Chemometrics.
%
% Qm = dcrossvalPls(x,y) % minimum call
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
%   first two LVs). By default, lvs = 0:rank(x)
%
% 'VarNumber': [1xK] Numbers of x-block variables selected. By default, VarNumber = M
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
% 'Repetitions': [1x1] number of repetitions for stability when 'MaxBlock' < N
%
% 'Plot': (bool) plot results if Repetition > 1
%       false: no plots.
%       true: plot (default)
%
% 'Selection': str
%   'Weights': filter method based on the PLS weights (W)
%   'AltWeights': filter method based on the PLS alternative weights (R)
%   'Regressors': filter method based on the PLS regression coefficients (beta)
%   'SR': filter method based on the selectivity ratio (by default)
%   'VIP': filter method based on Variance Importance in PLS Projection
%   'T2': wrapper method based on the Hotelling T2 statistic
%   'sPLS': embedded method based on sparse PLS
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
% keepXso: [blocksrx1] optimum number of keepXs in the inner loop
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,'LevelCorr',8);
% X = [X 0.1*randn(20,10) + X];
% Y = 0.1*randn(20,2) + X(:,1:2);
% lvs = 0:10;
% keepXs = 1:10;
% [Qm,Q,lvso,keepX] = dcrossvalPls(X,Y,'LVs',lvs,'VarNumber',keepXs,'MaxBlock',5)
% [Qmsimple,Qsimple,lvsosimple,keepXsimple] = dcrossvalPls(X,Y,'LVs',lvs,'VarNumber',keepXs,'Alpha',0.5,'MaxBlock',5)
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 31/Jan/2025
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
addParameter(p,'VarNumber',M);
addParameter(p,'Alpha',0);
addParameter(p,'MaxBlock',N);
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);
addParameter(p,'Repetitions',10);
addParameter(p,'Selection','SR'); 
addParameter(p,'Plot',true);    
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;
alpha = p.Results.Alpha;
keepXs = p.Results.VarNumber;
blocksr = p.Results.MaxBlock;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
rep = p.Results.Repetitions;
selection = p.Results.Selection;
opt = p.Results.Plot;

% Set repetitions to 1 for leave-one-out
if blocksr==N
    disp('Repetitions set to 1 for leave-one-out...')
    rep = 1;
end

% Extract LVs and Gamma length
A = length(lvs);
J =  length(keepXs);

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;
if size(keepXs,2) == 1, keepXs = keepXs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N O]), 'Dimension Error: parameter ''y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(keepXs), [1 J]), 'Dimension Error: parameter ''VarNumber'' must be 1-by-J. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(alpha), [1 1]), 'Dimension Error: parameter ''Alpha'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksr), [1 1]), 'Dimension Error: parameter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(rep), [1 1]), 'Dimension Error: parameter ''Repetition'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);
keepXs = unique(keepXs);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LVs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(keepXs), keepXs), 'Value Error: parameter ''VarNumber'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (alpha>=-1 & alpha<=1, 'Value Error: parameter ''Alpha'' must contain values in [-1, 1]. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocksr), blocksr), 'Value Error: parameter ''MaxBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr>3, 'Value Error: parameter ''MaxBlock'' must be above 3. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr<=N, 'Value Error: parameter ''MaxBlock'' must be at most N. Type ''help %s'' for more info.', routine(1).name);


%% Main code

for j=1:rep
    % Cross-validation
    
    rows = rand(1,N);
    [a,rind]=sort(rows);
    elemr=N/blocksr;
    
    for i=1:blocksr
        % disp(sprintf('Crossvalidation block %i of %i',i,blocksr))
        indi = rind(round((i-1)*elemr+1):round(i*elemr)); % Sample selection
        i2 = ones(N,1);
        i2(indi)=0;
        val = x(indi,:);
        rest = x(find(i2),:);
        valy = y(indi,:);
        resty = y(find(i2),:);
        
        [ccs,av,st] = preprocess2D(rest,'Preprocessing',prepx);
        [ccsy,avy,sty] = preprocess2D(resty,'Preprocessing',prepy);
        
        vcs = preprocess2Dapp(val,av,'Scale',st);
        vcsy = preprocess2Dapp(valy,avy,'Scale',sty);
        
        [cumpress,~,nze] =  crossvalPls(rest,resty,'LVs',lvs,'VarNumber',keepXs,'MaxBlock',blocksr-1,'PreprocessingX',prepx,'PreprocessingY',prepy,'Plot',false,'Selection',selection);
        
        cumpressb = (1-abs(alpha))*cumpress/max(max(cumpress)) + alpha*nze/max(max(nze));
        
        [l,k]=find(cumpressb==min(min(cumpressb)));
        lvso(j,i) = lvs(l(1));
        keepXso(j,i) = keepXs(k(1));
        
        if lvso(j,i)~=0
            model = vpls(ccs,ccsy,'LVs',1:lvso(j,i),'VarNumber',keepXso(j,i),'Selection',selection);
            srec = vcs*model.beta;
        else
            keepXso(j,i) = nan;
            srec = zeros(size(vcsy));
        end
        
        Qu(i) = sum(sum((vcsy-srec).^2));
        Qd(i) = sum(sum(vcsy.^2));
    end
    
    Q(j) = 1-sum(Qu)/sum(Qd);
    
end

Qm = mean(Q);


%% Show results if repetitions a larger than one

if opt && rep > 1
    plotVec(Q,'XYLabel',{'#Repetition','Goodness of Prediction'},'Plot','Lines'); 
end

