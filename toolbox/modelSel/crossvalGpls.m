function [cumpress,press,nze] = crossvalGpls(x,y,varargin)
% Row-wise k-fold (rkf) cross-validation for square-prediction-errors
% computing in GPLS. 
%
% [cumpress,press] = crossvalGpls(x,y) % minimum call
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
% 'Gamma': [1xJ] gamma values considered. By default, gammas = 0:0.1:1
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
% cumpress: [A x gammas x 1] Cumulative PRESS
%
% press: [A x gammas x O] PRESS per variable.
%
% nze: [A x gammas x 1] Average number of non-zero elements
%
%
% EXAMPLE OF USE: Random data, two examples of use.
%
% obs = 20;
% vars = 100;
% X = simuleMV(obs,vars,'LevelCorr',5);
% X = [0.1*randn(obs,5)+X(:,1)*ones(1,5) X(:,6:end)];
% Y = sum((X(:,1:5)),2);
% Y = 0.1*randn(obs,1)*std(Y) + Y;
% gammas=[0 0.5:0.1:1];
% lvs = 0:10;
% 
% % Mean Centering example with default gammas
% [cumpress,press,nze] = crossvalGpls(X,Y,'LVs',lvs,'PreprocessingX',1,'PreprocessingY',1);
% 
% % Auto scaling example with gammas
% [cumpress,press,nze] = crossvalGpls(X,Y,'LVs',lvs,'Gamma',gammas);
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
addParameter(p,'MaxBlock',N);
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);
addParameter(p,'Plot',true);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility

lvs = p.Results.LVs;
gammas = p.Results.Gamma;
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
assert (isequal(size(blocksr), [1 1]), 'Dimension Error: paramter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);
gammas = unique(gammas);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LVs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(gammas<0 | gammas>1)), 'Value Error: parameter ''Gamma'' must not contain values out of [0,1]. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocksr), blocksr), 'Value Error: parameter ''MaxBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr>2, 'Value Error: parameter ''MaxBlock'' must be above 2. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr<=N, 'Value Error: parameter ''MaxBlock'' must be at most N. Type ''help %s'' for more info.', routine(1).name);


%% Main code

% Initialization
cumpress = zeros(length(lvs),length(gammas),1);
press = zeros(length(lvs),length(gammas),O);
nze = zeros(length(lvs),length(gammas));

rows = rand(1,N);
[a,rind]=sort(rows);
elemr=N/blocksr;

% Cross-validation
        
for i=1:blocksr
    
    indi = rind(round((i-1)*elemr+1):round(i*elemr)); % Sample selection
    i2 = ones(N,1);
    i2(indi)=0;
    sample = x(indi,:);
    calibr = x(find(i2),:); 
    sampley = y(indi,:);
    calibry = y(find(i2),:); 

    [ccs,av,st] = preprocess2D(calibr,'Preprocessing',prepx);
    [ccsy,avy,sty] = preprocess2D(calibry,'Preprocessing',prepy);
        
    scs = preprocess2Dapp(sample,av,'Scale',st);
    scsy = preprocess2Dapp(sampley,avy,'Scale',sty);
     
    gammas2 = gammas;
    gammas2(find(gammas==0)) = [];  
    if ~isempty(gammas2)
        [kk,kk,kk,kk,kk,kk,stree] = gplsMeda(ccs,ccsy,'LVs',1:max(lvs),'Gamma',min(gammas2));
    else
        stree = [];
    end

    for gamma=1:length(gammas)
        
        [beta,W,P,Q,R] = gplsMeda(ccs,ccsy,'LVs',1:max(lvs),'Gamma',gammas(gamma),'Stree',stree);
            
        for lv=1:length(lvs)
                
            if lvs(lv) > 0
                Q2 = Q(:,1:min(lvs(lv),size(Q,2)));
                R2 = R(:,1:min(lvs(lv),size(Q,2)));
                beta=R2*Q2';
                srec = scs*beta;
                
                pem = scsy-srec;
                nze(lv,gamma) = nze(lv,gamma) + length(find(beta));
                            
            else % Modelling with the average
                pem = scsy;
            end
                
            press(lv,gamma,:) = squeeze(press(lv,gamma,:))' + sum(pem.^2,1);

        end
    end
    
end

cumpress = sum(press,3);

%% Show results

if opt
    plotVec(cumpress','EleLabel',gammas,'XYLabel',{'\gamma','PRESS'},'PlotType','Lines','VecLabel',lvs); 
end

