function [AUC,nze] = crossvalPlsDA(x,y,varargin)

% Row-wise k-fold (rkf) cross-validation in PLS-DA. We correct the 
% classification limit following Richard G. Brereton, J. Chemometrics 2014; 
% 28: 213–225. We extend to several classes by counting positives/negatives
% in each response dummy variable independently. 
%
% AUC = crossvalPlsDA(x,y) % minimum call
%
%
% See also: crossvalPls, crossvalPca
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of dummy variables (+1, -1)
%
%
% Optional INPUTS (parameter):
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x)
%
% 'VarNumber': [1xK] Numbers of x-block variables selected. By default, VarNumber = M
%
% 'MaxBlock': [1x1] maximum number of blocks of samples (the minimum number
%   of observations of a class divided by 2 by default)
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
% AAUC: [AxK] Macro-average Area Under the Curve in ROC
%
% AUC: [AxKxO] Area Under the Curve in ROC per variable
%
% nze: [AxK] Non-zero elements in the regression coefficient matrix.
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 2*(0.1*randn(20,1) + X(:,1)>0)-1;
% lvs = 0:10;
% keepXs = 1:10;
% [AUC,nze] = crossvalPlsDA(X,Y,'LVs',lvs,'MaxBlock',5);
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

ind = (size(y,2)+1)*ones(size(y,1),1);
[r,c]=find(y==1);
[r1,r2]=sort(r);
ind(r1) = c(r2);
vals = unique(ind);
rep = sort(histc(ind,vals),'ascend');
N2 = rep(1); % minimum length of a class


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
lat=0:rank(x);
addParameter(p,'LVs',lat'); 
addParameter(p,'VarNumber',M);
addParameter(p,'MaxBlock',2);
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);
addParameter(p,'Selection','SR'); 
addParameter(p,'Plot',true);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;
keepXs = p.Results.VarNumber;
blocksR = p.Results.MaxBlock;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
selection = p.Results.Selection;
opt = p.Results.Plot;

% Extract LVs and VarNumber length
A = length(lvs);
J =  length(keepXs);

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;
if size(keepXs,2) == 1, keepXs = keepXs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N O]), 'Dimension Error: parameter ''y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(keepXs), [1 J]), 'Dimension Error: parameter ''VarNumber'' must be 1-by-J. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksR), [1 1]), 'Dimension Error: parameter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);
keepXs = unique(keepXs);

% Validate values of input data
assert (isempty(find(y~=1 & y~=-1)), 'Value Error: parameter ''y'' must not contain values different to 1 or -1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LVs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(keepXs), keepXs), 'Value Error: parameter ''VarNumber'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksR), [1 1]), 'Dimension Error: parameter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocksR), blocksR), 'Value Error: parameter ''MaxBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocksR>1, 'Value Error: parameter ''MaxBlock'' must be above 1. Type ''help %s'' for more info.', routine(1).name);
assert (blocksR<=N2, 'Value Error: parameter ''MaxBlock'' must be at most %d. Type ''help %s'' for more info.', N2, routine(1).name);

%% Main code

% Initialization
AAUC = zeros(length(lvs),length(keepXs));
AUC = zeros(length(lvs),length(keepXs),O);
nze = zeros(length(lvs),length(keepXs));

ind = (size(y,2)+1)*ones(size(y,1),1);
[r,c]=find(y==1);
[r1,r2]=sort(r);
ind(r1) = c(r2);
vals = unique(ind);
for i=1:length(vals)
    y1{i} = find(ind==vals(i));
    rows = rand(1,length(y1{i}));
    [a,rIndn{i}]=sort(rows);
    elemR(i)=length(y1{i})/blocksR;
end

% Cross-validation
        
for i=1:blocksR
    
    cal = [];
    test = [];
    for j=1:length(vals)
        indIn1 = rIndn{j}(round((i-1)*elemR(j)+1):round(i*elemR(j))); % Sample selection
        i2 = ones(length(y1{j}),1);
        i2(indIn1)=0;
        cal = [cal;y1{j}(find(i2))];
        test = [test;y1{j}(indIn1)];
    end
    sample = x(test,:);
    calibr = x(cal,:);
    sampleY = y(test,:);
    calibrY = y(cal,:);  
    
    [ccs,av,st] = preprocess2D(calibr,'Preprocessing',prepx);    
    %[ccsY,avy,sty] = preprocess2D(calibrY,prepy);
    ccsY = calibrY;
    
    ind = (size(ccsY,2)+1)*ones(size(ccsY,1),1);
    [r,c]=find(ccsY==1);
    [r1,r2]=sort(r);
    ind(r1) = c(r2);
    vals = unique(ind);    
    for j=1:length(vals)
        ind2 = find(ind==vals(j));
        if ~isempty(ind2)
            [kk,m(j,:)] = preprocess2D(ccs(ind2,:),'Preprocessing',1);  % additional subtraction of class mean
        end
    end
    ccs = preprocess2Dapp(ccs,mean(m));
        
    scs = preprocess2Dapp(sample,av,'Scale',st);
    scs = preprocess2Dapp(scs,mean(m));

    %[ccs,PR] = reduce2(ccs,coef);
    
    if  ~isempty(find(lvs))
        
        for lv=1:length(lvs)

            for keepX=1:length(keepXs)
                
                if lvs(lv)
                    model = vpls(ccs,ccsY,'LVs',1:lvs(lv),'VarNumber',keepXs(keepX),'Selection',selection);
                   
                    srec1(test,lv,keepX,:) = scs*model.beta;%scs*PR*beta;
					nze(lv,keepX) = nze(lv,keepX) + length(find(model.beta)); 
                else
                    srec1(test,lv,keepX,:) = 0;
					nze(lv,keepX) = nze(lv,keepX) + M*O; 
                end
                
            end
            
        end
        
    else
        srec1(test,1,1,:) = 0;
		nze = nze + ones(length(keepXs),1)*M*O;
    end
    
end

for lv=1:size(srec1,2)
    for keepX=1:size(srec1,3)
        for o = 1:O
            [~,~,~,AUC(lv,keepX,o)] = perfcurve(y(:,o),srec1(:,lv,keepX,o),1);
        end
    end
end

AAUC =  mean(AUC,3);


%% Show results

if opt
    plotVec(AAUC','EleLabel',keepXs,'XYLabel',{'#NZV','AUC'},'PlotType','Lines','VecLabel',lvs,'Color','jet'); 
end

