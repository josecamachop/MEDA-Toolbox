function [AUCm,AUC,lvso,keepXso] = dcrossvalPlsDA(x,y,varargin)

% Row-wise k-fold (rkf) double cross-validation in PLS-DA, restricted to 
% one response categorical variable of two levels. The algorithm uses 
% repetitions of the dCV loop to estimate the stability: see Szymanska, E.,
% Saccenti, E., Smilde, A.K., Weterhuis, J. Metabolomics (2012) 8: 3. It 
% also corrects the classification limit following Richard G. Brereton, J. 
% Chemometrics 2014; 28: 213–225. The algorithm also considers the best
% trade-off when non-significant differences among solutions. in J. Camacho, 
% J. González-Martínez and E. Saccenti. Rethinking cross-validation in PLS 
% with Variable Selection. Submitted to Journal of Chemometrics.
%
% AUCm = dcrossvalPlsDA(x,y) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [Nx1] billinear data set of one categorical variable with two levels.
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
% 'MaxBlock': [1x1] maximum number of blocks of samples (the minimum number
%   of observations of a class by default)
%
% 'PreprocesingX': [1x1] preprocesing of the x-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)  
%
% 'PreprocesingY': [1x1] preprocesing of the y-block
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
% AUCm: [1x1] Mean Area Under the ROC 
%
% AUC: [rep x 1] Area Under the ROC
%
% lvso: [rep x blocksr] optimum number of LVs in the inner loop
%
% keepXso: [rep x blocksr] optimum number of keepXs in the inner loop
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,'LevelCorr',8);
% X = [X 0.1*randn(20,10) + X];
% Y = 2*(0.1*randn(20,1) + X(:,1)>0)-1;
% lvs = 0:10;
% keepXs = 1:10;
% [AUCm,AUC,lvso,keepX] = dcrossvalPlsDA(X,Y,'LVs',lvs,'VarNumber',keepXs,'MaxBlock',5)
% [AUCmsimple,AUCsimple,lvsosimple,keepXsimple] = dcrossvalPlsDA(X,Y,'LVs',lvs,'VarNumber',keepXs,'Alpha',0.5,'MaxBlock',5)
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
O = size(y, 2);
M = size(x, 2);

vals = unique(y);
rep2 = sort(histc(y,vals),'descend');
N2 = rep2(2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
lat=0:rank(x);
addParameter(p,'LVs',lat'); 
addParameter(p,'VarNumber',M);
addParameter(p,'Alpha',0);
addParameter(p,'MaxBlock', max(3,round(N2/2)));
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
assert (isequal(size(y), [N 1]), 'Dimension Error: parameter ''y'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
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

assert (isempty(find(y~=1 & y~=-1)), 'Value Error: parameter ''y'' must not contain values different to 1 or -1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LVs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(keepXs), keepXs), 'Value Error: parameter ''VarNumber'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (alpha>=-1 & alpha<=1, 'Value Error: parameter ''Alpha'' must contain values in [-1, 1]. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksr), [1 1]), 'Dimension Error: parameter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocksr), blocksr), 'Value Error: parameter ''MaxBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr>3, 'Value Error: parameter ''MaxBlock'' must be above 3. Type ''help %s'' for more info.', routine(1).name);
assert (blocksr<=N, 'Value Error: parameter ''MaxBlock'' must be at most N. Type ''help %s'' for more info.', routine(1).name);


%% Main code

for j=1:rep
    % Cross-validation
    
    y1 = find(y==1);
    yn1 = find(y==-1);
    
    rows = rand(1,length(y1));
    [a,rind1]=sort(rows);
    elemr1=length(y1)/blocksr;
    
    rows = rand(1,length(yn1));
    [a,rindn1]=sort(rows);
    elemrn1=length(yn1)/blocksr;
    
    % Cross-validation
    
    for i=1:blocksr
        
        indi1 = rind1(round((i-1)*elemr1+1):round(i*elemr1)); % Sample selection
        i2 = ones(length(y1),1);
        i2(indi1)=0;
        val = x(y1(indi1),:);
        rest = x(y1(find(i2)),:);
        valy = y(y1(indi1),:);
        resty = y(y1(find(i2)),:);
        
        indin1 = rindn1(round((i-1)*elemrn1+1):round(i*elemrn1)); % Sample selection
        i2 = ones(length(yn1),1);
        i2(indin1)=0;
        val = [val;x(yn1(indin1),:)];
        rest = [rest;x(yn1(find(i2)),:)];
        valy = [valy;y(yn1(indin1),:)];
        resty = [resty;y(yn1(find(i2)),:)];
        
        [ccs,av,st] = preprocess2D(rest,'Preprocessing',prepx);
        %[ccsy,avy,sty] = preprocess2D(resty,prepy);
        ccsy = resty;
        
        [kk,m1] = preprocess2D(ccs(find(resty==1),:),'Preprocessing',1);  % additional subtraction of class mean
        [kk,mn1] = preprocess2D(ccs(find(resty==-1),:),'Preprocessing',1);
        ccs = preprocess2Dapp(ccs,(m1+mn1)/2);
        
        vcs = preprocess2Dapp(val,av,'Scale',st);
        vcs = preprocess2Dapp(vcs,(m1+mn1)/2);
        
        %vcsy = preprocess2Dapp(valy,avy,'Scale',sty);
        vcsy = valy;
        
        [AUCt,nze] =  crossvalPlsDA(rest,resty,'LVs',lvs,'VarNumber',keepXs,'MaxBlock',blocksr-1,'PreprocessingX',prepx,'PreprocessingY',prepy,'Plot',false,'Selection',selection);
        
        cumpressb = (abs(alpha)-1)*AUCt/max(max(AUCt)) + alpha*nze/max(max(nze));
        
        [l,k]=find(cumpressb==min(min(cumpressb)));
        lvso(j,i) = lvs(l(1));
        keepXso(j,i) = keepXs(k(1));
        
        if lvso(j,i)~=0
            model = vpls(ccs,ccsy,'LVs',1:lvso(j,i),'VarNumber',keepXso(j,i),'Selection',selection);
            sr = vcs*model.beta;
            srec1(indi1') = sr(1:length(indi1));
            srecn1(indin1') = sr(length(indi1)+1:end);
            
        else
            keepXso(j,i) = nan;
            srec1(indi1') = 0;
            srecn1(indin1') = 0;
        end
        
    end
    
    [~,~,~,AUC(j)] = perfcurve(y([y1;yn1]),[srec1';srecn1'],1);
end

AUCm = mean(AUC);


%% Show results if repetitions a larger than one

if opt && rep > 1
    plotVec(AUC,'XYLabel',{'#Repetition','AUC'},'Plot','Lines');
end


