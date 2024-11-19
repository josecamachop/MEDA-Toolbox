function [AUCm,AUC,lvso] = dcrossvalPlsDA(x,y,varargin)

% Row-wise k-fold (rkf) double cross-validation in PLS-DA, restricted to one 
% response categorical variable of two levels. The algorithm uses
% repetitions of the dCV loop to estimate the stability: see Szymanska, E., 
% Saccenti, E., Smilde, A.K., Weterhuis, J. Metabolomics (2012) 8: 3. It
% also corrects the classification limit following Richard G. Brereton, 
% J. Chemometrics 2014; 28: 213–225
%
% AUCm = dcrossvalPlsDA(x,y) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [Nx1] billinear data set of one categorical variable with two levels
%
%
% Optional INPUTS (parameters):
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x)
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
% 'Repetition': [1x1] number of repetitinos for stability
%
% 'Option': [1x1] options for data plotting
%       0: no plots
%       1: bar plot (default)
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
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 2*(0.1*randn(20,1) + X(:,1)>0)-1;
% lvs = 0:10;
% [AUCm,AUC,lvso] = dcrossvalPlsDA(X,Y,'LVs',lvs,'Repetition',5)
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 22/Apr/24
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


%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
O = size(y, 2);
vals = unique(y);
repb = sort(histc(y,vals),'descend');
N2 = repb(2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
lat=0:rank(x);
blocksr = max(3,round(N2/2));
addParameter(p,'LVs',lat'); 
addParameter(p,'MaxBlock',blocksr);
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);
addParameter(p,'Repetitions',10);
addParameter(p,'Option',1);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility

lvs = p.Results.LVs;
blocksr = p.Results.MaxBlock;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
rep = p.Results.Repetitions;
opt = p.Results.Option;

% Extract LVs length
A = length(lvs);

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N 1]), 'Dimension Error: parameter ''y'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocksr), [1 1]), 'Dimension Error: parameter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(rep), [1 1]), 'Dimension Error: parameter ''Repetitions'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: parameter ''Option'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);

% Validate values of input data
assert (isempty(find(y~=1 & y~=-1)), 'Value Error: parameter ''y'' must not contain values different to 1 or -1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LVs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' contain integers. Type ''help %s'' for more info.', routine(1).name);
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
        
        %vcsy = preprocess2Dapp(valy,avy,Scale,sty);
        vcsy = valy;
        
        AUCt =  crossvalPlsDA(rest,resty,'LVs',lvs,'MaxBlock',blocksr-1,'PreprocessingX',prepx,'PreprocessingY',prepy,'Option',0);
        
        idx=find(AUCt==max(AUCt),1);
        lvso(j,i) = lvs(idx);
        
        if lvso(j,i)~=0
            
            X = ccs;
            Y = ccsy;
              
            model = simpls(X,Y,'LVs',1:lvso(j,i));
            
            sr = vcs*model.beta;
            srec1(indi1') = sr(1:length(indi1));
            srecn1(indin1') = sr(length(indi1)+1:end);
            
        else
             srec1(indi1') = 0;
             srecn1(indin1') = 0;
        end
        
    end
    
    [~,~,~,AUC(j)] = perfcurve(y([y1;yn1]),[srec1';srecn1'],1);
end

AUCm = mean(AUC);

%% Show results

if opt == 1
    figh = plotVec(AUC,'XYLabel',{'#Repetition','AUC'},'Option',11);
end

