function [AUCm,AUC,lvso] = dcrossval_pls_da(x,y,lvs,blocks_r,prepx,prepy,opt)

% Row-wise k-fold (rkf) double cross-validation in PLS-DA, restricted to one 
% response categorical variable of two levels. 
%
% AUCm = dcrossval_pls_da(x,y) % minimum call
% [AUCm,AUC,lvso] = dcrossval_spls_da(x,y,lvs,blocks_r,prepx,prepy,opt) % complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [Nx1] billinear data set of one categorical variable with two levels
%
% lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x)
%
% blocks_r: [1x1] maximum number of blocks of samples (the minimum number
%   of observations of a class divided by 2 by default)
%
% prepx: [1x1] preprocesing of the x-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)
%
% prepy: [1x1] preprocesing of the y-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)
%
% opt: [1x1] options for data plotting
%       0: no plots
%       1: bar plot (default)
%
%
% OUTPUTS:
%
% AUCm: [1x1] Mean Area Under the ROC 
%
% AUC: [blocks_rx1] Area Under the ROC
%
% lvso: [blocks_rx1] optimum number of LVs in the inner loop
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,8);
% Y = 2*(0.1*randn(20,1) + X(:,1)>0)-1;
% lvs = 0:10;
% [AUCm,AUC,lvso] = dcrossval_pls_da(X,Y,lvs,5)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Edoardo Saccenti (edoardo.saccenti@wur.nl )
% last modification: 04/Apr/18.
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez, Edoardo Saccenti
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
if nargin < 3 || isempty(lvs), lvs = 0:rank(x); end;
A = length(lvs);

vals = unique(y);
rep = sort(histc(y,vals),'descend');
N2 = rep(2);
if nargin < 4 || isempty(blocks_r), blocks_r = max(3,round(N2/2)); end;
if nargin < 5 || isempty(prepx), prepx = 2; end;
if nargin < 6 || isempty(prepy), prepy = 2; end;
if nargin < 7 || isempty(opt), opt = 1; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N 1]), 'Dimension Error: 2nd argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocks_r), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);

% Validate values of input data
assert (isempty(find(y~=1 & y~=-1)), 'Value Error: 2rd argument must not contain values different to 1 or -1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(lvs<0)), 'Value Error: 3rd argument must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: 3rd argumentmust contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocks_r), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocks_r), blocks_r), 'Value Error: 4th argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r>3, 'Value Error: 4th argument must be above 3. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r<=N, 'Value Error: 4th argument must be at most N. Type ''help %s'' for more info.', routine(1).name);


%% Main code

% Cross-validation

y1 = find(y==1);
yn1 = find(y==-1);

rows = rand(1,length(y1));
[a,r_ind1]=sort(rows);
elem_r1=length(y1)/blocks_r;

rows = rand(1,length(yn1));
[a,r_indn1]=sort(rows);
elem_rn1=length(yn1)/blocks_r;

% Cross-validation

for i=1:blocks_r,
    
    ind_i = r_ind1(round((i-1)*elem_r1+1):round(i*elem_r1)); % Sample selection
    i2 = ones(length(y1),1);
    i2(ind_i)=0;
    val = x(y1(ind_i),:);
    rest = x(y1(find(i2)),:);
    val_y = y(y1(ind_i),:);
    rest_y = y(y1(find(i2)),:);
    
    ind_i = r_indn1(round((i-1)*elem_rn1+1):round(i*elem_rn1)); % Sample selection
    i2 = ones(length(yn1),1);
    i2(ind_i)=0;
    val = [val;x(yn1(ind_i),:)];
    rest = [rest;x(yn1(find(i2)),:)];
    val_y = [val_y;y(yn1(ind_i),:)];
    rest_y = [rest_y;y(yn1(find(i2)),:)];  ;
    
    [ccs,av,st] = preprocess2D(rest,prepx);
    [ccs_y,av_y,st_y] = preprocess2D(rest_y,prepy);
    
    vcs = preprocess2Dapp(val,av,st);
    vcs_y = preprocess2Dapp(val_y,av_y,st_y);
    
    AUCt =  crossval_pls_da(rest,rest_y,lvs,blocks_r-1,prepx,prepy,0);
    
    idx=find(AUCt==max(AUCt),1);
    lvso(i) = lvs(idx);
    
    if lvso(i)~=0,
        
        X = ccs;
        Y = ccs_y;
        XY = X'*Y;
        XX = X'*X;
        
        beta = kernel_pls(XX,XY,1:lvso(i));
        
        srec = vcs*beta;
        [~,~,~,AUC(i)] = perfcurve(val_y,srec,1);
        
    else
        AUC(i) = [];
    end
    
end

AUCm = mean(AUC);

%% Show results

if opt == 1,
    fig_h = plot_vec(AUC,[],[],{'#Split','AUC'},[],1);
end

