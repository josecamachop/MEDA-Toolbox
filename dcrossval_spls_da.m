function [AUCm,AUC,lvso,keepXso] = dcrossval_spls_da(x,y,lvs,keepXs,alpha,blocks_r,prepx,prepy,opt)

% Row-wise k-fold (rkf) double cross-validation in SPLS-DA, restricted to 
% one response categorical variable of two levels. Reference:
% J. Camacho, J. González-Martínez and E. Saccenti.
% Rethinking cross-validation in SPLS. Submitted to Journal of Chemometrics. 
%
% AUCm = dcrossval_spls_da(x,y) % minimum call
% [AUCm,AUC,lvso,keepXso] = dcrossval_spls_da(x,y,lvs,keepXs,alpha,blocks_r,prepx,prepy,opt) % complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [Nx1] billinear data set of one categorical variable with two levels.
%
% lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x)
%
% keepXs: [1xK] Numbers of x-block variables kept per latent variable modeled. By default,
% 	keepXs = 1:M
%
% alpha: [1x1] Trade-off controlling parameter that goes from -1 (maximum 
%   completeness), through 0 (pure prediction, by default) to 1 (maximum 
%   parsimony) 
%
% blocks_r: [1x1] maximum number of blocks of samples (the minimum number
%   of observations of a class by default)
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
% keepXso: [blocks_rx1] optimum number of keepXs in the inner loop
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,8);
% X = [X 0.1*randn(20,10) + X];
% Y = 2*(0.1*randn(20,1) + X(:,1)>0)-1;
% lvs = 0:10;
% keepXs = 1:10;
% [AUCm,AUC,lvso,keepX] = dcrossval_spls_da(X,Y,lvs,keepXs,0,5)
% [AUCm_simple,AUC_simple,lvso_simple,keepX_simple] = dcrossval_spls_da(X,Y,lvs,keepXs,0.5,5)
% [AUCm_complete,AUC_complete,lvso_complete,keepX_complete] = dcrossval_spls_da(X,Y,lvs,keepXs,-0.5,5)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 04/Apr/18.
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez
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
if nargin < 4 || isempty(keepXs), keepXs = 1:M; end;
J =  length(keepXs);
if nargin < 5 || isempty(alpha), alpha = 0; end;

vals = unique(y);
rep = sort(histc(y,vals),'descend');
N2 = rep(2);
if nargin < 6 || isempty(blocks_r), blocks_r = max(3,round(N2/2)); end;
if nargin < 7 || isempty(prepx), prepx = 2; end;
if nargin < 8 || isempty(prepy), prepy = 2; end;
if nargin < 9 || isempty(opt), opt = 1; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;
if size(keepXs,2) == 1, keepXs = keepXs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N 1]), 'Dimension Error: 2nd argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(keepXs), [1 J]), 'Dimension Error: 4th argument must be 1-by-J. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(alpha), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocks_r), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 9th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);
keepXs = unique(keepXs);

% Validate values of input data

assert (isempty(find(y~=1 & y~=-1)), 'Value Error: 2rd argument must not contain values different to 1 or -1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(lvs<0)), 'Value Error: 3rd argument must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: 3rd argumentmust contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(keepXs), keepXs), 'Value Error: 4th argument must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (alpha>=-1 & alpha<=1, 'Value Error: 5th argument must contain values in [-1, 1]. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocks_r), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocks_r), blocks_r), 'Value Error: 6th argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r>3, 'Value Error: 6th argument must be above 3. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r<=N, 'Value Error: 6th argument must be at most N. Type ''help %s'' for more info.', routine(1).name);


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
    disp(sprintf('Crossvalidation block %i of %i',i,blocks_r))
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
        
    [AUCt,nze] =  crossval_spls_da(rest,rest_y,lvs,keepXs,blocks_r-1,prepx,prepy,0);
        
    cumpressb = (abs(alpha)-1)*AUCt/max(max(AUCt)) + alpha*nze/max(max(nze));
 
    [l,k]=find(cumpressb==min(min(cumpressb)));
    lvso(i) = lvs(l(1));
    keepXso(i) = keepXs(k(1));
    
    if lvso(i)~=0,
        
        model = sparsepls2(ccs, ccs_y, lvso(i), keepXso(i)*ones(size(1:lvso(i))), O*ones(size(1:lvso(i))), 500, 1e-10, 1, 0);
        beta = model.R*model.Q';

        srec = vcs*beta;
        [X,Y,T,AUC(i)] = perfcurve(val_y,srec,1);
        
    else
        keepXso(i) = nan;
        AUC(i) = [];
    end
    
end

AUCm = mean(AUC);

%% Show results

if opt == 1,
    fig_h = plot_vec(AUC,[],[],{'#Split','AUC'},[],1);
end


