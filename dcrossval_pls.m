function [Qm,Q,lvso] = dcrossval_pls(x,y,lvs,blocks_r,prepx,prepy,rep,opt)

% Row-wise k-fold (rkf) double cross-validation for PLS. The algorithm uses
% repetitions of the dCV loop to estimate the stability: see Szymanska, E., 
% Saccenti, E., Smilde, A.K., Weterhuis, J. Metabolomics (2012) 8: 3.
%
% Qm = dcrossval_pls(x,y) % minimum call
% [Qm,Q,lvso] = dcrossval_pls(x,y,lvs,blocks_r,prepx,prepy,rep,opt) % complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
% lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x)
%
% blocks_r: [1x1] maximum number of blocks of samples (N by default)
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
% rep: [1x1] number of repetitios for stability.
%
% opt: [1x1] options for data plotting
%       0: no plots
%       1: bar plot (default)
%
%
% OUTPUTS:
%
% Qm: [1x1] Mean Goodness of Prediction
%
% Q: [rep x 1] Goodness of Prediction
%
% lvso: [rep x blocks_r] optimum number of LVs in the inner loop
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% lvs = 0:10;
% Q = dcrossval_pls(X,Y,lvs);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/May/23
%
% Copyright (C) 2023  University of Granada, Granada
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
if nargin < 4 || isempty(blocks_r), blocks_r = N; end;
if nargin < 5 || isempty(prepx), prepx = 2; end;
if nargin < 6 || isempty(prepy), prepy = 2; end;
if nargin < 7 || isempty(rep), rep = 10; end;
if nargin < 8 || isempty(opt), opt = 1; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N O]), 'Dimension Error: 2nd argument must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocks_r), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(rep), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: 3rd argument must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: 3rd argumentmust contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocks_r), blocks_r), 'Value Error: 4th argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r>3, 'Value Error: 4th argument must be above 3. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r<=N, 'Value Error: 4th argument must be at most N. Type ''help %s'' for more info.', routine(1).name);


%% Main code

for j=1:rep
    % Cross-validation
    
    rows = rand(1,N);
    [a,r_ind]=sort(rows);
    elem_r=N/blocks_r;
    
    for i=1:blocks_r
        ind_i = r_ind(round((i-1)*elem_r+1):round(i*elem_r)); % Sample selection
        i2 = ones(N,1);
        i2(ind_i)=0;
        val = x(ind_i,:);
        rest = x(find(i2),:);
        val_y = y(ind_i,:);
        rest_y = y(find(i2),:);
        
        cumpress = crossval_pls(rest,rest_y,lvs,blocks_r-1,prepx,prepy,0);
        
        lvso(j,i) = lvs(find(cumpress==min(cumpress),1));
        
        [ccs,av,st] = preprocess2D(rest,prepx);
        [ccs_y,av_y,st_y] = preprocess2D(rest_y,prepy);
        
        vcs = preprocess2Dapp(val,av,st);
        vcs_y = preprocess2Dapp(val_y,av_y,st_y);
        
        beta = simpls(ccs,ccs_y,1:lvso(i));
        srec = vcs*beta;
        
        Qu(i) = sum(sum((vcs_y-srec).^2));
        Qd(i) = sum(sum(vcs_y.^2));
    end
    
    Q(j) = 1-sum(Qu)/sum(Qd);
    
end

Qm = mean(Q);

%% Show results

if opt == 1
   fig_h = plot_vec(Q,[],[],{'#Repetition','Goodness of Prediction'},[],1); 
end

