function [AUC,nze] = crossval_spls_da(x,y,lvs,keepXs,blocks_r,prepx,prepy,opt)

% Row-wise k-fold (rkf) cross-validation in SPLS-DA, restricted to one 
% response categorical variable of two levels.
%
% cumpress = crossval_spls_da(x,y) % minimum call
% [AUC,nze] =
% crossval_spls_da(x,y,lvs,keepXs,blocks_r,prepx,prepy,opt) % complete call
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
% keepXs: [1xK] Numbers of x-block variables kept per latent variable modeled. 
%   By default, keepXs = 1:M
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
% opt: (str or num) options for data plotting.
%       0: no plots.
%       1: plot (default)
%
%
% OUTPUTS:
%
% AUC: [AxK] Area Under the Curve in ROC
%
% nze: [AxK] Non-zero elements in the regression coefficient matrix.
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,8);
% Y = 2*(0.1*randn(20,1) + X(:,1)>0)-1;
% lvs = 0:10;
% keepXs = 1:10;
% [AUC,nze] = crossval_spls_da(X,Y,lvs,keepXs,5);
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
M = size(x, 2);
O = size(y, 2);

if nargin < 3 || isempty(lvs), lvs = 0:rank(x); end;
A = length(lvs);
if nargin < 4 || isempty(keepXs), keepXs = 1:M; end;
J =  length(keepXs);

vals = unique(y);
rep = sort(histc(y,vals),'descend');
N2 = rep(2);
if nargin < 5 || isempty(blocks_r), blocks_r = max(2,round(N2/2)); end;

if nargin < 6 || isempty(prepx), prepx = 2; end;
if nargin < 7 || isempty(prepy), prepy = 2; end;
if nargin < 8 || isempty(opt), opt = 1; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;
if size(keepXs,2) == 1, keepXs = keepXs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N 1]), 'Dimension Error: 2nd argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(keepXs), [1 J]), 'Dimension Error: 4th argument must be 1-by-J. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);
keepXs = unique(keepXs);

% Validate values of input data

assert (isempty(find(y~=1 & y~=-1)), 'Value Error: 2rd argument must not contain values different to 1 or -1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(lvs<0)), 'Value Error: 3rd argument must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: 3rd argument must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(keepXs), keepXs), 'Value Error: 4th argument must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocks_r), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocks_r), blocks_r), 'Value Error: 5th argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r>2, 'Value Error: 5th argument must be above 2. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r<=N2, 'Value Error: 5th argument must be at most %d. Type ''help %s'' for more info.', N2, routine(1).name);



%% Main code

% Initialization
AUC = zeros(length(lvs),length(keepXs));
nze = zeros(length(lvs),length(keepXs));

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
    sample = x(y1(ind_i),:);
    calibr = x(y1(find(i2)),:); 
    sample_y = y(y1(ind_i),:);
    calibr_y = y(y1(find(i2)),:); 
    
    ind_i = r_indn1(round((i-1)*elem_rn1+1):round(i*elem_rn1)); % Sample selection
    i2 = ones(length(yn1),1);
    i2(ind_i)=0;
    sample = [sample;x(yn1(ind_i),:)];
    calibr = [calibr;x(yn1(find(i2)),:)]; 
    sample_y = [sample_y;y(yn1(ind_i),:)];
    calibr_y = [calibr_y;y(yn1(find(i2)),:)];    
    
    [ccs,av,st] = preprocess2D(calibr,prepx);
    [ccs_y,av_y,st_y] = preprocess2D(calibr_y,prepy);
        
    scs = preprocess2Dapp(sample,av,st);
    
    if  ~isempty(find(lvs)),
        
        for lv=1:length(lvs),

            for keepX=1:length(keepXs),
                
                if lvs(lv),
                    model = sparsepls2(ccs, ccs_y, lvs(lv), keepXs(keepX)*ones(size(1:lvs(lv))), O*ones(size(1:lvs(lv))), 500, 1e-10, 1, 0);
                    beta = model.R*model.Q';

                    srec = scs*beta;
                    [X,Y,T,AUCt] = perfcurve(sample_y,srec,1);
                    
                    AUC(lv,keepX) = AUC(lv,keepX) + AUCt;
					nze(lv,keepX) = nze(lv,keepX) + length(find(beta)); 
                else
                    AUC(lv,keepX) = AUC(lv,keepX) + 0.5;
					nze(lv,keepX) = nze(lv,keepX) + M*O; 
                end
                
            end
            
        end
        
    else
        AUC = AUC + ones(length(keepXs),1)*0.5;
		nze = nze + ones(length(keepXs),1)*M*O;
    end
    
end
    
AUC = AUC/blocks_r;


%% Show results

if opt == 1,
    fig_h = plot_vec(AUC',keepXs,[],{'#NZV','AUC'},[],0,lvs); 
    legend('show')
end

