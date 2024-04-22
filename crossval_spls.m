function [cumpress,press,nze] = crossval_spls(x,y,varargin)

% Row-wise k-fold (rkf) cross-validation for square-prediction-errors computing in SPLS.
%
% cumpress = crossval_spls(x,y) % minimum call
% [cumpress,press,nze] =
% crossval_spls(x,y,'LatVars',lvs,'KeepXBlock',keepXs,'MaxBlock',blocks_r,'PreprocessingX',prepx,'PreprocessingY',prepy,'Option',opt) % complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
% Optional INPUTS:
%
% 'LatVars': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x)
%
% 'KeepXBlock': [1xK] Numbers of x-block variables kept per latent variable modeled. By default, keepXs = 1:M
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
% 'Option': (str or num) options for data plotting.
%       0: no plots.
%       1: plot (default)
%
%
% OUTPUTS:
%
% cumpress: [AxK] Cumulative PRESS
%
% press: [AxKxO] PRESS per variable
%
% nze: [AxK] Non-zero elements in the regression coefficient matrix.
%
%
% EXAMPLE OF USE: Random data with structural relationship
% 
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% keepXs = 1:10;
% [cumpress,press,nze] = crossval_spls(X,Y,'KeepXBlock',keepXs);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 22/Apr/24.
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
M = size(x, 2);
O = size(y, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
lat=0:rank(x);
addParameter(p,'LatVars',lat'); 
keep = 1:M;
addParameter(p,'KeepXBlock',keep);
addParameter(p,'MaxBlock',N);
addParameter(p,'PreprocessingX',2);   
addParameter(p,'PreprocessingY',2);
addParameter(p,'Option',1);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LatVars;
keepXs = p.Results.KeepXBlock;
blocks_r = p.Results.MaxBlock;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
opt = p.Results.Option;

% Extract LatVars and KeepXBlock length
A = length(lvs);
J =  length(keepXs);

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;
if size(keepXs,2) == 1, keepXs = keepXs'; end;

% Validate dimensions of input data
assert (isequal(size(y), [N O]), 'Dimension Error: parameter ''y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LatVars'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(keepXs), [1 J]), 'Dimension Error: parameter ''KeepXBlock'' must be 1-by-J. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(blocks_r), [1 1]), 'Dimension Error: parameter ''MaxBlock'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: parameter ''Option'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
lvs = unique(lvs);
keepXs = unique(keepXs);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: parameter ''LatVars'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(lvs), lvs), 'Value Error: parameter ''LatVars'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(keepXs), keepXs), 'Value Error: parameter ''KeepXBlock'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(blocks_r), blocks_r), 'Value Error: parameter ''MaxBlock'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r>2, 'Value Error: parameter ''MaxBlock'' must be above 2. Type ''help %s'' for more info.', routine(1).name);
assert (blocks_r<=N, 'Value Error: parameter ''MaxBlock'' must be at most N. Type ''help %s'' for more info.', routine(1).name);


%% Main code

% Initialization
press = zeros(length(lvs),length(keepXs),O);
nze = zeros(length(lvs),length(keepXs));


rows = rand(1,N);
[a,r_ind]=sort(rows);
elem_r=N/blocks_r;

% Cross-validation
        
for i=1:blocks_r,
    
    ind_i = r_ind(round((i-1)*elem_r+1):round(i*elem_r)); % Sample selection
    i2 = ones(N,1);
    i2(ind_i)=0;
    sample = x(ind_i,:);
    calibr = x(find(i2),:); 
    sample_y = y(ind_i,:);
    calibr_y = y(find(i2),:); 

    [ccs,av,st] = preprocess2D(calibr,'Preprocessing',prepx);
    [ccs_y,av_y,st_y] = preprocess2D(calibr_y,'Preprocessing',prepy);
        
    scs = preprocess2Dapp(sample,av,'SDivideTest',st);
    scs_y = preprocess2Dapp(sample_y,av_y,'SDivideTest',st_y);
    
    if  ~isempty(find(lvs)),
        
        for lv=1:length(lvs),

            for keepX=1:length(keepXs),
                
                if lvs(lv),
                    model = sparsepls2(ccs, ccs_y, lvs(lv), keepXs(keepX)*ones(size(1:lvs(lv))), O*ones(size(1:lvs(lv))), 500, 1e-10, 1, 0);
                    beta = model.R*model.Q';

                    srec = scs*beta;
                    pem = scs_y-srec;

                    press(lv,keepX,:) = squeeze(press(lv,keepX,:))' + sum(pem.^2,1);
					nze(lv,keepX) = nze(lv,keepX) + length(find(beta)); 
                else
                    press(lv,keepX,:) = squeeze(press(lv,keepX,:))' + sum(scs_y.^2,1);
					nze(lv,keepX) = nze(lv,keepX) + M*O; 
                end
                
            end
            
        end
        
    else
        pem = scs_y;
        press = press + ones(length(keepXs),1)*sum(pem.^2,1);
		nze = nze + ones(length(keepXs),1)*M*O;
    end
    
end

cumpress = sum(press,3);

%% Show results

if opt == 1,
    fig_h = plot_vec(cumpress','EleLabel',keepXs,'XYLabel',{'#NZV','PRESS'},'Option','01','VecLabel',lvs); 
end

