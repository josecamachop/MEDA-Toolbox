function [xp,average,scale,N] = preprocess2Di(x,varargin)

% Iteratively preprocess 2-way data.
%
% [xp,average,scale] = preprocess2Di(x)          % for mean centering
% [xp,average,scale] = preprocess2Di(x,'Preprocessing',prep)     % equivalent to preprocess2D
% [xp,average,scale,N] = preprocess2Di(x,'Preprocessing',prep,'Ndim',ndim,'Lambda',lambda,'Average',average,'Scale',scale,'N',N,'Weight',weight) % complete call
%
% INPUTS:
%
% x: [NxM] Two-way data matrix
%
% Optional INPUTS (parameters):
%
% 'Preprocessing': [1x1] preprocesing of the data
%       0: no preprocessing.
%       1: mean centering
%       2: auto-scaling (default, it centers and scales data so that each variable 
%          has variance 1)  
%       3: scaling (it scales previously centered data so that each variable 
%          has variance 1)  
%
% 'Ndim': [1x1] 0 observations (by default), otherwise variables.
%
% 'Lambda': [1x1] forgetting factor between 0 (fast adaptation) and 1 (long
%   history)
%
% 'Average': [1x(M or N)] previous sample average according to the preprocessing
%   method.
%
% 'Scale': [1x(M or N)] previous sample scale according to the preprocessing
%   method.
%
% 'N': [1x1] number of effective observations in the model.
%
% 'Weight': [1x(M or N)] weight applied after preprocessing scaling. Set to 1 
% by defect.
%
%
% OUTPUTS:
%
% xp: [NxM] preprocessed data.
%
% average: [1x(M or N)] sample average according to the preprocessing method.
%
% scale: [1x(M or N)] sample scale according to the preprocessing method.
%
%
% EXAMPLE OF USE: Random data:
%
% X = simuleMV(10,10,'LevelCorr',8);
% [Xcs,av,sc] = preprocess2Di(X,'Preprocessing',2,'Ndim',0);
% plotVec([av' sc'],'VecLabel',{'Average','Std Dev'},'PlotType','Lines')
%
% X = simuleMV(10,10,'LevelCorr',8);
% [Xcs,av,sc] = preprocess2Di(X,'Preprocessing',2,'Ndim',0,'Lambda',0.9,'Average',av,'Scale',sc);
% plotVec([av' sc'],'VecLabel',{'Average','Std Dev'},'PlotType','Lines')
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 28/May/2026
%
% Copyright (C) 2026  University of Granada, Granada
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
routine = dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
n = size(x, 1);
M = size(x, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Preprocessing',2);   
addParameter(p,'Ndim',0);
addParameter(p,'Lambda',0);
addParameter(p,'Average',[]);
addParameter(p,'Scale',[]);
addParameter(p,'N',0);
addParameter(p,'Weight',[]);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
prep = p.Results.Preprocessing;
ndim = p.Results.Ndim;
lambda = p.Results.Lambda;
average = p.Results.Average;
scale = p.Results.Scale;
N = p.Results.N;
weight = p.Results.Weight;

% Compute dependencies for missing parameters
prep2 = prep;
if prep2 == 3, prep2 = 2; end;
[xcs,av,sc] = preprocess2D(x,'Preprocessing',prep2);

if isempty(average), average = av; end
if isempty(scale), scale = sc; end

if isempty(weight) 
    if ndim
        weight = ones(1,n);
    else
        weight = ones(1,M);
    end
end
    
% Convert column arrays to row arrays
if size(weight,2) == 1, weight = weight'; end;

% Determine expected dimension for checks based on Ndim
if ndim
    expected_dim = n;
else
    expected_dim = M;
end

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error in parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ndim), [1 1]), 'Dimension Error in parameter ''Ndim'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lambda), [1 1]), 'Dimension Error in parameter ''Lambda'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(N), [1 1]), 'Dimension Error in parameter ''N'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(weight), [1 expected_dim]), 'Dimension Error in parameter ''Weight'' must be 1-by-%d. Type ''help %s'' for more info.', expected_dim, routine(1).name);

% Validate values of input data
assert (prep>=0 && prep<=3 && isequal(fix(prep), prep), 'Value Error: parameter ''Preprocessing'' must contain integers between 0 and 3. Type ''help %s'' for more info.', routine(1).name);
assert (ndim==0 || ndim == 1, 'Value Error: parameter ''Ndim'' must be 0 or 1. Type ''help %s'' for more info.', routine(1).name);
assert (lambda>=0 && lambda<=1, 'Value Error: parameter ''Lambda'' must be between 0 and 1. Type ''help %s'' for more info.', routine(1).name);
assert (N>=0, 'Value Error: parameter ''N'' must be positive. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(weight<0)) && isempty(find(weight==Inf)), 'Value Error: parameter ''Weight'' must contain positive values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

if ndim
    x = x';
    s = size(x);
end

acc = average*N;
acc2 = (scale.^2)*max(N-1,0);
N = lambda*N + n; % update number of elements
switch prep
    
    case 1 % mean centering
        acc = lambda*acc + sum(x,1); % update accumulate 
        average = acc/N;    
        scale = ones(1,M);
        xp = x - ones(n,1)*average;
        
    case 2 % auto-scaling
        acc = lambda*acc + sum(x,1); % update accumulate 
        average = acc/N;   
        xc = x - ones(n,1)*average; 
        acc2 = lambda*acc2 + sum(xc.^2,1);% update variability  
        scale = sqrt(acc2/(N-1));
        mS = min(scale(find(scale)));
        if isempty(mS)
            mS = 2;
        end
        scale(find(scale==0))=mS/2; % use 1 by default may reduce detection of anomalous events 
        xp = xc./(ones(n,1)*scale);
        
    case 3 % scaling  
        average = zeros(1,M); 
        acc2 = lambda*acc2 + sum(x.^2,1);% update variability  
        if acc2 < 0, pause, end
        scale = sqrt(acc2/(N-1));
        mS = min(scale(find(scale)));
        if isempty(mS)
            mS = 2;
        end
        scale(find(scale==0))=mS/2; % use 1 by default may reduce detection of anomalous events  
        xp = x./(ones(n,1)*scale);
        
    otherwise % No preprocessing 
        average = zeros(1,M);     
        scale = ones(1,M); 
        xp = x;
end

xp = xp.*(ones(n,1)*weight);

if ndim
    xp = xp';
end

