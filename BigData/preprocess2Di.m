function [xp,average,scale,N] = preprocess2Di(x,prep,ndim,lambda,average,scale,N,weight)

% Iteratively preprocess 2-way data.
%
% [xp,average,scale] = preprocess2Di(x)          % for mean centering
% [xp,average,scale] = preprocess2Di(x,prep)     % equivalent to preprocess2D
% [xp,average,scale,N] = preprocess2Di(x,prep,ndim,lambda,average,scale,N,weight) % complete call
%
% INPUTS:
%
% x: [NxM] Two-way data matrix
%
% prep: [1x1] preprocesing of the data
%       0: no preprocessing.
%       1: mean centering
%       2: auto-scaling (default, it centers and scales data so that each variable 
%           has variance 1)  
%       3: scaling (it scales previously centered data so that each variable 
%           has variance 1)  
%
% ndim: [1x1] 0 observations (by default), otherwise variables.
%
% lambda [1x1] forgetting factor between 0 (fast adaptation) and 1 (long
%   history)
%
% average: [1x(M or N)] previous sample average according to the preprocessing
%   method.
%
% scale: [1x(M or N)] previous sample scale according to the preprocessing
%   method.
%
% N: [1x1] number of effective observations in the model.
%
% weight: [1x(M or N)] weight applied after preprocessing scaling. Set to 1 
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
% X = simuleMV(10,10,8);
% [Xcs,av,sc] = preprocess2Di(X,2,0);
% plot_vec([av' sc'],[],[],{'Average & Std Dev' ,''},[], 0);
%
% X = simuleMV(10,10,8);
% [Xcs,av,sc] = preprocess2Di(X,2,0,0.9,av,sc);
% plot_vec([av' sc'],[],[],{'Average & Std Dev' ,''},[], 0);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 27/May/2017
%
% Copyright (C) 2017  University of Granada, Granada
% Copyright (C) 2017  Jose Camacho Paez
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
n = size(x, 1);
M = size(x, 2);

if nargin < 2 || isempty(prep), prep = 2; end;
prep2=prep;
if prep2 == 3, prep2=2; end;
[xcs,av,sc] = preprocess2D(x,prep2);

if nargin < 3 || isempty(ndim), ndim = 0; end;
if nargin < 4 || isempty(lambda), lambda = 0; end;
if nargin < 5 || isempty(average), average = av; end;
if nargin < 6 || isempty(scale), scale = sc; end;
if nargin < 7 || isempty(N), N = 0; end;
if nargin < 8 || isempty(weight), 
    if ndim,
        weight = ones(1,n);
    else
        weight = ones(1,M);
    end;
end
    
% Convert column arrays to row arrays
if size(weight,2) == 1, weight = weight'; end;

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ndim), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lambda), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(N), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(weight), [1 M]), 'Dimension Error: 8th argument must be 1-by-M. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (prep>=0 && prep<=3 && isequal(fix(prep), prep), 'Value Error: 2nd argument must contain integers between 0 and 3. Type ''help %s'' for more info.', routine(1).name);
assert (ndim==0 || ndim == 1, 'Value Error: 3rd argument must be 0 or 1. Type ''help %s'' for more info.', routine(1).name);
assert (lambda>=0 && lambda<=1, 'Value Error: 4th argument must be between 0 and 1. Type ''help %s'' for more info.', routine(1).name);
assert (N>=0, 'Value Error: 7th argument must be positive. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(weight<0)) && isempty(find(weight==Inf)), 'Value Error: 8th argument must contain positive values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

if ndim,
    x = x';
    s = size(x);
end

acc = average*N;
acc2 = (scale.^2)*max(N-1,0);
N = lambda*N + n; % update number of elements
switch prep,
    
    case 1, % mean centering
        acc = lambda*acc + sum(x,1); % update accumulate 
        average = acc/N;    
        scale = ones(1,M);
        xp = x - ones(n,1)*average;
        
    case 2, % auto-scaling
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
        
    case 3, % scaling  
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
        
    otherwise, % No preprocessing 
        average = zeros(1,M);     
        scale = ones(1,M); 
        xp = x;
end

xp = xp.*(ones(n,1)*weight);

if ndim,
    xp = xp';
end

