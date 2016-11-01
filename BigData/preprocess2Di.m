function [xp,average,scale,N] = preprocess2Di(x,prep,ndim,lambda,average,scale,N,weight)

% Iteratively preprocess 2-way data.
%
% [xp,average,scale] = preprocess2Di(x)          % for mean centering
% [xp,average,scale] = preprocess2Di(x,prep)     % equivalent to preprocess2D
% [xp,average,scale,N] = preprocess2Di(x,prep,ndim,lambda,average,scale,N,weight) % complete call
%
% INPUTS:
%
% x: (NxM) Two-way data matrix, N(observations) x M(variables)
%
% prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering (default) 
%       2: auto-scaling (centers and scales data so that each variable 
%           has variance 1)  
%       3: scaling (scales previously centered data so that each variable 
%           has variance 1)  
%
% ndim: (1x1) 0 observations (by default), otherwise variables.
%
% lambda (1x1) forgetting factor between 0 (fast adaptation) and 1 (long
%   history)
%
% average: (1x(M or N)) previous sample average according to the preprocessing
%   method.
%
% scale: (1x(M or N)) previous sample scale according to the preprocessing
%   method.
%
% N: (1x1) number of effective observations in the model.
%
% weight: (1x(M or N)) weight applied after preprocessing scaling. Set to 1 
% by defect.
%
% OUTPUTS:
%
% xp: (NxM) preprocessed data.
%
% average: (1x(M or N)) sample average according to the preprocessing method.
%
% scale: (1x(M or N)) sample scale according to the preprocessing method.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 03/Jun/15.
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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
    
%% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if nargin < 2 || isempty(prep), prep = 1; end;
if (prep<0||prep>3), error('Incorrect value of prep.'); end;
if nargin < 3 || isempty(ndim), ndim = 0; end;
if nargin < 4 || isempty(lambda), lambda = 0; end;
if nargin < 5 || isempty(average), average = 0; end;
if nargin < 6 || isempty(scale), scale = 0; end;
if nargin < 7 || isempty(N), N = 0; end;
if nargin < 8 || isempty(weight), 
    if ndim,
        weight = ones(1,s(1));
    else
        weight = ones(1,s(2));
    end;
end
    
% Computation

if ndim,
    x = x';
    s = size(x);
end

acc = average*N;
acc2 = (scale.^2)*max(N-1,0);
N = lambda*N + s(1); % update number of elements
switch prep,
    
    case 1, % mean centering
        acc = lambda*acc + sum(x,1); % update accumulate 
        average = acc/N;    
        scale = ones(1,s(2));
        xp = x - ones(s(1),1)*average;
        
    case 2, % auto-scaling
        acc = lambda*acc + sum(x,1); % update accumulate 
        average = acc/N;   
        xc = x - ones(s(1),1)*average; 
        acc2 = lambda*acc2 + sum(xc.^2,1);% update variability  
        scale = sqrt(acc2/(N-1));
        mS = min(scale(find(scale)));
        if isempty(mS)
            mS = 2;
        end
        scale(find(scale==0))=mS/2; % use 1 by default may reduce detection of anomalous events 
        xp = xc./(ones(s(1),1)*scale);
        
    case 3, % scaling  
        average = zeros(1,s(2)); 
        acc2 = lambda*acc2 + sum(x.^2,1);% update variability  
        if acc2 < 0, pause, end
        scale = sqrt(acc2/(N-1));
        mS = min(scale(find(scale)));
        if isempty(mS)
            mS = 2;
        end
        scale(find(scale==0))=mS/2; % use 1 by default may reduce detection of anomalous events  
        xp = x./(ones(s(1),1)*scale);
        
    otherwise, % No preprocessing 
        average = zeros(1,s(2));     
        scale = ones(1,s(2)); 
        xp = x;
end

xp = xp.*(ones(s(1),1)*weight);

if ndim,
    xp = xp';
end

