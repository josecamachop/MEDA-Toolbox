function [p,t] = pca_pp(xcs,pcs)

% Principal Component Analysis based on svd.
%
% p = pca_pp(xcs)     % minimum call
% [p,t] = pca_pp(xcs,pcs)     % complete call
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:min(size(xcs))
%
%
% OUTPUTS:
%
% p: [MxA] matrix of loadings.
%
% t: [NxA] matrix of scores.
%
%
% EXAMPLE OF USE: Random data:
%
% X = simuleMV(20,10,8);
% Xcs = preprocess2D(X,2);
% pcs = 1:3;
% [p,t] = pca_pp(Xcs,pcs);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 09/May/16.
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(xcs, 1);
M = size(xcs, 2);
if nargin < 2 || isempty(pcs), pcs = 0:size(xcs,2); end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
pcs(find(pcs>size(xcs,2))) = [];
A = length(pcs);

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 2nd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: 2nd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[u,d,p]=svd(xcs);
t = u*d;
p = p(:,pcs);
t = t(:,pcs);

        



