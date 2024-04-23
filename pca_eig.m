function model = pca_eig(xcs,varargin)

% Principal Component Analysis based on the eigendecompostion of XX.
%
% p = pca_pp(xcs)     % minimum call
% [p,t,model] = pca_pp(xcs,'Pcs',pcs)     % complete call
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% Optional INPUTS (parameter):
%
% 'Pcs': [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:min(size(xcs))
%
%
% OUTPUTS:
%
% p: [MxA] matrix of loadings.
%
% t: [NxA] matrix of scores.
%
% model: structure that contains model information.
%
%
% EXAMPLE OF USE: Random data:
%
% X = simuleMV(20,10,'LevelCorr',8);
% Xcs = preprocess2D(X,'Preprocessing',2);
% pcs = 1:3;
% model = pca_eig(Xcs,'Pcs',pcs)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Apr/2024
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(xcs, 1);
M = size(xcs, 2);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Pcs',0:rank(xcs));   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
pcs = p.Results.Pcs;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
pcs(find(pcs>size(xcs,2))) = [];
pcs(find(pcs>rank(xcs))) = [];
A = length(pcs);

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''Pcs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''Pcs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

if N>M,
    XX = xcs'*xcs;
    [p,D] = eig(XX);
    [kk,ind] = sort(real(diag(D)),'descend');
    p = p(:,ind);
    t = xcs*p;
else,
    XX = xcs*xcs';
    [t,D] = eig(XX);
    s = real(sqrt(real(diag(D))));
    [kk,ind] = sort(s,'descend');
    t = t(:,ind).*(ones(N,1)*s(ind)');
    p = xcs'*t;
    for i=1:size(p,2)
        p(:,i) = p(:,i)/sqrt(p(:,i)'*p(:,i));
    end
end

p = p(:,pcs);
t = t(:,pcs);

model.var = trace(XX);
model.lvs = 1:size(p,2);
model.loads = p;
model.scores = t;
model.type = 'PCA';
        



