function model = pcaEig(xcs,varargin)

% Principal Component Analysis based on the eigendecompostion of XX.
%
% model = pcaEig(xcs)     % minimum call
%
% See also: kernelpls, simpls, asca, missTsr2D
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
%
% Optional INPUTS (parameter):
%
% 'PCs': [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:min(size(xcs))
%
%
% OUTPUTS:
%
% model: structure that contains model information
%
%
% EXAMPLE OF USE: Random data:
%
% X = simuleMV(20,10,'LevelCorr',8);
% Xcs = preprocess2D(X,'Preprocessing',2);
% pcs = 1:3;
% model = pcaEig(Xcs,'PCs',pcs)
%
%
% EXAMPLE OF USE: Compare PCA algorithms for wide data
%
% X = randn(20,1e4);
% Xcs = preprocess2D(X,'Preprocessing',2);
% tic, model = pcaEig(Xcs,'PCs',1:2); toc
% tic, coefSvd = pca(Xcs,'Algorithm','svd','Centered',false,'NumComponents',2); toc
% tic, coefEig = pca(Xcs,'Algorithm','eig','Centered',false,'NumComponents',2); toc
%
% for pc = 1:2,
%   s = sign(model.loads(abs(model.loads(:,pc))==max(abs(model.loads(:,pc))),pc));
%   model.loads(:,pc) = s*model.loads(:,pc);
%   s = sign(coefSvd(abs(coefSvd(:,pc))==max(abs(coefSvd(:,pc))),pc));
%   coefSvd(:,pc) = s*coefSvd(:,pc);
%   s = sign(coefEig(abs(coefEig(:,pc))==max(abs(coefEig(:,pc))),pc));
%   coefEig(:,pc) = s*coefEig(:,pc);
% end
%
% norm(model.loads-coefSvd)
% norm(coefEig-coefSvd)
% norm(model.loads-coefEig)
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 12/Jul/2026
% Dependencies: Matlab R2024b, MEDA v1.14
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
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(xcs, 1);
M = size(xcs, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'PCs',1:min(N,M));   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
pcs = p.Results.PCs;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs>min(N,M))) = [];
A = length(pcs);

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''Pcs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(pcs<=0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''Pcs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

if N>M
    XX = xcs'*xcs;
    [p,D] = eig(1/(N-1)*XX);
    [lambda,ind] = sort(diag(D),'descend');
    p = p(:,ind);
    t = xcs*p;
else
    XX = xcs*xcs';
    [t,D] = eig(1/(N-1)*XX);
    [lambda,ind] = sort(diag(D),'descend');
    t = t(:,ind).*(sqrt((N-1))*ones(N,1)*real(sqrt(lambda')));
    p = xcs'*t;
    for i=1:size(p,2)
        p(:,i) = p(:,i)/sqrt(p(:,i)'*p(:,i));
    end
end

p = p(:,pcs);
t = t(:,pcs);

model.lambda = lambda;
model.var = trace(XX);
model.lvs = 1:size(p,2);
model.loads = p;
model.scores = t;
model.residuals = xcs - t*p';
model.type = 'PCA';
        



