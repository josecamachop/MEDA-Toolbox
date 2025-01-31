function model = simpls(xcs,ycs,varargin)

% Simpls algorithm for Partial Least Squares. De Jong, Sijmen. "SIMPLS: an 
% alternative approach to partial least squares regression." Chemometrics 
% and intelligent laboratory systems 18.3 (1993): 251-263.
%
% model = simpls(xcs,ycs)     % minimum call
%
% See also: kernelpls, vpls, pcaEig, asca
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% ycs: [NxO] preprocessed billinear data set of responses
%
%
% Optional INPUTS (parameter):
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:size(XX)
%
%
% OUTPUTS:
%
% model: structure that contains model information
%   var: [1x1] xcs sum of squares
%   lvs: [1xA] latent variable numbers
%   loads: [MxA] matrix of x-loadings P
%   yloads: [OxA] matrix of y-loadings Q
%   weights: [MxA] matrix of weights W
%   altweights: [MxA] matrix of alternative weights R
%   scores: [NxA] matrix of x-scores T
%   beta: [MxO] matrix of regressors
%   type: 'PLS'
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Xcs = preprocess2D(X);
% Ycs = preprocess2D(Y);
% lvs = 1:10;
% model= simpls(Xcs,Ycs,'LVs',lvs);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 18/Nov/2024
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
[N,M] = size(xcs);
O = size(ycs, 2);

% Introduce optional inputs as parameters
p = inputParser;
addParameter(p,'LVs',0:rank(xcs));   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
lvs(find(lvs>M)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (isequal(size(xcs), [N M]), 'Dimension Error: parameter ''X'' must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ycs), [N O]), 'Dimension Error: parameter ''Y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

ncomp = max(lvs);
[n,dx] = size(xcs);
dy = size(ycs,2);

% Preallocate outputs
Xloadings = zeros(dx,ncomp);
Yloadings = zeros(dy,ncomp);
Weights = zeros(dx,ncomp);

% An orthonormal basis for the span of the X loadings, to make the successive
% deflation X'*Y simple - each new basis vector can be removed from Cov
% separately.
V = zeros(dx,ncomp);

Cov = xcs'*ycs;
for i = 1:ncomp
    % Find unit length ti=X*ri and ui=Y*ci whose covariance, ri'*X'*Y*ci, is
    % jointly maximized, subject to ti'*tj=0 for j=1:(i-1).
    [ri,si,ci] = svd(Cov,'econ'); 
    ri = ri(:,1); ci = ci(:,1); si = si(1);
    ti = xcs*ri;
    normti = norm(ti); ti = ti ./ normti; % ti'*ti == 1
    Xloadings(:,i) = xcs'*ti;
    
    qi = si*ci/normti; % = Y'*ti
    Yloadings(:,i) = qi;
    
    Weights(:,i) = ri ./ normti; % rescaled to make ri'*X'*X*ri == ti'*ti == 1

    % Update the orthonormal basis with modified Gram Schmidt (more stable),
    % repeated twice (ditto).
    vi = Xloadings(:,i);
    for repeat = 1:2
        for j = 1:i-1
            vj = V(:,j);
            vi = vi - (vj'*vi)*vj;
        end
    end
    vi = vi ./ norm(vi);
    V(:,i) = vi;

    % Deflate Cov, i.e. project onto the ortho-complement of the X loadings.
    % First remove projections along the current basis vector, then remove any
    % component along previous basis vectors that's crept in as noise from
    % previous deflations.
    Cov = Cov - vi*(vi'*Cov);
    Vi = V(:,1:i);
    Cov = Cov - Vi*(Vi'*Cov);
end

% Postprocessing
R = Weights*inv(Xloadings'*Weights);
R = R(:,lvs);
W = Weights(:,lvs);
P = Xloadings(:,lvs);
Q = Yloadings(:,lvs);
beta=R*Q';

for i=1:size(R,2)
    P(:,i) = P(:,i)*norm(R(:,i));
    Q(:,i) = Q(:,i)*norm(R(:,i));
    W(:,i) = W(:,i)/norm(R(:,i));
    R(:,i) = R(:,i)/norm(R(:,i));
end
T = xcs*R;

model.var = sum(sum(xcs.^2));
model.lvs = 1:size(P,2);
model.loads = P;
model.yloads = Q;
model.weights = W;
model.altweights = R;
model.scores = T;
model.beta = beta;
model.type = 'PLS';
        