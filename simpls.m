function [beta,W,P,Q,R,model] = simpls(X,Y,varargin)

% Simpls algorithm for Partial Least Squares. De Jong, Sijmen. "SIMPLS: an 
% alternative approach to partial least squares regression." Chemometrics 
% and intelligent laboratory systems 18.3 (1993): 251-263.
%
% beta = simpls(xcs,ycs)     % minimum call
% [beta,W,P,Q,R,model] = simpls(xcs,ycs,'LatVars',lvs)     % complete call
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% ycs: [NxO] preprocessed billinear data set of responses
%
% Optional INPUTS (parameter):
%
% 'LatVars': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:size(XX)
%
%
% OUTPUTS:
%
% beta: [MxO] matrix of regression coefficients: W*inv(P'*W)*Q'
%
% W: [MxA] matrix of weights
%
% P: [MxA] matrix of x-loadings
%
% Q: [OxA] matrix of y-loadings
%
% R: [MxA] matrix of modified weights: W*inv(P'*W)
%
% model: structure that contains model information
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Xcs = preprocess2D(X);
% Ycs = preprocess2D(Y);
% lvs = 1:10;
% [beta,W,P,Q,R] = simpls(Xcs,Ycs,'LatVars',lvs);
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
[N,M] = size(X);
O = size(Y, 2);

% Introduce optional inputs as parameters
p = inputParser;
addParameter(p,'LatVars',0:rank(X));   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LatVars;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
lvs(find(lvs>M)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (isequal(size(X), [N M]), 'Dimension Error: parameter ''X'' must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Y), [N O]), 'Dimension Error: parameter ''Y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LatVars'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LatVars'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

ncomp = max(lvs);
[n,dx] = size(X);
dy = size(Y,2);

% Preallocate outputs
Xloadings = zeros(dx,ncomp);
Yloadings = zeros(dy,ncomp);
Weights = zeros(dx,ncomp);


% An orthonormal basis for the span of the X loadings, to make the successive
% deflation X'*Y simple - each new basis vector can be removed from Cov
% separately.
V = zeros(dx,ncomp);

Cov = X'*Y;
for i = 1:ncomp
    % Find unit length ti=X*ri and ui=Y*ci whose covariance, ri'*X'*Y*ci, is
    % jointly maximized, subject to ti'*tj=0 for j=1:(i-1).
    [ri,si,ci] = svd(Cov,'econ'); 
    ri = ri(:,1); ci = ci(:,1); si = si(1);
    ti = X*ri;
    normti = norm(ti); ti = ti ./ normti; % ti'*ti == 1
    Xloadings(:,i) = X'*ti;
    
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
T = X*R;


model.var = sum(sum(X.^2));
model.lvs = 1:size(P,2);
model.loads = P;
model.yloads = Q;
model.weights = W;
model.scores = T;
model.type = 'PLS';
        