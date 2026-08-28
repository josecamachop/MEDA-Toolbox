function model = kernelpls(xcs,ycs,varargin)

% Kernel algorithm for Partial Least Squares. References:
% F. Lindgren, P. Geladi and S. Wold, J. Chemometrics, 7, 45 (1993).
% S. De Jong and C.J.F. Ter Braak, J. Chemometrics, 8, 169 (1994).
% B.S. Dayal and J.F. MacGregor. J. Chemometrics, 11, 73–85 (1997). Main
% code is extracted from the last reference (Algorithm v1)
%
% model = kernelpls(X,Y)     % minimum call
%
% See also: kernelpls2, pcaEig, asca, gpls, sparsepls2
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
%   first two LVs). By default, lvs = 0:size(xcs,2). 
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
%   residuals: [NxM] matrix of x-residuals
%   type: 'PLS'
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Xcs = preprocess2D(X,'Preprocessing',2);
% Ycs = preprocess2D(Y,'Preprocessing',2);
% lvs = 1:10;
% model = kernelpls(Xcs,Ycs,'LVs',lvs);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 21/Aug/2026
% Dependencies: Matlab R2024b, MEDA v1.15
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
[N,M] = size(xcs);
O = size(ycs, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'LVs',0:size(xcs,2));  
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end

% Preprocessing
lvs = unique(lvs);
lvs(lvs<0) = [];
lvs(lvs>rank(xcs)) = [];
A = max(lvs);

% Validate dimensions of input data
assert (isequal(size(xcs), [N M]), 'Dimension Error: parameter ''X'' must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ycs), [N O]), 'Dimension Error: parameter ''Y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 size(lvs,2)]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

T = zeros(N,A);
W = zeros(M,A);
P = zeros(M,A);
Q = zeros(O,A);
R = zeros(M,A);
XY = xcs'*ycs;
for i=1:A
    if O==1 
        W(:,i) = XY; 
    else 
        [C,D] = eig(XY'*XY); 
        q = C(:,diag(D)==max(diag(D))); 
        W(:,i)=(XY*q); 
    end
    W(:,i) = W(:,i)/sqrt(W(:,i)'*W(:,i)); 
    R(:,i) = W(:,i); 
    for j=1:i-1
        R(:,i) = R(:,i)-(P(:,j)'*W(:,i))*R(:,j);
    end
    T(:,i) = xcs*R(:,i); 
    tt = T(:,i)'*T(:,i);
    P(:,i) = xcs'*T(:,i)/tt; 
    Q(:,i) = (R(:,i)'*XY)'/tt; 
    XY = XY-(P(:,i)*Q(:,i)')*tt;
end

lvs(lvs==0) = [];
P = P(:,lvs);
Q = Q(:,lvs);
W = W(:,lvs);
T = T(:,lvs);
R = W*pinv(P'*W);

beta = R*Q';

model.var = trace(xcs'*xcs);
model.lvs = 1:size(P,2);
model.loads = P;
model.yloads = Q;
model.weights = W;
model.altweights = R;
model.scores = T;
model.beta = beta;
model.residuals = xcs - T*P';
model.type = 'PLS';
