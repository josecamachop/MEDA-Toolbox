function model = spcaZou(X,Gram,K,stop,varargin)

% SPCA Zou, Hastie, Tibshirani 2006 with enhancement to select the number
% of non-zero elements. Modified version of the code of the SpaSM Toolbox 
% (https://www2.imm.dtu.dk/projects/spasm/), Sjöstrand, Karl, et al. 
% "Spasm: A matlab toolbox for sparse statistical modeling." Journal of 
% Statistical Software 84.10 (2018). The input parameters are consistent 
% with spca_zouhastie in that toolbox, but the output is modified to
% follow the interpretation approach of "All sparse PCA models are wrong, 
% but some are useful. Part III: model interpretation", Chemometrics and
% Intelligent Laboratory Systems. Furthermore, this implementation only
% considers soft-thresholding, especially useful for wide matrices (M>>N). 
% Please, use SpaSM or the original code by Zou et al. in R for other 
% situations.  
%
% [P,Q] = spcaZou(X,Gram,K,stop)     % minimum call
%
%
% INPUTS:
%
% X: [NxM] preprocessed billinear data set 
%
% Gram: [NxN] Gram matrix from X (either X or Gram have to be inputted) 
%
% K: [1x1] Number of sparse components.
%
% stop: [1x1] is the stopping criterion. If STOP is negative, its absolute
%    (integer) value corresponds to the desired number of non-zero
%    variables. If STOP is positive, it corresponds to an upper bound on
%    the L1-norm of the BETA coefficients. STOP = 0 results in a regular
%    PCA.
%
% Optional INPUTS (parameters):
%
% 'Tolerance': [1x1] tolerance value. By default, 1e-15.
%
% 'MaxIters': [1x1] maximum iterations. By default, 1e3.
%
%
% OUTPUTS:
%
% model: structure that contains model information
%   var: [1x1] xcs sum of squares
%   expvar: [1x1] explained sum of squares per component: tr(Q*A2'*A2*Q')
%   lvs: [1xA] latent variable numbers
%   loads: [MxA] matrix of x-loadings Q
%   weights: [MxA] matrix of weights P
%   altweights: [MxA] matrix of alternative weights R
%   scores: [NxA] matrix of x-scores T
%   type: 'sPCA'
%
%
% EXAMPLE OF USE: Pitprops
%
% var_l = {'topdiam' 'length'  'moist' 'testsg' 'ovensg' 'ringtop' 'ringbut' 'bowmax' 'bowdist' 'whorls'  'clear'  'knots' 'diaknot'};
%
% XX=[ 1.000  0.954  0.364  0.342 -0.129   0.313   0.496  0.424   0.592  0.545  0.084 -0.019   0.134
%      0.954  1.000  0.297  0.284 -0.118   0.291   0.503  0.419   0.648  0.569  0.076 -0.036   0.144
%      0.364  0.297  1.000  0.882 -0.148   0.153  -0.029 -0.054   0.125 -0.081  0.162  0.220   0.126
%      0.342  0.284  0.882  1.000  0.220   0.381   0.174 -0.059   0.137 -0.014  0.097  0.169   0.015
%     -0.129 -0.118 -0.148  0.220  1.000   0.364   0.296  0.004  -0.039  0.037 -0.091 -0.145  -0.208
%      0.313  0.291  0.153  0.381  0.364   1.000   0.813  0.090   0.211  0.274 -0.036  0.024  -0.329
%      0.496  0.503 -0.029  0.174  0.296   0.813   1.000  0.372   0.465  0.679 -0.113 -0.232  -0.424
%      0.424  0.419 -0.054 -0.059  0.004   0.090   0.372  1.000   0.482  0.557  0.061 -0.357  -0.202
%      0.592  0.648  0.125  0.137 -0.039   0.211   0.465  0.482   1.000  0.526  0.085 -0.127  -0.076
%      0.545  0.569 -0.081 -0.014  0.037   0.274   0.679  0.557   0.526  1.000 -0.319 -0.368  -0.291
%      0.084  0.076  0.162  0.097 -0.091  -0.036  -0.113  0.061   0.085 -0.319  1.000  0.029   0.007
%     -0.019 -0.036  0.220  0.169 -0.145   0.024  -0.232 -0.357  -0.127 -0.368  0.029  1.000   0.184
%      0.134  0.144  0.126  0.015 -0.208  -0.329  -0.424 -0.202  -0.076 -0.291  0.007  0.184   1.000];
%
% model = spcaZou([], XX, 6, -2);
% Q = model.loads;
% R = model.altweights;
% f = plotMap([Q*R'*XX*R*Q'],'VarsLabel',var_l);
% a = get(f,'Children');
% set(a(2),'XTickLabelRotation',45);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 4/Feb/2025
%
% Copyright (C) 2025  University of Granada, Granada
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
assert (nargin >= 4, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
    
% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Tolerance',1e-15);
addParameter(p,'MaxIters',1e3);          
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
tol = p.Results.Tolerance;
max_iter = p.Results.MaxIters;

% Replicate STOP value for all components if necessary
if length(stop) ~= K
  stop = stop(1)*ones(1,K);
end

if isempty(X) && isempty(Gram)
  error('SpaSM:spca', 'Must supply a data matrix or a Gram matrix or both.');
end

%% Main code

if isempty(X)
  % Infer X from X'*X
  [Vg Dg] = eig(Gram);
  X = Vg*sqrt(abs(Dg))*Vg';
end

[u,s,V] = svd(X, 'econ');

[n p] = size(X);

A = V(:,1:K);

B = [];
for k = 1:K
    if isempty(Gram)
        AXX = (A(:,k)'*X')*X;
    else
        AXX = A(:,k)'*Gram;
    end
    if stop(k) < 0 && -stop(k) < p
        sortedAXX = sort(abs(AXX), 'descend');
        B(:,k) = ( sign(AXX).*max(0, abs(AXX) - sortedAXX(-floor(stop(k)) + 1)) )';
    else
        B(:,k) = ( sign(AXX).*max(0, abs(AXX) - stop(k)) )';
    end

    % Normalize coefficients such that loadings has Euclidean length 1
    B_norm = sqrt(B(:,k)'*B(:,k));
    if B_norm == 0
        B_norm = 1;
    end
    B(:,k) = B(:,k)/B_norm;
    
end


stopit = false;
iter = 0;
while ~stopit

    if isempty(Gram)
        [u,s,v] = svd(X'*(X*B),0);
    else
        [u,s,v] = svd(Gram*B,0);
    end 
    A = u*v';
    
    Bk_old = B;
    for k = 1:K
        % Soft thresholding, calculate beta directly
        if isempty(Gram)
            AXX = (A(:,k)'*X')*X;
        else
            AXX = A(:,k)'*Gram;
        end
        if stop(k) < 0 && -stop(k) < p
            sortedAXX = sort(abs(AXX), 'descend');
            B(:,k) = ( sign(AXX).*max(0, abs(AXX) - sortedAXX(-floor(stop(k)) + 1)) )';
        else
            B(:,k) = ( sign(AXX).*max(0, abs(AXX) - stop(k)) )';
        end
        
        % Normalize coefficients such that loadings has Euclidean length 1
        B_norm = sqrt(B(:,k)'*B(:,k));
        if B_norm == 0
            B_norm = 1;
        end
        B(:,k) = B(:,k)/B_norm;
    end

    iter = iter+1;

    if norm(B-Bk_old) < tol || iter > max_iter
        stopit = true;
    end
end


P = B;
Q = A;
R =  Q*pinv(P'*Q);
A2 = X*R; 

for i=1:K 
    v(i) = trace(Q(:,i)*A2(:,i)'*A2(:,i)*Q(:,i)');
end

model.expvar = v;
model.var = sum(sum(X.^2));
model.lvs = 1:K;
model.loads = Q;
model.weights = P;
model.altweights = R;
model.scores = A2;
model.type = 'sPCA';
