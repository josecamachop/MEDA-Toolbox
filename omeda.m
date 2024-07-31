
function omeda_vec = omeda(testcs,dummy,R,varargin)

% Observation-based Missing data methods for Exploratory Data Analysis 
% (oMEDA). The original paper is Journal of Chemometrics, 2011, 25 
% (11): 592-600. This algorithm follows the direct computation for
% Known Data Regression (KDR) missing data imputation.
%
% omeda_vec = omeda(testcs,dummy,R) % minimum call
% omeda_vec = omeda(testcs,dummy,R,'OutSubspace',Q) % complete call
%
%
% INPUTS:
%
% testcs: [NxM] preprocessed billinear data set with the observations to be 
%   compared.
%
% dummy: [Nx1] dummy variable containing weights for the observations to 
%   compare, and 0 for the rest of observations.
%
% R: [MxA] Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (testcs = T*P'), this is the matrix of loadings 
%   P. For PLS (Y = testcs*W*inv(P'*W)*Q), this matrix is W*inv(P'*W). For the 
%   original space (default) the identity matrix is used.  
%
% Optional INPUTS (parameter):
%
% 'OutSubspace': [MxA] Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (testcs = T*P'), this is the matrix of 
%   loadings P. For PLS (Y = testcs*W*inv(P'*W)*Q), this matrix is also P. 
%   For the original space the identity matrix is used. Q=R is used by 
%   default. 
%
%
% OUTPUTS:
%
% omeda_vec: [Mx1] oMEDA vector.
%
%
%
% EXAMPLE OF USE: oMEDA on PCA, anomaly on first observation and first 2
% variables.
%
% n_obs = 100;
% n_vars = 10;
% n_PCs = 10;
% X = simuleMV(n_obs,n_vars,'LevelCorr',6);
% [Xcs, m, sc] = preprocess2D(X,'Preprocessing',2);
% pcs = 1:n_PCs;
% p = pca_pp(Xcs,'Pcs',pcs);
% 
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,'LevelCorr',6,'Covar',cov(X)*(n_obst-1));
% test(1,1:2) = 10*max(abs(X(:,1:2))); 
% dummy = zeros(10,1);
% dummy(1) = 1;
% testcs = preprocess2Dapp(test,m,'SDivideTest',sc);
% 
% omeda_vec = omeda(testcs,dummy,p);
% 
% plot_vec(omeda_vec);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 22/Apr/2024.
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
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(testcs, 1);
M = size(testcs, 2);
A = size(R, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'OutSubspace',[]);     
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
Q = p.Results.OutSubspace;
if isempty(Q), Q = R; end;

% Convert row arrays to column arrays
if size(dummy,1) == 1, dummy = dummy'; end;

% Validate dimensions of input data
assert (isequal(size(dummy), [N 1]), 'Dimension Error: parameter ''dummy'' must be 1-by-N. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(R), [M A]), 'Dimension Error: parameter ''R'' must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Q), [M A]), 'Dimension Error: parameter ''OutSubspace'' must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);



%% Main code

ind=find(dummy>0);
dummy(ind) = dummy(ind)/max((dummy(ind)));
ind=find(dummy<0);
dummy(ind) = -dummy(ind)/min((dummy(ind)));

xA = testcs*R*Q';
sumA = xA'*dummy;
sum = testcs'*dummy;

omeda_vec = ((2*sum-sumA).*abs(sumA))./sqrt(dummy'*dummy);

    


        