
function [Dst,Qst] = mspc(testcs,invCT,R,Q)

% Multivariate Statistical Process Control statistics
%
% Dst = mspc(testcs) % minimum call (only Q-st)
% Dst = mspc(testcs,invCT,R) % minimum call for D-st and Q-st
% [Dst,Qst] = mspc(testcs,invCT,R,Q) % complete call
%
%
% INPUTS:
%
% testcs: [NxM] preprocessed billinear data set with the observations to be 
%   monitored.
%
% invCT: [AxA] inverse of covariance matrix of T, where T are the 
%   calibration scores.
%
% R: [MxA] Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (testcs = T*P'), this is the matrix of loadings 
%   P. For PLS (Y = testcs*W*inv(P'*W)*Q), this matrix is W*inv(P'*W). For the 
%   original space (default) the identity matrix is used.  
%
% Q: [MxA] Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (testcs = T*P'), this is the matrix of 
%   loadings P. For PLS (Y = testcs*W*inv(P'*W)*Q), this matrix is also P. 
%   For the original space the identity matrix is used. Q=R is used by 
%   default. 
%
%
% OUTPUTS:
%
% Dst: [Nx1] D-statistic or Hotelling T2
%
% Qst: [Nx1] Q-statistic
%
%
%
% EXAMPLE OF USE: PCA-based MSPC on NOC test data and anomalies.
%
% n_obs = 100;
% n_vars = 10;
% n_PCs = 1;
% X = simuleMV(n_obs,n_vars,6);
% [Xcs, m, sc] = preprocess2D(X,2);
% 
% pcs = 1:n_PCs;
% [p,t] = pca_pp(Xcs,pcs);
% e = Xcs - t*p';
% UCLd = hot_lim(n_PCs,n_obs,0.05,2);
% UCLq = spe_lim(e,0.01);
% 
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,6,corr(X)*(n_obst-1)/(n_obs-1));
% test(6:10,:) = 3*test(6:10,:);
% testcs = preprocess2Dapp(test,m,sc);
% 
% [Dst,Qst] = mspc(testcs,inv(cov(t)),p);
% 
% plot_scatter([Dst,Qst],[],[ones(5,1);2*ones(5,1)],{'D-st','Q-st'},{UCLd,UCLq}); 
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 22/Jan/17.
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
M = size(testcs, 2);
if nargin < 2, invCT = []; end;
A = size(invCT, 1);
if nargin < 3 || isempty(R), R = zeros(M,A); end;
if nargin < 4 || isempty(Q), Q = R; end;

% Validate dimensions of input data
assert (isequal(size(invCT), [A A]), 'Dimension Error: 2nd argument must be LVs-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(R), [M A]), 'Dimension Error: 3rd argument must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Q), [M A]), 'Dimension Error: 4th argument must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);


%% Main code

t = testcs * R;
e = testcs - t * Q';

Dst = sum((t * invCT) .* t,2);
Qst = sum(e.^2,2);


    


        