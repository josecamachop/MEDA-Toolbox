
function meda_map = meda(XX,R,Q)

% Missing data methods for Exploratory Data Analysis (MEDA). The original
% paper is Chemometrics and Intelligent Laboratory Systems 103(1), 2010, pp.
% 8-18. This algorithm follows the suggested computation by Arteaga in his
% technical report "A Note on MEDA", attached to the toolbox, which
% makes use of the covariance matrices.
%
% meda_map = meda(XX,R)   % minimum call
% meda_map = meda(XX,R,Q) % complete call
%
%
% INPUTS:
%
% XX: [MxM] cross-product X'*X 
%
% R: [MxA] Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (X = T*P'), this is the matrix of loadings P. 
%   For PLS (X = App*W*inv(P'*W)*Q'), this matrix is W*inv(P'*W). For the 
%   original space (default) the identity matrix is used. 
%
% Q: [MxA] Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (X = T*P'), this is the matrix of loadings 
%   P. For PLS (Y = X*W*inv(P'*W)*Q), this matrix is also P. For the 
%   original space the identity matrix is used. Q=R is used by default. 
%
%
% OUTPUTS:
%
% meda_map: [MxM] MEDA matrix.
%
%
% EXAMPLE OF USE: MEDA on PCA
%
% X = simuleMV(20,10,8);
% Xcs = preprocess2D(X,2);
% pcs = 1:3;
% p = pca_pp(Xcs,pcs);
%
% meda_map = meda(Xcs'*Xcs,p,p);
%
% [meda_map,ord] = seriation(meda_map);
% plot_map(meda_map,ord);
%
%
% EXAMPLE OF USE: MEDA on PLS
%
% X = simuleMV(20,10,8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Xcs = preprocess2D(X,2);
% Ycs = preprocess2D(Y,2);
% lvs = 1:10;
% [beta,W,P,Q,R] = kernel_pls(Xcs'*Xcs,Xcs'*Ycs,lvs);
%
% meda_map = meda(Xcs'*Xcs,R,P);
%
% [meda_map,ord] = seriation(meda_map);
% plot_map(meda_map,ord);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Mar/16.
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
M = size(XX, 1);
A = size(R, 2);
if nargin < 3 || isempty(Q), Q = R; end;

% Validate dimensions of input data
assert (isequal(size(XX), [M M]), 'Dimension Error: 1st argument must be M-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(R), [M A]), 'Dimension Error: 2nd argument must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Q), [M A]), 'Dimension Error: 3rd argument must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);


%% Main code

XXA = XX*R*Q';
dXX = diag(XX)*diag(XX)';

meda_map = (2*XX.*abs(XXA) - XXA.*abs(XXA))./dXX;

meda_map(find(~dXX))=0;


    


        