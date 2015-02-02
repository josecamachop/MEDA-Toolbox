
function [meda_map,meda_dis] = meda(XX,R,Q,thres)

% Missing data methods for Exploratory Data Analysis (MEDA). The original
% paper is Chemometrics and Intelligent Laboratory Systems 103(1), 2010, pp.
% 8-18. This algorithm follows the suggested computation by Arteaga in his
% technical report "A Note on MEDA", attached to the toolbox, which
% makes use of the covariance matrices.
%
% [meda_map,meda_dis] = meda(XX,R)   % minimum call
% [meda_map,meda_dis] = meda(XX,R,Q,thres) % complete call
%
%
% INPUTS:
%
% XX: (MxM) cross-product matrix (XX = xp'*xp) from the preprocessed 
%   billinear data set under analysis (xp). 
%
% R: (MxLVs) Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (X = T*P'), this is the matrix of loadings P. 
%   For PLS (Y = X*W*inv(P'*W)*Q), this matrix is W*inv(P'*W).  
%
% Q: (MxLVs) Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (X = T*P'), this is the matrix of loadings 
%   P. For PLS (Y = X*W*inv(P'*W)*Q), this matrix is also P. (Q=R by
%   default)
%
% thres: (1x1) threshold for the discretized MEDA matrix (0.1 by default)
%
%
% OUTPUTS:
%
% meda_map: (MxM) MEDA matrix.
%
% meda_dis: (MxM) discretized MEDA matrix.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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

%% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
s = size(XX);
if size(R,2)<1 || size(R,1)~=s(2), error('Error in the dimension of the arguments.'); end;
if nargin < 3, Q = R; end;
if size(R) ~= size(Q), error('Error in the dimension of the arguments.'); end;
if nargin < 4, thres=0.1; end; 

%% Main code

XXA = XX*R*Q';
dXX = diag(XX)*diag(XX)';

meda_map = (2*XX.*abs(XXA) - XXA.*abs(XXA))./dXX;

meda_map(find(~dXX))=0;

meda_dis = meda_map;
ind=find(abs(meda_map)<=thres);
meda_dis(ind) = 0;

    


        