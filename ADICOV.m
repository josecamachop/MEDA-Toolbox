function A = ADICOV(XX,L,rC,R,Q,multn)

% Approximation of a DIstribution for a given COVariance (ADICOV). The original
% paper is Chemometrics and Intelligent Laboratory Systems 105(2), 2011, pp.
% 171-180.
%
% A = ADICOV(XX,L,rC) % minimum call
% A = ADICOV(XX,L,rC,R,Q,multn) % complete call
%
%
% INPUTS:
%
% XX: (MxM) cross-product or covariance matrix to approximate. 
%   M has to be at least rC.
%
% L: (NxM) data set with the distribution of the observations to approximate.
%
% rC: (1x1) number of eigenvectors-eigenvalues of CV which are maintained 
%   in the approaximation.
%
% R: (MxLVs) Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (A = T*P'), this is the matrix of loadings P. 
%   For PLS (Y = A*W*inv(P'*W)*Q), this matrix is W*inv(P'*W). For the 
%   approximation in the original space (default) the identity matrix is
%   used. The number of director vectors in R (LVs) should be at least rC.
%
% Q: (MxLVs) Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (A = T*P'), this is the matrix of loadings 
%   Q. For PLS (Y = A*W*inv(P'*W)*Q), this matrix is also P. For the 
%   approximation in the original space the identity matrix is used. (Q=R 
%   by default) The number of director vectors in Q (LVs) should be at 
%   least rC.
%
% multn: (Nx1) multiplicity of the observations (ones by default)
%
%
% OUTPUTS:
%
% A: (NxM) approximation matrix.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 26/Sep/11.
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  José Camacho Páez
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

if nargin < 3, error('Error in the number of arguments.'); end;
s = size(XX,1);
if nargin < 4, R = eye(s); end;
if size(R,2)<rC || size(R,1)~=s, error('Error in the dimension of the arguments.'); end;
if nargin < 5, Q = R; end;
if size(R) ~= size(Q), error('Error in the dimension of the arguments.'); end;
if nargin < 6, multn = ones(size(L,1),1); end;

%% Main code

Gx = R(:,1:rC)'*XX*R(:,1:rC);
r = rank(Gx);
[Vx,Dx] = eig(Gx);
dd = abs(real(diag(Dx)));
[dd,indv] = sort(dd);
Vx = Vx(:,indv);
Sx = diag(sqrt(dd));
vind = size(Sx,1):-1:size(Sx,1)-r+1; 
Cmult = (sqrt(multn)*ones(1,s)).*L;
Tl = Cmult*R(:,1:rC);
M = Tl*(Sx(vind,vind)*Vx(:,vind)')';
[Um,Sm,Vm] = svd(M,0);
dime = min(size(Um,2),size(Vm,2));
Ua = Um(:,1:dime)*Vm(:,1:dime)';
Ta = Ua*Sx(vind,vind)*Vx(:,vind)';
A = (Ta*Q(:,1:rC)')./(sqrt(multn)*ones(1,s));
