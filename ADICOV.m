function App = ADICOV(XX,L,Neig,R,Q,multn)

% Approximation of a DIstribution for a given COVariance (ADICOV). The original
% paper is Chemometrics and Intelligent Laboratory Systems 105(2), 2011, pp.
% 171-180.
%
% App = ADICOV(XX,L,Neig) % minimum call
% App = ADICOV(XX,L,Neig,R,Q,multn) % complete call
%
%
% INPUTS:
%
% XX: [MxM] cross-product or covariance matrix to approximate. M has to be 
%   at least Neig
%
% L: [NxM] data set with the distribution of the observations to approximate
%
% Neig: [1x1] number of eigenvectors-eigenvalues of XX which are maintained 
%   in the approximation.
%
% R: [MxA] Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (App = T*P'), this is the matrix of loadings P. 
%   For PLS (Y = App*W*inv(P'*W)*Q'), this matrix is W*inv(P'*W). For the 
%   approximation in the original space (default) the identity matrix is
%   used. The number of director vectors in R (LVs) should be at least Neig 
%
% Q: [MxA] Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (App = T*P'), this is the matrix of loadings 
%   P. For PLS (Y = App*W*inv(P'*W)*Q), this matrix is also P. For the 
%   approximation in the original space the identity matrix is used. Q=R is 
%   used by default. The number of director vectors in Q (LVs) should be at 
%   least Neig
%
% multn: [Nx1] multiplicity of the observations (ones by default)
%
%
% OUTPUTS:
%
% App: [NxM] approximation matrix.
%
%
% EXAMPLE OF USE: To obtain a matrix 100x10 with random covariance matrix, 
% use the following call:
%
% X = real(ADICOV(randn(10,10),randn(100,10),10));
%
% The call to real is necessary to avoid imaginary part in the results
% since we did not constrain XX to be positive definite. 
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 05/May/16.
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
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
M = size(XX, 1);
N = size(L, 1);
if nargin < 4 || isempty(R), R = eye(M); end;
A = size(R, 2);
if nargin < 5 || isempty(Q), Q = R; end;
if nargin < 6 || isempty(multn),  multn = ones(size(L,1),1); end;

% Convert row arrays to column arrays
if size(multn,1) == 1,     multn = multn'; end;

if Neig>rank(XX)
    Neig = rank(XX);
end

% Validate dimensions of input data
assert (isequal(size(XX), [M M]), 'Dimension Error: 1st argument must be M-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(L), [N M]), 'Dimension Error: 2nd argument must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Neig), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(R), [M A]), 'Dimension Error: 4th argument must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Q), [M A]), 'Dimension Error: 5th argument must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(multn), [N 1]), 'Dimension Error: 6th argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (Neig>0, 'Value Error: 3rd argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (Neig<=rank(XX), 'Value Error: 3rd argument must not be above the rank of XX. Type ''help %s'' for more info.', routine(1).name);
assert (Neig<=A, 'Value Error: 3rd argument must not be above the number of columns of 4th and 5th arguments. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(Neig), Neig), 'Value Error: 3rd argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(multn<0)), 'Value Error: 6th argument must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(multn), multn), 'Value Error: 6th argument must contain integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

Gx = R(:,1:Neig)'*XX*R(:,1:Neig);
r = rank(Gx);
[Vx,Dx] = eig(Gx);
dd = abs(real(diag(Dx)));
[dd,indv] = sort(dd);
Vx = Vx(:,indv);
Sx = diag(sqrt(dd));
vind = size(Sx,1):-1:size(Sx,1)-r+1; 
Cmult = (sqrt(multn)*ones(1,M)).*L;
Tl = Cmult*R(:,1:Neig);
[Um,Sm,Vm] = svd(Tl*(Sx(vind,vind)*Vx(:,vind)')',0);
dime = min(size(Um,2),size(Vm,2));
Ua = Um(:,1:dime)*Vm(:,1:dime)';
Ta = Ua*Sx(vind,vind)*Vx(:,vind)';
App = (Ta*Q(:,1:Neig)')./(sqrt(multn)*ones(1,M));
