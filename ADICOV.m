function App = ADICOV(XX,L,Neig,varargin)

% Approximation of a DIstribution for a given COVariance (ADICOV). App is a 
% matrix that approximates L but with covariance equal to XX. Reference: 
% Camacho, J., Padilla, P., Díaz-Verdejo, J., Smith, K., Lovett, D. 
% Least-squares approximation of a space distribution for a given 
% covariance and latent sub-space. Chemometrics and Intelligent Laboratory 
% Systems, 2011, 105 (2): 171-180.
%
% App = ADICOV(XX,L,Neig) % minimum call
% App = ADICOV(XX,L,Neig,'InSubspace',R,'OutSubspace',Q,'Multiplicity',multn) % complete call
%
%
% See also: ADindex, simuleMV, MSPC_ADICOV
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
% Optional INPUTS:
%
% 'InSubspace': [MxA] Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (App = T*P'), this is the matrix of loadings P. 
%   For PLS (Y = App*W*inv(P'*W)*Q'), this matrix is W*inv(P'*W). For the 
%   approximation in the original space (default) the identity matrix is
%   used. The number of director vectors in R (LVs) should be at least Neig 
%
% 'OutSubspace': [MxA] Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (App = T*P'), this is the matrix of loadings 
%   P. For PLS (Y = App*W*inv(P'*W)*Q), this matrix is also P. For the 
%   approximation in the original space the identity matrix is used. Q=R is 
%   used by default. The number of director vectors in Q (LVs) should be at 
%   least Neig
%
% 'Multiplicity': [Nx1] multiplicity of the observations (ones by default)
%
%
% OUTPUTS:
%
% App: [NxM] approximation matrix.
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   We approximate the distribution of L for a fixed cross-product XX, both
%   simulated at random.
%
% X = randn(100,2);
% XX = X'*X
% L = randn(100,2); L(1,:) = 10* L(1,:); % data with an outlier
% appL = ADICOV(XX,L,2);
% 
% appL'*appL % same covariance as X
% plot_scatter(appL); % similar distribution as L
% plot_scatter(L); 
%
%
% EXAMPLE OF USE: (copy and paste the code in the command line)
%   Similar example as before in a PCA subspace
%
% X = randn(100,10);
% XX = X'*X;
% L = randn(100,10); L(1,:) = 10* L(1,:); % data with an outlier
% model = pca_eig(L,'Pcs',1:2);
% appL = ADICOV(XX,L,2,'InSubspace',model.loads);
% 
% model.loads'*X'*X*model.loads % same covariance as the scores of X
% model.loads'*appL'*appL*model.loads
% 
% plot_scatter(appL*model.loads); % similar distribution as the scores of L
% plot_scatter(L*model.loads); 
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 15/Apr/2024
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

% Check minimum call
routine=dbstack;
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
M = size(XX, 1);
N = size(L, 1);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'InSubspace',eye(M));  
addParameter(p,'OutSubspace',[]);
addParameter(p,'Multiplicity',ones(size(L,1),1));            
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
R = p.Results.InSubspace;
A = size(R, 2);
Q = p.Results.OutSubspace;
if isempty(Q), Q = R; end;
multn = p.Results.Multiplicity;

% Convert row arrays to column arrays
if size(multn,1) == 1,     multn = multn'; end;

if Neig>rank(XX)
    Neig = rank(XX);
end

% Validate dimensions of input data
assert (isequal(size(XX), [M M]), 'Dimension Error: ''XX'' must be M-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(L), [N M]), 'Dimension Error: ''L'' must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Neig), [1 1]), 'Dimension Error: ''Neig'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(R), [M A]), 'Dimension Error: ''InSubspace'' must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Q), [M A]), 'Dimension Error: ''OutSubspace'' must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(multn), [N 1]), 'Dimension Error: ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (Neig>0, 'Value Error: ''Neig'' must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (Neig<=rank(XX), 'Value Error: ''Neig'' must not be above the rank of XX. Type ''help %s'' for more info.', routine(1).name);
assert (Neig<=A, 'Value Error: ''Neig'' must not be above the number of columns of 4th and 5th arguments. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(Neig), Neig), 'Value Error: ''Neig'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(multn<0)), 'Value Error: ''Multiplicity'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(multn), multn), 'Value Error: ''Multiplicity'' must contain integers. Type ''help %s'' for more info.', routine(1).name);


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
App = real((Ta*Q(:,1:Neig)')./(sqrt(multn)*ones(1,M)));
