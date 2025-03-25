function App = AFamily(XX,L,varargin)

% Family of approximation methods to generate multivariate data with control 
% over the structure of columns and rows using two-sided orthogonal Procrustes. 
% Reference:  
%
% App = AFamily(XX,L,neig) % minimum call
%
%
% See also: ADICOV, ADindex, simuleMV, MSPC_ADICOV
%
%
% INPUTS:
%
% XX: [MxM] cross-product or covariance matrix to approximate. M has to be 
%   at least neig
%
% L: [NxM] data set with the distribution of the observations to approximate
%
%
% Optional INPUTS: 
%
% 'neig': [1x1] number of eigenvectors-eigenvalues of XX which are maintained 
%   in the approximation.
%
% 'InSubspace': [MxA] Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (App = T*P'), this is the matrix of loadings P. 
%   For PLS (Y = App*W*inv(P'*W)*Q'), this matrix is W*inv(P'*W). For the 
%   approximation in the original space (default) the identity matrix is
%   used. The number of director vectors in R (LVs) should be at least neig 
%
% 'OutSubspace': [MxA] Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (App = T*P'), this is the matrix of loadings 
%   P. For PLS (Y = App*W*inv(P'*W)*Q), this matrix is also P. For the 
%   approximation in the original space the identity matrix is used. Q=R is 
%   used by default. The number of director vectors in Q (LVs) should be at 
%   least neig
%
% 'Multiplicity': [Nx1] multiplicity of the observations (ones by default)
%
% 'Method': [1x1] One of the followings:
%   - 'ADIEN': Approximation of a DIstribution for given EigeNvalues (ADIEN)
%   - 'ARDIEN': Approximation of a Row DIstribution for given EigeNvalues (ARDIEN)
%   - 'ARDICO': Approximation of a Row DIstribution for a given Covariance (ARDICO, by default)
%   - 'ACORDI': Approximation of a COvariance for a given Row Distribution (ACORDI)
%   - 'Continuum': Continuum ARDICO – ACORDI
%   
% 'Gamma': [1x1] tradeoff for continuum, from 0 (ARDICO) to 1 (ACORDI).
%   Value of 0.5 by default.
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
% appL = AFamily(XX,L);
% 
% appL'*appL % same covariance as X
% plot_scatter(appL); % similar distribution as L
% plot_scatter(L); 
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 15/Sep/2024
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
M = size(XX, 1);
N = size(L, 1);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'neig',M);
addParameter(p,'InSubspace',eye(M));  
addParameter(p,'OutSubspace',[]);
addParameter(p,'Multiplicity',ones(size(L,1),1)); 
addParameter(p,'Method','ARDICO');  
addParameter(p,'Gamma',0.5);             
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
R = p.Results.InSubspace;
A = size(R, 2);
Q = p.Results.OutSubspace;
if isempty(Q), Q = R; end;
multn = p.Results.Multiplicity;
method = p.Results.Method;
gamma = p.Results.Gamma;
neig = p.Results.neig;

% Convert row arrays to column arrays
if size(multn,1) == 1,     multn = multn'; end;


% Validate dimensions of input data
assert (isequal(size(XX), [M M]), 'Dimension Error: ''XX'' must be M-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(L), [N M]), 'Dimension Error: ''L'' must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(neig), [1 1]), 'Dimension Error: ''Neig'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(R), [M A]), 'Dimension Error: ''InSubspace'' must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Q), [M A]), 'Dimension Error: ''OutSubspace'' must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(multn), [N 1]), 'Dimension Error: ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(gamma), [1 1]), 'Dimension Error: ''Gamma'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (neig>0, 'Value Error: ''Neig'' must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (neig<=A, 'Value Error: ''Neig'' must not be above the number of columns of the subspaces. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(neig), neig), 'Value Error: ''Neig'' must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(multn<0)), 'Value Error: ''Multiplicity'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(multn), multn), 'Value Error: ''Multiplicity'' must contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (gamma>=0 & gamma<=1, 'Value Error: ''Gamma'' must be between 0 and 1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

Gx = R(:,1:neig)'*XX*R(:,1:neig);
[Vc,D2c] = eig(Gx);
dd = abs(real(diag(D2c)));
[dd,indv] = sort(dd,'descend');
Vc = Vc(:,indv);
Dc = diag(sqrt(dd));

Cmult = (sqrt(multn)*ones(1,M)).*L;
Tl = Cmult*R(:,1:neig);
[Uy,Sy,Vy] = svd(Tl,0);

switch lower(method)
    case {'adien'}
        disp('Method is ADIEN')
        App = sqrt(size(Uy,1)-1)*Uy*Dc*Vy';
        
    case 'ardien'
        disp('Method is ARDIEN')
        App = sqrt(size(Uy,1)-1)*Uy*Dc*Vy';
        
    case {'ardico'}
        disp('Method is ARDICO')
        App = sqrt(size(Uy,1)-1)*Uy*Dc*Vc';
        
    case 'acordi'
        disp('Method is ACORDI')
        App = Uy*Sy*Vc';
        
    case 'continuum'
        disp(sprintf('Method is Continuum for Gamma = %d',gamma))
        App = Uy*abs(gamma*Sy+(1-gamma)*sqrt(size(Uy,1)-1)*Dc)*Vc';
        
    otherwise
        disp('Unknown method.')
        App = L;
        
end
    
App = App./(sqrt(multn)*ones(1,M));
