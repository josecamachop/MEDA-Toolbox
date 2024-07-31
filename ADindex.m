function [ind,diff] = ADindex(L,App,varargin)

% ADICOV similarity index. This index allows us to estiamte the distance
% between a matrix and its corresponding approximation with ADICOV. 
% Reference: Camacho, J., Padilla, P., Díaz-Verdejo, J., Smith, K., Lovett, 
% D. Least-squares approximation of a space distribution for a given 
% covariance and latent sub-space. Chemometrics and Intelligent Laboratory 
% Systems, 2011, 105 (2): 171-180.
%
% ind = ADindex(L,App) % minimum call
% [ind,diff] = ADindex(L,App,'InSubspace',R,'Index',index) % complete call
%
%
% See also: ADICOV, simuleMV, MSPC_ADICOV
%
%
% INPUTS:
%
% L: [NxM] original data set.
%
% App: [NxM] data set approximated by ADICOV.
%
%
% Optional INPUTS:
%
% 'InSubspace': [MxA] proyection matrix.
%
% 'Index': [1x1] MSPC index definition
%       0: ADICOV similarity index according to Chemometrics and Intelligent 
%           Laboratory Systems 105, 2011, pp. 171-180 (default)
%       1: Modified index 
%
%
% OUTPUTS:
%
% ind: [Nx1] similarity index.
%
% diff: [NxM] difference betwen original and approximation in the
%   proyection space.
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   We approximate the distribution of L for a fixed cross-product XX, both
%   simulated at random.
%
% X = randn(100,2);
% XX = X'*X;
% L = randn(100,2); L(1,:) = 10* L(1,:); % data with an outlier
% appL = ADICOV(XX,L,2);
% 
% ind = ADindex(L,appL)
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
% ind = ADindex(L,appL,'InSubspace',model.loads)
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

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(L, 1);
M = size(L, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'InSubspace',eye(M));  
addParameter(p,'Index',1);          
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
R = p.Results.InSubspace;
index = p.Results.Index;

% Validate dimensions of input data
assert (isequal(size(L), [N M]), 'Dimension Error: ''L'' must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(App), [N M]), 'Dimension Error: ''App'' argument must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(R,1), M), 'Dimension Error: ''InSubspace'' must be M-by-PCs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(index), [1 1]), 'Dimension Error: ''Index'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); 

% Validate values of input data
assert (index==0 || index==1, 'Value Error: ''Index'' must be 0 or 1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

ux = L*R;
uy = App*R;
diff = ux-uy;

if index==1, diff(find(diff(:)<0)) = 0; end

ind = sum(diff.^2,2)/size(diff,2);
 
ind = sum(ind)/size(diff,1);