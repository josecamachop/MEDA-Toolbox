function [ind,diff] = ADindex(L,App,R,index)

% ADICOV similarity index according to Chemometrics and Intelligent 
% Laboratory Systems 105, 2011, pp. 171-180
%
% ind = ADindex(L,App,R) % minimum call
% [ind,diff] = ADindex(L,App,R,index) % complete call
%
%
% INPUTS:
%
% L: [NxM] original data set.
%
% App: [NxM] data set approximated by ADICOV.
%
% R: [MxA] proyection matrix.
%
% index: (1x1) MSPC index definition
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
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 24/May/17
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
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(L, 1);
M = size(L, 2);
if nargin < 4 || isempty(index), index = 1; end;

% Validate dimensions of input data
assert (isequal(size(L), [N M]), 'Dimension Error: 1st argument must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(App), [N M]), 'Dimension Error: 2nd argument must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(R,1), M), 'Dimension Error: 3rd argument must be M-by-PCs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(index), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); 

% Validate values of input data
assert (index==0 || index==1, 'Value Error: 4th argument must be 0 or 1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

ux = L*R;
uy = App*R;
diff = ux-uy;

if index==1, diff(find(diff(:)<0)) = 0; end

ind = sum(diff.^2,2)/size(diff,2);
 
ind = sum(ind)/size(diff,1);