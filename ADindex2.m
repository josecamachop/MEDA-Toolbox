function ind = ADindex2(L,App,R)

% Modified ADICOV similarity index based on division
%
% ind = ADindex(L,App,R) % complete call
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
%
% OUTPUTS:
%
% ind: [Nx1] similarity index.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 22/Mar/16
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin == 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine.name);
N = size(L, 1);
M = size(L, 2);

% Validate dimensions of input data
assert (isequal(size(L), [N M]), 'Dimension Error: 1st argument must be N-by-M. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(App), [N M]), 'Dimension Error: 2nd argument must be N-by-M. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(R,1), M), 'Dimension Error: 3rd argument must be M-by-PCs. Type ''help %s'' for more info.', routine.name);


%% Main code

ux = L*R;
uy = App*R;
r = abs(ux) - abs(uy);
r(find(r(:)<0)) = 0;

ind = sum(r.^2,2)/size(r,2);
 
ind = sum(ind)/size(r,1);