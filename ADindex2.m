function ind = ADindex2(L,A,R)

% Modified ADICOV similarity index based on division
%
% ind = ADindex2(L,A,R) % complete call
%
%
% INPUTS:
%
% L: {NxM} original data set.
%
% A: {NxM} data set approximated by ADICOV.
%
% R: (Mxpcs) proyection matrix.
%
%
% OUTPUTS:
%
% ind: (Nx1) similarity index.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 04/Sep/15
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


%% Main code

ux = L*R;
uy = A*R;
r = abs(ux) - abs(uy);
r(find(r(:)<0)) = 0;

ind = sum(r.^2,2)/size(r,2);
 
ind = sum(ind)/size(r,1);