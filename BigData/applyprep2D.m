function xp = applyprep2D(x,average,scale,weight)

% Apply preprocessing.
%
% xp = applyprep2D(x,average) % minumim call
% xp = applyprep2D(x,average,scale,weight) % complete call
%
%
% INPUTS:
%
% x: (NxM) Two-way data matrix, N(observations) x M(variables)
%
% average: (1xM) average to subtract.
%
% scale: (1xM) scale to divide.
%
% weight: (1xM) weight to multiply. Vector of ones by default.
%
%
% OUTPUTS:
%
% xp: (NxM) preprocessed data.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Mar/15.
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
    
%% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if nargin < 3, scale = ones(1,s(2)); end;
if nargin < 4, weight = ones(1,s(2)); end;
    
% Computation

xc = x - ones(s(1),1)*average; 
xp = xc./(ones(s(1),1)*scale);
xp = xp.*(ones(s(1),1)*weight);