function Lmodel2 = selectVars(Lmodel,ind)

% Reordering/selection of variables in Lmodel.
%
% Lmodel2 = selectVars(Lmodel,ind)     % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) input Lmodel
%   model:
%       Lmodel.centr: [NxM] centroids of the clusters of observations.
%       Lmodel.vclass: [Mx1] class associated to each variable.
%       Lmodel.XX: [MxM] X-block cross-product matrix.
%       Lmodel.av: [1xM] sample average according to the preprocessing method.
%       Lmodel.sc: [1xM] sample scale according to the preprocessing method.
%       Lmodel.weight: [1xM] weight applied after the preprocessing method.
%       Lmodel.var_l: {Mx1} label of each variable.
%       Lmodel.mat: [MxA] projection matrix for distance computation.
%       Lmodel.XY: [MxL] sample cross-product matrix of X and Y.
%
% ind: [1xL]: ordering/selection of variables.
%
%
% OUTPUTS:
%
% Lmodel2: (struct Lmodel) model after variable selection/reordering.
%
%
% EXAMPLE OF USE: Random data:
%
% X = simuleMV(20,10,'LevelCorr',8);
% Lmodel = iniLmodel(X);
% Lmodel2 = selectVars(Lmodel, [10 1:3])
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 19/Nov/2024
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

% Validate values of input data
[ok, Lmodel] = checkLmodel(Lmodel);

%% Main code

Lmodel2 = Lmodel;
Lmodel2.vclass = Lmodel.vclass(ind);
Lmodel2.XX = Lmodel.XX(ind,ind);
Lmodel2.av = Lmodel.av(ind);
Lmodel2.sc = Lmodel.sc(ind);
Lmodel2.weight = Lmodel.weight(ind);
if isfield(Lmodel,'mat') && ~isempty(Lmodel.mat), Lmodel2.mat = Lmodel.mat(ind,:); end
if isfield(Lmodel,'XY') && ~isempty(Lmodel.XY),Lmodel2.XY = Lmodel.XY(ind,:); end
Lmodel2.centr = Lmodel.centr(:,ind);
Lmodel2.var_l = Lmodel.var_l(ind);
