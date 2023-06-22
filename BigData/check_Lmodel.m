function [ok,Lmodel] = check_Lmodel(Lmodel)

% Check Lmodel integrity
%
% ok = check_Lmodel(Lmodel) % complete call
%
%
% INPUT:
%
% Lmodel.centr: [NxM] centroids of the clusters of observations.
%
% Lmodel.centrY: [NxL] responses of centroids of the clusters of
% observations.
%
% Lmodel.nc: [1x1] number of clusters in the model.
%
% Lmodel.multr: [ncx1] multiplicity of each cluster.
%
% Lmodel.class: [ncx1] class associated to each cluster.
%
% Lmodel.vclass: [Mx1] class associated to each variable.
%
% Lmodel.N: [1x1] number of effective observations in the model.
%
% Lmodel.type: [1x1] PCA (1) o PLS (2)
%
% Lmodel.update: [1x1] EWMA (1) or ITERATIVE (2)
%
% Lmodel.XX: [MxM] sample cross-product matrix of X.
%
% Lmodel.lvs: [1x1] number of latent variables (e.g. lvs = 1:2 selects the
%   first two LVs). By default, Lmodel.lvs = 1:rank(xcs)
%
% Lmodel.prep: [1x1] preprocesing of the data
%       0: no preprocessing.
%       1: mean centering (default) 
%       2: auto-scaling (centers and scales data so that each variable 
%           has variance 1)
%
% Lmodel.av: [1xM] sample average according to the preprocessing method.
%
% Lmodel.sc: [1xM] sample scale according to the preprocessing method.
%
% Lmodel.weight: [1xM] weight applied after the preprocessing method.
%
% Lmodel.updated: [ncx1] specifies whether a data point is new.
%
% Lmodel.obs_l: {ncx1} label of each cluster.
%
% Lmodel.var_l: {Mx1} label of each variable.
%
% Lmodel.mat: [MxA] projection matrix for distance computation.
%
% Lmodel.prepy: [1x1] preprocesing of the data
%       0: no preprocessing.
%       1: mean centering (default) 
%       2: auto-scaling (centers and scales data so that each variable 
%           has variance 1)
%
% Lmodel.avy: [1xL] sample average according to the preprocessing method.
%
% Lmodel.scy: [1xL] sample scale according to the preprocessing method.
%
% Lmodel.weighty: [1xL] weight applied after the preprocessing method.
%
% Lmodel.XY: [MxL] sample cross-product matrix of X and Y.
%
% Lmodel.YY: [LxL] sample cross-product matrix of Y.
%
% Lmodel.index_fich: {ncx1} file system with the original observations in
%   each cluster for ITERATIVE models.
%
% Lmodel.path: (str) path to the file system for ITERATIVE models.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 15/Jun/2023
%
% Copyright (C) 2023  University of Granada, Granada
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

routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
assert (isfield(Lmodel,'centr'), 'Content Error: Lmodel without centroids. Type ''help %s'' for more info.', routine(1).name);
M = size(Lmodel.centr, 2);

% Set default values
if ~isfield(Lmodel,'centrY') || isempty(Lmodel.centrY), Lmodel.centrY = []; end
L = size(Lmodel.centrY, 2);
if ~isfield(Lmodel,'type') || isempty(Lmodel.type), Lmodel.type = 1; end
if ~isfield(Lmodel,'update') || isempty(Lmodel.update), Lmodel.update = 2; end
if ~isfield(Lmodel,'nc') || isempty(Lmodel.nc) || Lmodel.nc==0, Lmodel.nc = max(100,size(Lmodel.centr,1)); end
if ~isfield(Lmodel,'N') || isempty(Lmodel.N) || Lmodel.N==0, Lmodel.N = size(Lmodel.centr,1); end
if ~isfield(Lmodel,'lvs') || isempty(Lmodel.lvs), Lmodel.lvs = 1:rank(Lmodel.XX); end
if ~isfield(Lmodel,'prep') || isempty(Lmodel.prep), Lmodel.prep = 0; end
if Lmodel.nc>0
    if ~isfield(Lmodel,'multr') || isempty(Lmodel.multr), Lmodel.multr = ones(size(Lmodel.centr,1),1); end
    if ~isfield(Lmodel,'class') || isempty(Lmodel.class), Lmodel.class = ones(size(Lmodel.centr,1),1); end
    if ~isfield(Lmodel,'updated') || isempty(Lmodel.updated), Lmodel.updated = ones(size(Lmodel.centr,1),1); end
    if ~isfield(Lmodel,'obs_l') || isempty(Lmodel.obs_l) 
        if size(Lmodel.centr,1)>1
            Lmodel.obs_l = cellstr(num2str((1:size(Lmodel.centr,1))')); 
        else
            Lmodel.obs_l = {};
        end
    end
else
    if ~isfield(Lmodel,'multr'), Lmodel.multr = []; end
    if ~isfield(Lmodel,'class'), Lmodel.class = []; end
    if ~isfield(Lmodel,'updated'), Lmodel.updated = []; end
    if ~isfield(Lmodel,'obs_l'), Lmodel.obs_l = {}; end
end
if ~isfield(Lmodel,'XX') || isempty(Lmodel.XX)
    if isempty(Lmodel.centr)
        Lmodel.XX = [];
    else
        X = (Lmodel.multr * ones(1,M)) .* Lmodel.centr;
        Lmodel.XX = X'*X;
    end
end
if M>0
    if ~isfield(Lmodel,'av') || isempty(Lmodel.av), Lmodel.av = zeros(1,M); end
    if ~isfield(Lmodel,'sc') || isempty(Lmodel.sc), Lmodel.sc = ones(1,M); end
    if ~isfield(Lmodel,'vclass') || isempty(Lmodel.vclass), Lmodel.vclass = ones(1,M); end
    if ~isfield(Lmodel,'weight') || isempty(Lmodel.weight), Lmodel.weight = ones(1,M); end
    if ~isfield(Lmodel,'var_l') || isempty(Lmodel.var_l), Lmodel.var_l = cellstr(num2str((1:M)')); end
else
    if ~isfield(Lmodel,'av'), Lmodel.av = []; end
    if ~isfield(Lmodel,'sc'), Lmodel.sc = []; end
    if ~isfield(Lmodel,'vclass'), Lmodel.vclass = []; end
    if ~isfield(Lmodel,'weight'), Lmodel.weight = []; end
    if ~isfield(Lmodel,'var_l'), Lmodel.var_l = {}; end
end
if ~isfield(Lmodel,'YY') || isempty(Lmodel.YY) 
    if isempty(Lmodel.centrY) 
        Lmodel.YY = [];
    else
        Y = (Lmodel.multr * ones(1,L)) .* Lmodel.centrY;
        Lmodel.YY = Y'*Y;
    end
end
if ~isfield(Lmodel,'XY') || isempty(Lmodel.XY)
    if isempty(Lmodel.centr) || isempty(Lmodel.centrY) 
        Lmodel.XY = [];
    else
        Lmodel.XY = X'*Y;
    end
end

% Convert column arrays to row arrays
if size(Lmodel.lvs,2) == 1, Lmodel.lvs = Lmodel.lvs'; end;
if size(Lmodel.av,2) == 1, Lmodel.av = Lmodel.av'; end;
if size(Lmodel.sc,2) == 1, Lmodel.sc = Lmodel.sc'; end;
if size(Lmodel.weight,2) == 1, Lmodel.weight = Lmodel.weight'; end;
if size(Lmodel.avy,2) == 1, Lmodel.avy = Lmodel.avy'; end;
if size(Lmodel.scy,2) == 1, Lmodel.scy = Lmodel.scy'; end;
if size(Lmodel.weighty,2) == 1, Lmodel.weighty = Lmodel.weighty'; end;

% Preprocessing
Lmodel.lvs = unique(Lmodel.lvs);
Lmodel.lvs(find(Lmodel.lvs==0)) = [];
A = length(Lmodel.lvs);

% Validate dimensions of input data
assert (isequal(size(Lmodel.XX), [M M]), 'Dimension Error: Lmodel.XX must be M-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Lmodel.lvs), [1 A]) | isequal(size(Lmodel.lvs), [0 1]), 'Dimension Error: Lmodel.lvs must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(Lmodel.type~=1 & Lmodel.type~=2 & Lmodel.type~=3)), 'Value Error: Lmodel.type must contain 1, 2 or 3. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(Lmodel.update~=1 & Lmodel.update~=2)), 'Value Error: Lmodel.update must contain 1 or 2. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(Lmodel.lvs<0)), 'Value Error: Lmodel.lvs must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(Lmodel.lvs), Lmodel.lvs), 'Value Error: Lmodel.lvs must contain integers. Type ''help %s'' for more info.', routine(1).name);

ok = true;