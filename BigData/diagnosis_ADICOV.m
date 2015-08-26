function [diff,ord]= diagnosis_ADICOV(Lmodel,pcs,file,opt)

% Diagnosis using ADICOV.
%
% Diagnosis_ADICOV(Lmodel,pcs,file) % minimum call
% Diagnosis_ADICOV(Lmodel,pcs,file,opt) % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) Lmodel for monitoring.
%
% pcs: (1x1) number of components.
%
% file: {str} name of the file with x for diagnosis.
%
% opt: (1x1) options for data plotting.
%       0: no plots
%       1: Differential map (default)
%
%
% OUTPUTS:
%
% diff: (MxM) Differential map
%
% ord: (Mx1) ordering of variables in the map
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Aug/15.
%
% Copyright (C) 2015  University of Granada, Granada
% Copyright (C) 2015  Jose Camacho Paez
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

if nargin < 3, error('Error in the number of arguments.'); end;
if nargin < 4, opt = 1; end;

Lm = Lmodel;
load(file,'x','obs_l');
xs = applyprep2D(x,Lm.av,Lm.sc,Lm.weight);

opt.plot = 0;
opt.seriated = 1;
opt.discard = 0;
[meda_map_L,kk,ord] = meda_Lpca(Lm,pcs,0.1,opt);
opt.seriated = 0;
meda_map_x = meda_pca(xs(:,ord),pcs,0,0.1,opt);

diff = meda_map_x - meda_map_L;

if opt, 
    plot_map(diff, var_l(ord));
end

