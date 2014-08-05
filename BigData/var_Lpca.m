function x_var = var_Lpca(Lmodel,maxpcs,opt)

% Variability captured in terms of the number of PCs.
%
% var_Lpca(Lmodel,maxpcs) % minimum call
% var_Lpca(Lmodel,maxpcs,opt) %complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%
% maxpcs: (1x1) Principal Components considered (e.g. maxpcs = 2 selects the
%   first two PCs)
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: bar plot (default)
%
%
% OUTPUTS:
%
% x_var: ((maxpcs+1)x1) Percentage of captured variance of X.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 08/May/13.
%
% Copyright (C) 2014  José Camacho Páez
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
if nargin < 3, opt = 1; end;

%% Main code

Lmodel.lv = maxpcs;   
P = Lpca(Lmodel);

totalVx = sum(eig(Lmodel.XX));
x_var = ones(maxpcs+1,1);
for i=1:maxpcs,
    x_var(i+1) = x_var(i+1) - sum(eig(P(:,1:i)'*Lmodel.XX*P(:,1:i)))/totalVx;
end
    
%% Show results

if opt == 1,
    fig_h = plot_vec(x_var,num2str((0:maxpcs)')','% Residual Variance',[],1);
end
