
function x_var = var_pca(cal,maxpcs,prep,opt)

% Variability captured in terms of the number of PCs.
%
% var_pca(cal,maxpcs) % minimum call
% var_pca(cal,maxpcs,prep,opt) %complete call
%
%
% INPUTS:
%
% cal: (LxM) billinear data set
%
% maxpcs: (1x1) Principal Components considered (e.g. maxpcs = 2 selects the
%   first two PCs)
%
% prep: (1x1) preprocesing 
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)  
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
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
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

%% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
s = size(cal);
if s(1) < 1 || s(2) < 1 || ndims(cal)~=2, error('Error in the dimension of the arguments.'); end;
if nargin < 3, prep = 2; end;
if nargin < 4, opt = 1; end;

%% Main code

cal_prep = preprocess2D(cal,prep); 

[p,T] = pca_pp(cal_prep,maxpcs);

totalVx = sum(sum(cal_prep.^2));
x_var = ones(maxpcs+1,1);
for i=1:maxpcs,
    x_var(i+1) = x_var(i+1) - sum(eig(T(:,1:i)'*T(:,1:i)))/totalVx;
end
    
%% Show results

if opt == 1,
    fig_h = plot_vec(x_var,num2str((0:maxpcs)')','% Residual Variance',[],1);
end

        