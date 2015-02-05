function [y_var,t_var] = var_Lpls(Lmodel,maxlvs,opt)

% Variability captured in terms of the number of LVs.
%
% var_Lpls(Lmodel,maxlvs) % minimum call
% var_Lpls(Lmodel,maxlvs,opt) %complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PLS
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.XY: (MxL) cross-product matrix between the x-block and the
%           y-block.
%       Lmodel.YY: (LxL) Y-block cross-product matrix.
%
% maxlvs: (1x1) Latent Variables considered (e.g. maxlvs = 2 selects the
%   first two lvs)
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: bar plot (default)
%
%
% OUTPUTS:
%
% y_var: ((maxlvs+1)x1) Percentage of captured variance of Y.
%
% t_var: ((maxlvs+1)x1) Percentage of captured variance of the scores.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 08/May/13.
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
if nargin < 3, opt = 1; end;

%% Main code

Lmodel.lv = maxlvs;   
[beta,W,P,Q,R] = Lpls(Lmodel);

totalVt = sum(eig(Lmodel.XX));
t_var = ones(maxlvs+1,1);
totalVy = sum(eig(Lmodel.YY));
y_var = ones(maxlvs+1,1);
for i=1:maxlvs,
    t_var(i+1) = t_var(i+1) - sum(eig(R(:,1:i)'*Lmodel.XX*R(:,1:i)))/totalVt;
    y_var(i+1) = y_var(i+1) - sum(eig(Q(:,1:i)*R(:,1:i)'*Lmodel.XX*R(:,1:i)*Q(:,1:i)'))/totalVy;
end
    
%% Show results

if opt == 1,
    fig_h = plot_vec(y_var,num2str((0:maxlvs)')','% Residual Variance',[],1);
    fig_h = plot_vec(t_var,num2str((0:maxlvs)')','% Residual Variance',[],1,'r--',fig_h,{'Y','Scores'});
end
