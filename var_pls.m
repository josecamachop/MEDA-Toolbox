
function [y_var,t_var] = var_pls(cal,y,maxlvs,prepx,prepy,opt)

% Variability captured in terms of the number of LVs.
%
% var_pls(cal,y,maxlvs) % minimum call
% var_pls(cal,y,maxlvs,prepx,prepy,opt) %complete call
%
%
% INPUTS:
%
% cal: (LxM) billinear data set of predictor variables
%
% y: (LxO) billinear data set of predicted variables
%
% maxlvs: (1x1) Latent Variables considered (e.g. maxlvs = 2 selects the
%   first two lvs)
%
% prepx: (1x1) preprocesing of the x-block
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)  
%
% prepy: (1x1) preprocesing of the y-block
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
% y_var: ((maxlvs+1)x1) Percentage of captured variance of Y.
%
% t_var: ((maxlvs+1)x1) Percentage of captured variance of the scores.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
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

if nargin < 3, error('Error in the number of arguments.'); end;
s = size(cal);
if s(1) < 1 || s(2) < 1 || ndims(cal)~=2, error('Error in the dimension of the arguments.'); end;
if nargin < 4, prepx = 2; end;
if nargin < 5, prepy = 2; end; 
if nargin < 6, opt = 1; end;

%% Main code

cal_prep = preprocess2D(cal,prepx); 
y_prep = preprocess2D(y,prepy); 

XX = cal_prep'*cal_prep;
XY = cal_prep'*y_prep;
[beta,W,P,Q,R] = kernel_pls(XX,XY,maxlvs);

totalVt = sum(sum(cal_prep.^2));
t_var = ones(maxlvs+1,1);
totalVy = sum(sum(y_prep.^2));
y_var = ones(maxlvs+1,1);
for i=1:maxlvs,
    t_var(i+1) = t_var(i+1) - sum(eig(R(:,1:i)'*XX*R(:,1:i)))/totalVt;
    y_var(i+1) = y_var(i+1) - sum(eig(Q(:,1:i)*R(:,1:i)'*XX*R(:,1:i)*Q(:,1:i)'))/totalVy;
end
    
%% Show results

if opt == 1,
    fig_h = plot_vec(y_var,num2str((0:maxlvs)')','% Residual Variance',[],1);
    fig_h = plot_vec(t_var,num2str((0:maxlvs)')','% Residual Variance',[],1,'r--',fig_h,{'Y','Scores'});
end

        