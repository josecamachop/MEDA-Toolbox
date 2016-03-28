
function [x_var,cumpress] = var_pca(cal,maxpcs,prep,opt)

% Variability captured in terms of the number of PCs. It includes the ckf
% algorithm.
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
%       1: % Residual Variance in X 
%       2: % Residual Variance in X and ckf (default)
%
%
% OUTPUTS:
%
% x_var: ((maxpcs+1)x1) Percentage of captured variance of X.
%
% cumpress: ((maxpcs+1)x1) ckf curve.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 28/Mar/16.
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
if nargin < 4, opt = 2; end;

if opt >2 opt = 2; end;

%% Main code

cal_prep = preprocess2D(cal,prep); 

[p,T] = pca_pp(cal_prep,1:maxpcs);

totalVx = sum(sum(cal_prep.^2));
x_var = ones(maxpcs+1,1);
for i=1:maxpcs,
    x_var(i+1) = x_var(i+1) - sum(eig(T(:,1:i)'*T(:,1:i)))/totalVx;
end
    
if opt ==2,

    cumpress = zeros(maxpcs+1,1);
    press = zeros(maxpcs+1,s(2));

    xcs = cal_prep;
    
    if ~prep,
        avs_prep=ones(s(1),1)*mean(xcs);
    else
        avs_prep=zeros(s);
    end
    
    [p,t_est] = pca_pp(xcs,maxpcs);
    
    for i=0:maxpcs,
        
        if i > 0, % PCA Modelling
            
            p2 = p(:,1:min(i,end));
            srec = t_est(:,1:min(i,end))*p2';
            erec = xcs - srec;
            term3_p = erec;
            term1_p = (xcs-avs_prep).*(ones(s(1),1)*(sum(p2.*p2,2))');
            
        else % Modelling with the average
            term1_p = zeros(size(xcs));
            term3_p = xcs;
        end
        
        term1 = sum(term1_p.^2,1);
        term2 = sum(2*term1_p.*term3_p,1);
        term3 = sum(term3_p.^2,1);
        
        press(i+1,:) = sum([term1;term2;term3]);
        
        cumpress(i+1) = sum(press(i+1,:));
    end
end
    
%% Show results

if opt >= 1,
    fig_h = plot_vec(x_var,num2str((0:maxpcs)')','% Residual Variance',[],1);
    if opt == 2,
        fig_h = plot_vec(cumpress/cumpress(1),num2str((0:maxpcs)')','% Residual Variance',[],1,'r--',fig_h,{'X','ckf'});
    end
end

        