
function [x_var,cumpress] = var_pca(x,pcs,prep,opt)

% Variability captured in terms of the number of PCs. It includes the ckf
% algorithm.
%
% var_pca(x,pcs) % minimum call
% var_pca(x,pcs,prep,opt) %complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:rank(x)
%
% prep: [1x1] preprocesing
%       0: no preprocessing 
%       1: mean-centering 
%       2: auto-scaling (default)  
%
% opt: [1x1] options for data plotting
%       0: no plots.
%       1: Residual Variance in X 
%       otherwise: Residual Variance in X and ckf (default)
%
%
% OUTPUTS:
%
% x_var: [Ax1] Percentage of captured variance of X.
%
% cumpress: [Ax1] ckf curve.
%
%
% EXAMPLE OF USE: Random data
%
% X = real(ADICOV(randn(10,10).^9,randn(100,10),10));
% pcs = 0:10;
% x_var = var_pca(X,pcs);
%
%
% codified by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 29/Mar/16.
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine.name);
N = size(x, 1);
M = size(x, 2);
if nargin < 2 || isempty(pcs), pcs = 0:rank(x); end;
A = length(pcs);
if nargin < 3 || isempty(prep), prep = 2; end;
if nargin < 4 || isempty(opt), opt = 2; end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 2nd argument must be 1-by-A. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);

% Preprocessing
pcs = unique(pcs);

% Validate values of input data
assert (isempty(find(pcs<0)), 'Value Error: 2nd argument must not contain negative values. Type ''help %s'' for more info.', routine.name);
assert (isequal(fix(pcs), pcs), 'Value Error: 2nd argumentmust contain integers. Type ''help %s'' for more info.', routine.name);


%% Main code

xcs = preprocess2D(x,prep); 

[p,T] = pca_pp(xcs,pcs);

totalVx = sum(sum(xcs.^2));
x_var = ones(length(pcs),1);

for i = pcs,
    x_var(i+1) = x_var(i+1) - sum(eig(T(:,1:i)'*T(:,1:i)))/totalVx;
end
    
if opt ==2,

    cumpress = zeros(length(pcs),1);
    press = zeros(length(pcs),M);
    
    if ~prep,
        avs_prep=ones(N,1)*mean(xcs);
    else
        avs_prep=zeros(N,M);
    end
    
    [p,t_est] = pca_pp(xcs,pcs);
    
    for i=pcs,
        
        if i > 0, % PCA Modelling
            
            p2 = p(:,1:min(i,end));
            srec = t_est(:,1:min(i,end))*p2';
            erec = xcs - srec;
            term3_p = erec;
            term1_p = (xcs-avs_prep).*(ones(N,1)*(sum(p2.*p2,2))');
            
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

if opt,
    switch opt,
        case 1
            plot_vec(x_var,pcs,[],{'% Residual Variance','PCs'},[],1);
        otherwise
            plot_vec([x_var cumpress/cumpress(1)],pcs,[],{'% Residual Variance & ckf','PCs'},[],1);
    end
end

        