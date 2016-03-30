
function [y_var,t_var] = var_pls(x,y,lvs,prepx,prepy,opt)

% Variability captured in terms of the number of LVs.
%
% var_pls(x,y,maxlvs) % minimum call
% var_pls(x,y,maxlvs,prepx,prepy,opt) %complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
% lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(x)
%
% prepx: [1x1] preprocesing of the x-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)  
%
% prepy: [1x1] preprocesing of the y-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)   
%
% opt: [1x1] options for data plotting.
%       0: no plots.
%       1: Residual Variance in Y 
%       2: Residual Variance in Y and Scores (default)
%
%
% OUTPUTS:
%
% y_var: [Ax1] Percentage of captured variance of Y.
%
% t_var: [Ax1] Percentage of captured variance of the scores.
%
%
% EXAMPLE OF USE: Random data
%
% X = real(ADICOV(randn(10,10).^9,randn(100,10),10));
% Y = randn(100,2) + X(:,1:2);
% lvs = 0:10;
% x_var = var_pls(X,Y,lvs);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 30/Mar/16.
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine.name);
N = size(x, 1);
M = size(x, 2);
O = size(y, 2);
if nargin < 3 || isempty(lvs), lvs = 0:rank(x); end;
A = length(lvs);
if nargin < 4 || isempty(prepx), prepx = 2; end;
if nargin < 5 || isempty(prepy), prepy = 2; end;
if nargin < 6 || isempty(opt), opt = 2; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Validate dimensions of input data
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);

% Preprocessing
lvs = unique(lvs);

% Validate values of input data
assert (isempty(find(lvs<0)), 'Value Error: 3rd argument must not contain negative values. Type ''help %s'' for more info.', routine.name);
assert (isequal(fix(lvs), lvs), 'Value Error: 3rd argumentmust contain integers. Type ''help %s'' for more info.', routine.name);


%% Main code

xcs = preprocess2D(x,prepx); 
ycs = preprocess2D(y,prepy); 

[beta,W,P,Q,R] = kernel_pls(xcs'*xcs,xcs'*ycs,1:lvs);

totalVt = sum(sum(xcs.^2));
t_var = ones(length(lvs),1);
totalVy = sum(sum(ycs.^2));
y_var = ones(length(lvs),1);
for i=1:length(lvs),
    t_var(i) = t_var(i) - sum(eig(R(:,1:lvs(i))'*xcs'*xcs*R(:,1:lvs(i))))/totalVt;
    y_var(i) = y_var(i) - sum(eig(Q(:,1:lvs(i))*R(:,1:lvs(i))'*xcs'*xcs*R(:,1:lvs(i))*Q(:,1:lvs(i))'))/totalVy;
end
    
%% Show results
        
if opt,
    switch opt,
        case 1
            plot_vec(y_var,lvs,[],{'% Residual Variance in Y','PCs'},[],1);
        otherwise
            plot_vec([y_var t_var],lvs,[],{'% Residual Variance','PCs'},[],1,{'Y','Scores'});
            legend('show');
    end
end
