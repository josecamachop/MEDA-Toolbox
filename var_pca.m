
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
%   first two PCs). By default, pcs = 0:rank(x). The value for 0 PCs is
%   added at the begining if not specified.
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
% last modification: 05/Apr/16.
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
if nargin < 2 || isempty(pcs), pcs = 0:rank(x); end;
A = length(pcs);
if nargin < 3 || isempty(prep), prep = 2; end;
if nargin < 4 || isempty(opt), opt = 2; end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 2nd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
pcs = unique([0 pcs]);

% Validate values of input data
assert (isempty(find(pcs<0)), 'Value Error: 2nd argument must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(pcs), pcs), 'Value Error: 2nd argumentmust contain integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,prep); 

[P,T] = pca_pp(xcs,1:max(pcs));
pcs(find(pcs>size(P,2))) = [];

totalVx = sum(sum(xcs.^2));
x_var = ones(length(pcs),1);

for i = 1:length(pcs),
    x_var(i) = x_var(i) - sum(eig(T(:,1:pcs(i))'*T(:,1:pcs(i))))/totalVx;
end

cumpress = zeros(length(pcs),1);
if nargout>1 | (opt~=0 & opt~=1),
    for i = 1:length(pcs),
         c = ckf(xcs,T(:,1:pcs(i)),P(:,1:pcs(i)),0);
         cumpress(i) = c(end);
    end
end

    
%% Show results

if opt,
    vecn = pcs;
    lv = length(vecn(2:end));
    div = 1:lv;
    div = div(rem(lv,div)==0);
    stepN = div(find(div>lv/20,1));
    vec = 1:(lv+1);
    vec = vec(round(1:stepN:end));
    for i=vec,
        label{i} = num2str(vecn(i));
    end
    switch opt,
        case 1
            plot_vec(x_var,label,[],{'% Residual Variance','PCs'},[],1);
        otherwise
            plot_vec([x_var cumpress/cumpress(1)],label,[],{'% Residual Variance','PCs'},[],1,{'X','ckf'});
            legend('show');
    end
end

        