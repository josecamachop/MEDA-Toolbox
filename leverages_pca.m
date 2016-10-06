
function L = leverages_pca(x,pcs,prep,opt,label,classes)

% Compute and plot the leverages of variables in the PCA model
%
% L = leverages_pca(x) % minimum call
% L = leverages_pca(x,pcs,prep,opt,label,classes) % complete call
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 1:rank(xcs)
%
% prep: [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
%
% opt: (str or num) options for data plotting
%       0: no plots.
%       1: plot bar plot of leverages (default) 
%
% label: [Mx1] name of the variables (numbers are used by default)
%
% classes: [Mx1] groups for different visualization (a single group 
%   by default)
%
%
% OUTPUTS:
%
% L: [Mx1] leverages of the variables
%
%
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(20,10,8);
% L = leverages_pca(X,1:3);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Apr/2016
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
if nargin < 2 || isempty(pcs), pcs = 1:rank(x); end;
if nargin < 3 || isempty(prep), prep = 2; end;
if nargin < 4 || isempty(opt), opt = 1; end; 
if nargin < 5 || isempty(label), label = [1:M]; end
if nargin < 6 || isempty(classes), classes = ones(M,1); end

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
A = length(pcs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: 2nd argument with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 2nd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: 5th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: 6th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: 2nd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 4th argument must contain a binary value. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,prep);
P = pca_pp(xcs,pcs);

L = diag(P*P');


%% Show results

if opt == '1', 
    plot_vec(L, label, classes, {'Variables','Leverages'});
end
        