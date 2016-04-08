
function [T,TT] = scores_pca(x,pcs,test,prep,opt,label,classes)

% Compute and plot scores in PCA
%
% T = scores_pca(x) % minimum call
% [T,TT] = scores_pca(x,pcs,test,prep,opt,label,classes) % complete call
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 1:rank(xcs)
%
% test: [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% prep: [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
%
% opt: [1x1] options for data plotting
%       0: no plots
%       1: scatter plot of pairs of PCs (default)
%       otherwise: bar plot of each single PC
%
% label: [Kx1] K=N+L, name of the observations (numbers are used by default)
%
% classes: [Kx1] K=N+L, groups for different visualization (a single group 
%   by default per calibration and test)
%
%
% OUTPUTS:
%
% T: [NxA] calibration scores
%
% TT: [NxA] test scores
%
%
% EXAMPLE OF USE: Random scores
%
% X = real(ADICOV(randn(10,10).^19,randn(100,10),10));
% T = scores_pca(X,1:3);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% n_obs = 100;
% n_vars = 10;
% XX = randn(n_vars,n_vars).^19; 
% X = real(ADICOV(n_obs*XX,randn(n_obs,n_vars),n_vars));
%
% n_obst = 10;
% test = real(ADICOV(n_obst*XX,randn(n_obst,n_vars),n_vars));
%
% scores_pca(X,1,test);
% scores_pca(X,1:2,test);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 05/Apr/2016
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
if nargin < 2 || isempty(pcs), pcs = 1:rank(x); end;
if nargin < 3, test = []; end;
L = size(test, 1);
K = N+L;
if nargin < 4 || isempty(prep), prep = 2; end;
if nargin < 5 || isempty(opt), opt = 1; end; 
if nargin < 6 || isempty(label), label = [1:N 1:L]; end
if nargin < 7 || isempty(classes), classes = [ones(N,1);2*ones(L,1)]; end

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
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 2nd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: 3rd argument must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (isequal(size(prep), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: 6th argument must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: 7th argument must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: 2nd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[xcs,m,sd] = preprocess2D(x,prep);
[P,T] = pca_pp(xcs,pcs);

if ~isempty(test),
    testcs = preprocess2Dapp(test,m,sd);
    TT = testcs*P;
else
    TT = [];
end


%% Show results

if opt,
    Tt = [T;TT];
    if length(pcs) == 1 || opt ~=1,
        for i=1:length(pcs),
            plot_vec(Tt(:,i), label, classes, {'',sprintf('Scores PC %d',pcs(i))});
        end
    else
        for i=1:length(pcs)-1,
            for j=i+1:length(pcs),
                plot_scatter([Tt(:,i),Tt(:,j)], label, classes, {sprintf('Scores PC %d',pcs(i)),sprintf('Scores PC %d',pcs(j))}');
            end      
        end
    end
end
        