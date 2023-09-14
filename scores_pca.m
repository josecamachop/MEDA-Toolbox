
function [T,TT] = scores_pca(x,pcs,test,prep,opt,label,classes,blur)


% Compute and plot scores in PCA. This routine is deprecated and superseded 
% by scores.m (please, use the latter)
%
% T = scores_pca(x) % minimum call
% [T,TT] = scores_pca(x,pcs,test,prep,opt,label,classes,blur) % complete call
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
% opt: (str or num) options for data plotting: binary code of the form 'abcd' for:
%       a:
%           0: no plots
%           1: plot scores
%       b:
%           0: scatter plot of pairs of PCs 
%           1: bar plot of each single PC
%       c:
%           0: plot calibration and test data
%           1: plot only test data 
%       d:
%           0: plot for categorical classes (consistent with a legend)
%           1: plot for numerical classes (consistent with a colorbar)
%
%   By deafult, opt = '1000'. If less than 4 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0, c=0 and 
%   d=0. If a=0, then b, c  and d are ignored.
%
% label: [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% classes: [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
%
% blur: [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default).
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
% X = simuleMV(20,10,8);
% T = scores_pca(X,1:3);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% n_obs = 100;
% n_vars = 10;
% X = simuleMV(n_obs,n_vars,8);
%
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,6,corr(X)*(n_obst-1)/(n_obs-1));
%
% scores_pca(X,1,test);
% scores_pca(X,1:2,test);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 21/Apr/2023
%
% Copyright (C) 2023  University of Granada, Granada
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
if nargin < 4 || isempty(prep), prep = 2; end;
if nargin < 5 || isempty(opt), opt = '100'; end; 

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
while length(opt)<4, opt = strcat(opt,'0'); end

if opt(4) == '0', opt(4) = '1'; else,  opt(4) = '0'; end
if opt(3) == 1 || opt(3) == '1'
    K = L;
else
    K = N+L;
end

if nargin < 6 || isempty(label) 
    if opt(3) == 1 || opt(3) == '1'
        label = 1:L;
    else
        label = [1:N 1:L]; 
    end
end
if nargin < 7 || isempty(classes)
    if opt(3) == 1 || opt(3) == '1' 
        classes = ones(L,1); 
    else
        classes = [ones(N,1);2*ones(L,1)];  
    end
end
if nargin < 8 || isempty(blur),    blur    = 1;       end;

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
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: 3rd argument must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (isequal(size(prep), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==4, 'Dimension Error: 5th argument must be a string or num of 4 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: 6th argument must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: 7th argument must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: 2nd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 5th argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[xcs,m,sd] = preprocess2D(x,prep);
[P,T] = pca_pp(xcs,pcs);

if ~isempty(test)
    testcs = preprocess2Dapp(test,m,sd);
    TT = testcs*P;
else
    TT = [];
end


%% Show results

if opt(1) == '1'
    
     if opt(3) == '0'
        Tt = [T;TT];
    else
        Tt = TT;
    end
    
    if length(pcs) == 1 || opt(2) == '1'
        for i=1:length(pcs)
            plot_vec(Tt(:,i), label, classes, {'',sprintf('Scores PC %d (%.0f%%)',pcs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2)))});
        end
    else
        for i=1:length(pcs)-1
            for j=i+1:length(pcs)
                plot_scatter([Tt(:,i),Tt(:,j)], label, classes, {sprintf('Scores PC %d (%.0f%%)',pcs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2))),sprintf('Scores PC %d (%.0f%%)',pcs(j),100*sum(T(:,j).^2)/sum(sum(xcs.^2)))}',[],opt(4),[],[],blur);
            end      
        end
    end
end
        