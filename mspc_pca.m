
function [Dst,Qst,Dstt,Qstt] = mspc_pca(x,pcs,test,prep,opt,label,classes,p_valueD,p_valueQ)

% Compute D-st and Q-st in PCA-based Multivariate Statistical Process 
% Control
%
% [Dst,Qst] = mspc_pca(x) % minimum call
% [Dst,Qst,Dstt,Qstt] = mspc_pca(x,pcs,test,prep,opt,label,classes,p_valueD,p_valueQ) % complete call
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
%       1: scatter plot (default)
%       otherwise: bar plot of each single statistic
%
% label: [Kx1] K=N+L, name of the observations (numbers are used by default)
%
% classes: [Kx1] K=N+L, groups for different visualization (a single group 
%   by default per calibration and test)
%
% p_valueD: [Ldx1] p-values for control limits in the D-st, in (0,1]. 
%   Values equal to 0.01 and 0.05 are used by default in bar plots, and
%   0.01 in scatter plots.
%
% p_valueQ: [Lqx1] p-values for control limits in the Q-st, in (0,1]. 
%   Values equal to 0.01 and 0.05 are used by default in bar plots, and
%   0.01 in scatter plots.
%
%
% OUTPUTS:
%
% Dst: [Nx1] D-statistic or Hotelling T2 of calibration
%
% Qst: [Nx1] Q-statistic of calibration
%
% Dstt: [Nx1] D-statistic or Hotelling T2 of test
%
% Qst: [Nx1] Q-statistic of test
%
%
% EXAMPLE OF USE: Random scores
%
% X = real(ADICOV(randn(10,10).^19,randn(100,10),10));
% [Dst,Qst] = mspc_pca(X,1:3);
%
%
% EXAMPLE OF USE: PCA-based MSPC on NOC test data and anomalies.
%
% n_obs = 100;
% n_vars = 10;
% n_PCs = 1;
% XX = randn(n_vars,n_vars).^19; 
% X = real(ADICOV(n_obs*XX,randn(n_obs,n_vars),n_vars));
% 
% pcs = 1:n_PCs;
% 
% n_obst = 10;
% test = real(ADICOV(n_obst*XX,randn(n_obst,n_vars),n_vars));
% test(6:10,:) = (1 + 1)*test(6:10,:);
% 
% [Dst,Qst,Dstt,Qstt] = mspc_pca(X,pcs,test,2,1,[],[ones(100,1);2*ones(5,1);3*ones(5,1)]);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 31/Mar/2016
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
if nargin < 8, 
    if opt == 1,
        p_valueD = 0.01; 
    else
        p_valueD = [0.01 0.05]; 
    end
end;
if nargin < 9, 
    if opt == 1,
        p_valueQ = 0.01; 
    else
        p_valueQ = [0.01 0.05]; 
    end
end;

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;
if ~isempty(p_valueD) && size(p_valueD,1) == 1,     p_valueD = p_valueD'; end;
if ~isempty(p_valueQ) && size(p_valueQ,1) == 1, p_valueQ = p_valueQ'; end;
Ld = size(p_valueD,1);
Lq = size(p_valueQ,1);

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
if ~isempty(p_valueD), assert (isequal(size(p_valueD), [Ld 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(p_valueQ), assert (isequal(size(p_valueQ), [Lq 1]), 'Dimension Error: 9th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;

% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: 2nd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(p_valueD), assert (isempty(find(p_valueD<0 || p_valueD>1)), 'Value Error: 8th argument must contain values in (0,1]. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(p_valueQ), assert (isempty(find(p_valueQ<0 || p_valueQ>1)), 'Value Error: 9th argument must contain values  in (0,1]. Type ''help %s'' for more info.', routine(1).name); end;


%% Main code

[xcs,m,sc] = preprocess2D(x,prep);
[P,T] = pca_pp(xcs,pcs);

[Dst,Qst] = mspc(xcs,inv(cov(T)),P);

if ~isempty(test)
    testcs = preprocess2Dapp(test,m,sc);
    [Dstt,Qstt] = mspc(testcs,inv(cov(T)),P);
else
    Dstt = [];
    Qstt = [];
end


%% Show results

if opt,

    UCLd = [];
    for i=1:Ld,
        if isempty(test),
            UCLd(i) = hot_lim(A,N,p_valueD(i),1);
        else
            UCLd(i) = hot_lim(A,N,p_valueD(i),2);
        end
    end
    
    E = xcs - T*P'; 
    UCLq = [];   
    for i=1:Lq,
        UCLq(i) = spe_lim(E,p_valueQ(i));
    end
    
    Dsttt = [Dst;Dstt];
    Qsttt = [Qst;Qstt];
    
    if opt ~=1,
        plot_vec(Dsttt, label, classes, {[],'D-st'}, UCLd);
        plot_vec(Qsttt, label, classes, {[],'Q-st'}, UCLq);
    else
        plot_scatter([Dsttt,Qsttt], label, classes, {'D-st','Q-st'}, {UCLd,UCLq});

    end
end
        