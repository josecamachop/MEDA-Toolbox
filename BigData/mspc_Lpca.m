
function [Dst,Qst,Dstt,Qstt,UCLd,UCLq] = mspc_Lpca(Lmodel,test,opt,label,classes,p_valueD,p_valueQ,limtype)

% Compute D-st and Q-st in PCA-based Multivariate Statistical Process 
% Control
%
% [Dst,Qst] = mspc_Lpca(Lmodel) % minimum call
% [Dst,Qst,Dstt,Qstt,UCLd,UCLq] = mspc_Lpca(Lmodel,test,opt,label,classes,p_valueD,p_valueQ,limtype) % complete call
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: [MxM] X-block cross-product matrix.
%       Lmodel.lvs: [1x1] number of PCs. 
%       Lmodel.centr: [NxM] centroids of the clusters of observations.
%       Lmodel.multr: [ncx1] multiplicity of each cluster.
%       Lmodel.class: [ncx1] class associated to each cluster.
%       Lmodel.av: [1xM] sample average according to the preprocessing method.
%       Lmodel.sc: [1xM] sample scale according to the preprocessing method.
%       Lmodel.weight: [1xM] weight applied after the preprocessing method.
%       Lmodel.obs_l: {ncx1} label of each cluster.
%
% test: [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% opt: (str or num) options for data plotting: binary code of the form 'abc' for:
%       a:
%           0: no plots
%           1: plot MSPC charts
%       b:
%           0: scatter plot
%           1: bar plot of each single statistic
%       c:
%           0: plot calibration and test data
%           1: plot only test data 
%   By deafult, opt = '100'. If less than 3 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0 and c=0. 
%   If a=0, then b and c are ignored.
%
% label: [Lx1] name of the test observations (numbers are used by default)
%
% classes: [Lx1] groups in test for different visualization (a single group 
%   by default)
%
% p_valueD: [Ldx1] p-values for control limits in the D-st, in (0,1]. 
%   Values equal to 0.01 and 0.05 are used by default in bar plots, and
%   0.01 in scatter plots
%
% p_valueQ: [Lqx1] p-values for control limits in the Q-st, in (0,1]. 
%   Values equal to 0.01 and 0.05 are used by default in bar plots, and
%   0.01 in scatter plots
%
% limtype: [1x1] type of control limit
%       0: theoretical (statistical distribution based, by default)
%       otherwise: percentiles
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
% Qstt: [Nx1] Q-statistic of test
%
% UCLd: [Ldx1] Control limits in D-statistic
%
% UCLq: [Lqx1] Control limits in Q-statistic
%
%
% EXAMPLE OF USE: PCA-based MSPC on NOC test data and anomalies.
%
% n_obs = 100;
% n_vars = 10;
% n_PCs = 1;
% Lmodel = Lmodel_ini(simuleMV(n_obs,n_vars,6));
% Lmodel.multr = 100*rand(n_obs,1); 
% Lmodel.lvs = 1:n_PCs;
% 
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,6,corr(Lmodel.centr)*(n_obst-1)/(Lmodel.N-1));
% test(6:10,:) = 3*test(6:10,:);
% 
% [Dst,Qst,Dstt,Qstt] = mspc_Lpca(Lmodel,test,100,[],[1*ones(5,1);2*ones(5,1)]);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 27/Aug/18.
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez
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

[ok, Lmodel] = check_Lmodel(Lmodel);

N = min(Lmodel.nc,size(Lmodel.centr,1));
M = size(Lmodel.XX, 2);

if nargin < 2, test = []; end;
L = size(test, 1);

if nargin < 3 || isempty(opt), opt = '100'; end; 

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'00'); end
if length(opt)<3, opt = strcat(opt,'0'); end
if opt(3) == '1'
    K = L;
else
    K = N+L;
end

if nargin < 4 || isempty(label)
    if  opt(3) == '1'
        label = cellstr(num2str((1:L)'));
    elseif isempty(Lmodel.obs_l)
        label = cellstr(num2str([1:N 1:L]'));
    else
        if L
            lb1 = cellstr(num2str((1:L)'));
            label = {Lmodel.obs_l{:} lb1{:}};
        else
            label = Lmodel.obs_l;
        end
    end
else
    if  opt(3) == '0'
        if isempty(Lmodel.obs_l)
            lb1 = cellstr(num2str((1:N)'));
            label = {lb1{:} label{:}};
        else
            label = {Lmodel.obs_l{:} label{:}};
        end
    end
end

if nargin < 5 || isempty(classes)
    if opt(3) == '1' 
        classes = ones(L,1); 
    else
        classes = [Lmodel.class;2*ones(L,1)];  
    end
elseif opt(3) == '0' && length(classes)==L
        classes = [Lmodel.class;2*classes];
end

if nargin < 6 || isempty(p_valueD)
    if opt(2) == 0 || opt(2) == '0'
        p_valueD = 0.01; 
    else
        p_valueD = [0.01 0.05]; 
    end
end;
if nargin < 7 || isempty(p_valueQ) 
    if opt(2) == 0 || opt(2) == '0'
        p_valueQ = 0.01; 
    else
        p_valueQ = [0.01 0.05]; 
    end
end;
if nargin < 8, limtype = 0; end;

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;
if ~isempty(p_valueD) && size(p_valueD,1) == 1,     p_valueD = p_valueD'; end;
if ~isempty(p_valueQ) && size(p_valueQ,1) == 1, p_valueQ = p_valueQ'; end;
Ld = size(p_valueD,1);
Lq = size(p_valueQ,1);

A = length(Lmodel.lvs);

% Validate dimensions of input data
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: 2nd argument must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (ischar(opt) && length(opt)==3, 'Dimension Error: 3rd argument must be a string or num of maximum 3 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: 4th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: 5th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(p_valueD), assert (isequal(size(p_valueD), [Ld 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(p_valueQ), assert (isequal(size(p_valueQ), [Lq 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
assert (isequal(size(limtype), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3rd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(p_valueD), assert (isempty(find(p_valueD<=0 | p_valueD>1)), 'Value Error: 6th argument must contain values in (0,1]. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(p_valueQ), assert (isempty(find(p_valueQ<=0 | p_valueQ>1)), 'Value Error: 7th argument must contain values  in (0,1]. Type ''help %s'' for more info.', routine(1).name); end;


%% Main code

[P,~,Lmodel] = Lpca(Lmodel);
iTT = diag((Lmodel.N-1)./Lmodel.LVvar(Lmodel.lvs));

[Dst,Qst] = mspc(Lmodel.centr,iTT,P);

if ~isempty(test)
    testcs = preprocess2Dapp(test,Lmodel.av,Lmodel.sc,Lmodel.weight);
    [Dstt,Qstt] = mspc(testcs,iTT,P);
else
    Dstt = [];
    Qstt = [];
end

if limtype==0
    UCLd = [];
    for i=1:Ld
        if isempty(test)
            UCLd(i) = hot_lim(A,N,p_valueD(i),1);
        else
            UCLd(i) = hot_lim(A,N,p_valueD(i),2);
        end
    end
    
    E = Lmodel.centr - Lmodel.centr*P*P'; 
    UCLq = [];   
    for i=1:Lq
        UCLq(i) = spe_lim(E,p_valueQ(i));
    end
else
    UCLd = [];   
    for i=1:Ld
        UCLd(i) = prctile(repelem(Dst,Lmodel.multr),100*(1-p_valueD(i)));
    end
    
    UCLq = [];   
    for i=1:Lq
        UCLq(i) = prctile(repelem(Qst,Lmodel.multr),100*(1-p_valueQ(i)));
    end
end

%% Show results

if opt(1) == '1'
     
    if opt(3) == '0'
        Dsttt = [Dst;Dstt];
        Qsttt = [Qst;Qstt];
        mult = [Lmodel.multr;ones(size(Dstt))];
    else
        Dsttt = Dstt;
        Qsttt = Qstt;
        mult = ones(size(Dstt));
    end
    
    if opt(2) == '0'
        plot_scatter([Dsttt,Qsttt], label, classes, {'D-st','Q-st'}, {UCLd,UCLq}, 11, mult);
    else
        plot_vec(Dsttt, label, classes, {[],'D-st'}, UCLd, 0, mult);
        plot_vec(Qsttt, label, classes, {[],'Q-st'}, UCLq, 0, mult);
    end
end
        