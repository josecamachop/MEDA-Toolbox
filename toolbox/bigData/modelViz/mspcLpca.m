
function [Dst,Qst,Dstt,Qstt,UCLd,UCLq] = mspcLpca(Lmodel,varargin)

% Compute D-st and Q-st in PCA-based Multivariate Statistical Process 
% Control
%
% [Dst,Qst] = mspcLpca(Lmodel) % minimum call
% [Dst,Qst,Dstt,Qstt,UCLd,UCLq] = mspcLpca(Lmodel,test,opt,label,classes,pvalueD,pvalueQ,limtype) % complete call
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
%       Lmodel.obsl: {ncx1} label of each cluster.
%
% Optional INPUTS (parameters):
%
% 'Test': [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% 'Option': (str or num) options for data plotting: binary code of the form 'abc' for:
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
% 'ObsLabel': [Lx1] name of the test observations (numbers are used by default)
%
% 'ObsClass': [Lx1] groups in test for different visualization (a single group 
%   by default)
%
% 'pvalueD': [Ldx1] p-values for control limits in the D-st, in (0,1]. 
%   Values equal to 0.01 and 0.05 are used by default in bar plots, and
%   0.01 in scatter plots
%
% 'pvalueQ': [Lqx1] p-values for control limits in the Q-st, in (0,1]. 
%   Values equal to 0.01 and 0.05 are used by default in bar plots, and
%   0.01 in scatter plots
%
% 'limType': [1x1] type of control limit
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
% nobs = 100;
% nvars = 10;
% nPCs = 1;
% Lmodel = iniLmodel(simuleMV(nobs,nvars,'LevelCorr',6));
% Lmodel.multr = 100*rand(nobs,1); 
% Lmodel.lvs = 1:nPCs;
% 
% nobst = 10;
% test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(Lmodel.centr)*(nobst-1)/(Lmodel.N-1));
% test(6:10,:) = 3*test(6:10,:);
% 
% [Dst,Qst,Dstt,Qstt] = mspcLpca(Lmodel,'Test',test,'Option',100,'ObsLabel',[],'ObsClass',[1*ones(5,1);2*ones(5,1)]);
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 03/Feb/2026
%
% Copyright (C) 2026  University of Granada, Granada
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

[ok, Lmodel] = checkLmodel(Lmodel);

N = min(Lmodel.nc,size(Lmodel.centr,1));
M = size(Lmodel.XX, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Test',[]);
addParameter(p,'Option','100');
addParameter(p,'ObsLabel',[]);
addParameter(p,'ObsClasses',[]);
addParameter(p,'pvalueD',[]);
addParameter(p,'pvalueQ',[]);
addParameter(p,'limType',0);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
test = p.Results.Test;
opt = p.Results.Option;
label = p.Results.ObsLabel;
classes = p.Results.ObsClasses;
pvalueQ = p.Results.pvalueQ;
pvalueD = p.Results.pvalueD;
limtype = p.Results.limType;

L = size(test, 1);

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

if isempty(label)
    if  opt(3) == '1'
        label = cellstr(num2str((1:L)'));
    elseif isempty(Lmodel.obsl)
        label = cellstr(num2str([1:N 1:L]'));
    else
        if L
            lb1 = cellstr(num2str((1:L)'));
            label = {Lmodel.obsl{:} lb1{:}};
        else
            label = Lmodel.obsl;
        end
    end
else
    if  opt(3) == '0'
        if isempty(Lmodel.obsl)
            lb1 = cellstr(num2str((1:N)'));
            label = {lb1{:} label{:}};
        else
            label = {Lmodel.obsl{:} label{:}};
        end
    end
end

if isempty(classes)
    if opt(3) == '1' 
        classes = ones(L,1); 
    else
        classes = [Lmodel.class;2*ones(L,1)];  
    end
elseif opt(3) == '0' && length(classes)==L
        classes = [Lmodel.class;2*classes];
end

if isempty(pvalueD)
    if opt(2) == 0 || opt(2) == '0'
        pvalueD = 0.01; 
    else
        pvalueD = [0.01 0.05]; 
    end
end;
if isempty(pvalueQ) 
    if opt(2) == 0 || opt(2) == '0'
        pvalueQ = 0.01; 
    else
        pvalueQ = [0.01 0.05]; 
    end
end;

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;
if ~isempty(pvalueD) && size(pvalueD,1) == 1, pvalueD = pvalueD'; end;
if ~isempty(pvalueQ) && size(pvalueQ,1) == 1, pvalueQ = pvalueQ'; end;
Ld = size(pvalueD,1);
Lq = size(pvalueQ,1);

A = length(Lmodel.lvs);

% Validate dimensions of input data
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: 2nd argument must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (ischar(opt) && length(opt)==3, 'Dimension Error: 3rd argument must be a string or num of maximum 3 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: 4th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: 5th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(pvalueD), assert (isequal(size(pvalueD), [Ld 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(pvalueQ), assert (isequal(size(pvalueQ), [Lq 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
assert (isequal(size(limtype), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3rd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(pvalueD), assert (isempty(find(pvalueD<=0 | pvalueD>1)), 'Value Error: 6th argument must contain values in (0,1]. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(pvalueQ), assert (isempty(find(pvalueQ<=0 | pvalueQ>1)), 'Value Error: 7th argument must contain values  in (0,1]. Type ''help %s'' for more info.', routine(1).name); end;


%% Main code

Lmodel = Lpca(Lmodel);
P = Lmodel.loads;
iTT = diag((Lmodel.N-1)./Lmodel.sdT(Lmodel.lvs));

[Dst,Qst] = mspc(Lmodel.centr,'InvCovarT',iTT,'InSubspace',P);

if ~isempty(test)
    testcs = preprocess2Dapp(test,Lmodel.av,'Scale',Lmodel.sc,'Weight',Lmodel.weight);
    [Dstt,Qstt] = mspc(testcs,'InvCovarT',iTT,'InSubspace',P);
else
    Dstt = [];
    Qstt = [];
end

if limtype==0
    UCLd = [];
    for i=1:Ld
        if isempty(test)
            UCLd(i) = hotLim(A,N,pvalueD(i),'Phase',1);
        else
            UCLd(i) = hotLim(A,N,pvalueD(i),'Phase',2);
        end
    end
    
    E = Lmodel.centr - Lmodel.centr*P*P'; 
    UCLq = [];   
    for i=1:Lq
        UCLq(i) = speLim(E,pvalueQ(i));
    end
else
    UCLd = [];   
    for i=1:Ld
        UCLd(i) = prctile(repelem(Dst,Lmodel.multr),100*(1-pvalueD(i)));
    end
    
    UCLq = [];   
    for i=1:Lq
        UCLq(i) = prctile(repelem(Qst,Lmodel.multr),100*(1-pvalueQ(i)));
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
        plotScatter([Dsttt,Qsttt], 'EleLabel', label, 'ObsClass', classes, 'XYLabel', {'D-st','Q-st'}, 'LimCont', {UCLd,UCLq}, 'Multiplicity', mult);
    else
        plotVec(Dsttt, 'EleLabel', label, 'ObsClass', classes, 'XYLabel', {[],'D-st'}, 'LimCont', UCLd, 'Multiplicity', mult);
        plotVec(Qsttt, 'EleLabel', label, 'ObsClass', classes, 'XYLabel', {[],'Q-st'}, 'LimCont', UCLq, 'Multiplicity', mult);
    end
end
        