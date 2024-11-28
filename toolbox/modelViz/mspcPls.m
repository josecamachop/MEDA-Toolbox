
function [Dst,Qst,Dstt,Qstt,UCLd,UCLq] = mspcPls(x,y,varargin)
% Compute D-st and Q-st in PLS-based Multivariate Statistical Process 
% Control
%
% [Dst,Qst] = mspcPls(x,y) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
%
% Optional INPUTS (parameters):
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 1:rank(x)
%
% 'ObsTest': [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% 'PreprocessingX': [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
%
% 'PreprocessingY': [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
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
% 'ObsLabel': [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% 'ObsClass': [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
%
% 'PValueD': [Ldx1] p-values for control limits in the D-st, in (0,1]. 
%   Values equal to 0.01 and 0.05 are used by default in bar plots, and
%   0.01 in scatter plots
%
% 'PValueQ': [Lqx1] p-values for control limits in the Q-st, in (0,1]. 
%   Values equal to 0.01 and 0.05 are used by default in bar plots, and
%   0.01 in scatter plots
%
% 'LimType': [1x1] type of control limit
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
% Qst: [Nx1] Q-statistic of test
%
% UCLd: [Ldx1] Control limits in D-statistic
%
% UCLq: [Lqx1] Control limits in Q-statistic
%
%
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(100,10);
% Y = 0.1*randn(100,2) + X(:,1:2);
% [Dst,Qst] = mspcPls(X,Y,'LVs',1:1);
%
%
% EXAMPLE OF USE: PLS-based MSPC on NOC test data and anomalies.
%
% nobs = 100;
% nvars = 10;
% nPCs = 1;
% X = simuleMV(nobs,nvars,'LevelCorr',6);
% Y = 0.1*randn(100,2) + X(:,1:2);
% 
% lvs = 1:nPCs;
% 
% nobst = 10;
% test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',cov(X)*(nobst-1));
% test(6:10,:) = 3*test(6:10,:);
% 
% [Dst,Qst,Dstt,Qstt] = mspcPls(X,Y,'LVs',lvs,'ObsTest',test);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 20/Nov/2024
%
% Copyright (C) 2024  University of Granada, Granada
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'LVs',1:rank(x)); 
addParameter(p,'ObsTest',[]);
addParameter(p,'PreprocessingX',2);
addParameter(p,'PreprocessingY',2);
addParameter(p,'Option','100');  
addParameter(p,'ObsLabel',[]);  
addParameter(p,'ObsClass',[]);  
addParameter(p,'PValueD',0.01);  
addParameter(p,'PValueQ',0.01);  
addParameter(p,'LimType',0);  
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;
test = p.Results.ObsTest;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
opt = p.Results.Option;
label = p.Results.ObsLabel;
classes = p.Results.ObsClass;
pvalueD = p.Results.PValueD;
pvalueQ = p.Results.PValueQ;
limtype = p.Results.LimType;

L = size(test, 1);
% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'00'); end
if length(opt)<3, opt = strcat(opt,'0'); end
if opt(3) == 1 || opt(3) == '1'
    K = L;
else
    K = N+L;
end

if  isempty(label) 
    if opt(3) == 1 || opt(3) == '1'
        label = 1:L;
    else
        label = [1:N 1:L]; 
    end
end
if isempty(classes)
    if opt(3) == 1 || opt(3) == '1' 
        classes = ones(L,1); 
    else
        classes = [ones(N,1);2*ones(L,1)];  
    end
end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;
if ~isempty(pvalueD) && size(pvalueD,1) == 1,     pvalueD = pvalueD'; end;
if ~isempty(pvalueQ) && size(pvalueQ,1) == 1, pvalueQ = pvalueQ'; end;
Ld = size(pvalueD,1);
Lq = size(pvalueQ,1);

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''LVs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: parameter ''ObsTest'' must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==3, 'Dimension Error: parameter ''Option'' must be a string or num of 3 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: parameter ''ObsLabel'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: parameter ''ObsClass'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(pvalueD), assert (isequal(size(pvalueD), [Ld 1]), 'Dimension Error: parameter ''PValueD'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(pvalueQ), assert (isequal(size(pvalueQ), [Lq 1]), 'Dimension Error: parameter ''PValueQ'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
assert (isequal(size(limtype), [1 1]), 'Dimension Error: parameter ''LimType'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(pvalueD), assert (isempty(find(pvalueD<0 | pvalueD>1)), 'Value Error: parameter ''PValueD'' must contain values in (0,1]. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(pvalueQ), assert (isempty(find(pvalueQ<0 | pvalueQ>1)), 'Value Error: parameter ''PValueQ'' must contain values  in (0,1]. Type ''help %s'' for more info.', routine(1).name); end;
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);

%% Main code

[xcs,m,sc] = preprocess2D(x,'Preprocessing',prepx);
ycs = preprocess2D(y,'Preprocessing',prepy);

model = simpls(xcs,ycs,'LVs',lvs);
R = model.altweights;
P = model.loads;
T = xcs*R;

[Dst,Qst] = mspc(xcs,'InvCovarT',inv(cov(T)),'InSubspace',R,'OutSubspace',P);

if ~isempty(test)
    testcs = preprocess2Dapp(test,m,'Scale',sc);
    [Dstt,Qstt] = mspc(testcs,'InvCovarT',inv(cov(T)),'InSubspace',R,'OutSubspace',P);
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
    
    E = xcs - T*P'; 
    UCLq = [];   
    for i=1:Lq
        UCLq(i) = speLim(E,pvalueQ(i));
    end
else
    UCLd = [];   
    for i=1:Ld
        UCLd(i) = prctile(Dst,100*(1-pvalueD(i)));
    end
    
    UCLq = [];   
    for i=1:Lq
        UCLq(i) = prctile(Qst,100*(1-pvalueQ(i)));
    end
end

%% Show results

if opt(1) == '1'
    
    if opt(3) == '0'
        Dsttt = [Dst;Dstt];
        Qsttt = [Qst;Qstt];
    else
        Dsttt = Dstt;
        Qsttt = Qstt;
    end
    
    if opt(2) == '0'
        plotScatter([Dsttt,Qsttt], 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'D-st','Q-st'}, 'LimCont',{UCLd,UCLq});
    else
        plotVec(Dsttt, 'EleLabel', label, 'ObsClass', classes, 'XYLabel', {[],'D-st'}, 'LimCont', UCLd);
        plotVec(Qsttt, 'EleLabel', label, 'ObsClass', classes, 'XYLabel', {[],'Q-st'}, 'LimCont', UCLq);
    end
end
        