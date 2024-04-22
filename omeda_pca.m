
function [omeda_vec,lim] = omeda_pca(x,pcs,test,dummy,varargin)

% Observation-based Missing data methods for Exploratory Data Analysis 
% (oMEDA) for PCA. The original paper is Journal of Chemometrics, 2011, 25 
% (11): 592-600. This algorithm follows the direct computation for
% Known Data Regression (KDR) missing data imputation.
%
% omeda_vec = omeda_pca(x,pcs,test,dummy) % minimum call
% [omeda_vec,lim] = omeda_pca(x,pcs,test,dummy,'Preprocessing',prep,'Option',opt,'VarsLabel',label,'VarsClass',classes) %complete call
%
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
% dummy: [Lx1] dummy variable containing weights for the observations to 
%   compare, and 0 for the rest of observations. By default all test
%   observations are set to 1.
%
% Optional INPUTS (parameters):
%
% 'Preprocessing': [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
%
% 'Option': (str or num) options for data plotting: binary code of the form 'abc' for:
%       a:
%           0: no plots
%           1: plot oMEDA vector
%       b:
%           0: no control limits
%           1: plot control limits 
%       c:
%           0: no normalization
%           1: normalize by control limits
%   By deafult, opt = '100'. If less than 3 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0 and c=0. 
%   If a=0, then b and c are ignored.
%
% 'VarsLabel': [Mx1] name of the variables (numbers are used by default)
%
% 'VarsClass': [Mx1] groups of variables (one group by default)
%
%
% OUTPUTS:
%
% omeda_vec: [Mx1] oMEDA vector.
%
% lim: [Mx1] oMEDA limits.
%
%
% EXAMPLE OF USE: Anomaly on first observation and first 2 variables.
% 
% n_obs = 100;
% n_vars = 10;
% n_PCs = 10;
% X = simuleMV(n_obs,n_vars,'LevelCorr',6);
% 
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,'LevelCorr',6,'Covar',cov(X)*(n_obst-1));
% test(1,1:2) = 10*max(abs(X(:,1:2))); 
% dummy = zeros(10,1);
% dummy(1) = 1;
% 
% pcs = 1:n_PCs;
% 
% omeda_vec = omeda_pca(X,pcs,test,dummy,'Preprocessing',1);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 18/Apr/2024
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
if isempty(pcs), pcs = 1:rank(x); end;
if isempty(test), test = x; end;
L = size(test, 1);
if nargin < 4 || isempty(dummy), dummy = ones(L,1); end;


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Preprocessing',2);     
addParameter(p,'Option','100');  
addParameter(p,'VarsLabel',1:M);  
addParameter(p,'VarsClass',ones(M,1));  
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
prep = p.Results.Preprocessing;
opt = p.Results.Option;
label = p.Results.VarsLabel;
classes = p.Results.VarsClass;

% Convert row arrays to column arrays
if size(label,1) == 1, label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
A = length(pcs);

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'00'); end
if length(opt)<3, opt = strcat(opt,'0'); end

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''pcs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''pcs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(test), [L M]), 'Dimension Error: parameter ''test'' must be L-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(dummy), [L 1]), 'Dimension Error: parameter ''dummy''must be L-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing''must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==3, 'Dimension Error: parameter ''Option'' must be a string or num of 3 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: parameter ''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(classes), [M 1]), 'Dimension Error: parameter ''VarsClass'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 

% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''pcs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[xcs,m,sd] = preprocess2D(x,'Preprocessing',prep);

P = pca_pp(xcs,'Pcs',pcs);
    
testcs = preprocess2Dapp(test,m,'SDivideTest',sd);
omeda_vec = omeda(testcs,dummy,P);

% heuristic: 95% limit for one-observation-dummy
xr = xcs*P*P';
omeda_x = abs((2*xcs-xr).*(xr));
lim = prctile(omeda_x,95)';
    

%% Show results

if opt(1) == '1'
    
    vec = omeda_vec;
 
    if opt(2) == '1'
        limp = lim;
    else
        limp = [];
    end
    
    if opt(3) == '1'
        ind = find(lim>1e-10);
        vec(ind) = vec(ind)./lim(ind);
    	if ~isempty(limp)
            limp(ind) = limp(ind)./lim(ind);
        end
    end
    
    plot_vec(vec,'EleLabel',label,'ObsClass',classes,'XYLabel',{[],'d^2_A'},'LimCont',[limp -limp]);
    
end

        