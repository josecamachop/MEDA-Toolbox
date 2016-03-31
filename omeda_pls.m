
function [omeda_vec,lim] = omeda_pls(x,y,lvs,test,dummy,prepx,prepy,opt,label)

% Observation-based Missing data methods for Exploratory Data Analysis 
% (oMEDA) for PLS. The original paper is Journal of Chemometrics, 2011, 25 
% (11): 592-600. This algorithm follows the direct computation for
% Known Data Regression (KDR) missing data imputation.
%
% omeda_vec = omeda_pls(x,y,lvs,test,dummy) % minimum call
% [omeda_vec,lim] = omeda_pls(x,y,lvs,test,dummy,prepx,prepy,opt,label) %complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
% lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 1:rank(x)
%
% test: [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% dummy: [Lx1] dummy variable containing weights for the observations to 
%   compare, and 0 for the rest of observations
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
% opt: (str) binary code of the form 'cba' for:
%       a:
%           0: no plots
%           1: plot oMEDA vector (default)
%       b:
%           0: no control limits (default)
%           1: plot control limits 
%       c:
%           0: no normalization (default)
%           1: normalize by control limits)
%
%
% label: [Mx1] name of the variables (numbers are used by default)
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
% n_LVs = 10;
% XX = randn(n_vars,n_vars).^19; 
% X = real(ADICOV(n_obs*XX,randn(n_obs,n_vars),n_vars));
% Y = randn(n_obs,2) + X(:,1:2);
%
% n_obst = 10;
% test = real(ADICOV(n_obst*XX,randn(n_obst,n_vars),n_vars));
% test(1,1:2) = 10*max(abs(X(:,1:2))); 
% dummy = zeros(10,1);
% dummy(1) = 1;
%
% lvs = 1:n_LVs;
%
% omeda_vec = omeda_pls(X,Y,lvs,test,dummy);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
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
assert (nargin >= 5, 'Error in the number of arguments. Type ''help %s'' for more info.', routine.name);
N = size(x, 1);
M = size(x, 2);
O = size(y, 2);
if isempty(lvs), lvs = 1:rank(x); end;
if isempty(test), test = x; end;
L = size(test, 1);
if isempty(dummy), dummyones(L,1); end;
if nargin < 6 || isempty(prepx), prepx = 2; end;
if nargin < 7 || isempty(prepy), prepy = 2; end;
if nargin < 8 || isempty(opt), opt = 1; end; 
if nargin < 9 || isempty(label), label = 1:M; end

% Convert row arrays to column arrays
if size(label,1) == 1, label = label'; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
A = length(lvs);
if isstruct(opt) % opt backward compatibility
    opt = opt.plot + 10*opt.seriated + 100*opt.discard;
end

% Convert int arrays to str
if isnumeric(opt), opt=fliplr(num2str(opt,'%.3d')); end

% Validate dimensions of input data
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(test), [L M]), 'Dimension Error: 4th argument must be L-by-M. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(dummy), [L 1]), 'Dimension Error: 5th argument must be L-by-1. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (ischar(opt), 'Dimension Error: 8th argument must be a string. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(label), [M 1]), 'Dimension Error: 9th argument must be M-by-1. Type ''help %s'' for more info.', routine.name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: 3rd argument must contain positive integers. Type ''help %s'' for more info.', routine.name);
assert (isempty(find(lvs>rank(x))), 'Value Error: 3rd argument must contain values below the rank of the data. Type ''help %s'' for more info.', routine.name);


%% Main code

[xcs,m,sd] = preprocess2D(x,prepx);
ycs = preprocess2D(y,prepy);

[beta,W,P,Q,R] = kernel_pls(xcs'*xcs,xcs'*ycs,lvs);
        
testcs = preprocess2Dapp(test,m,sd);
omeda_vec = omeda(testcs,dummy,R,P);
    
% heuristic: 95% limit for one-observation-dummy
xr = xcs*R*P';
omeda_x = abs((2*xcs-xr).*(xr));
lim = prctile(omeda_x,95)';


%% Show results

if opt(1) == '1',
    
    vec = omeda_vec;
 
    if opt(2) == '1',
        limp = lim;
    else
        limp = [];
    end
    
    if opt(3) == '1',
        vec = vec./lim;
    	if ~isempty(limp),
            limp = limp./lim;
        end
    end
    
    plot_vec(vec,label,[],{[],'d^2_A'},[limp -limp]);
    
end

        