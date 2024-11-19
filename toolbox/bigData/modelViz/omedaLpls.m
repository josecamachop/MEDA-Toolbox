function omedavec = omedaLpls(Lmodel,test,dummy,opt)

% Observation-based Missing data methods for Exploratory Data Analysis 
% (oMEDA) for PCA. The original paper is Journal of Chemometrics, 
% DOI: 10.1002/cem.1405. This algorithm follows the direct computation for
% Known Data Regression (KDR) missing data imputation.
%
% omedavec = omedaLpls(Lmodel,test,dummy) % minimum call
% [omedavec,lim] = omedaLpls(Lmodel,test,dummy,opt) %complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: [MxM] X-block cross-product matrix.
%       Lmodel.XY: [MxO] cross-product matrix between the x-block and the
%           y-block.
%       Lmodel.lvs: [1x1] number of PCs. 
%       Lmodel.centr: [NxM] centroids of the clusters of observations.
%       Lmodel.av: [1xM] sample average according to the preprocessing method.
%       Lmodel.sc: [1xM] sample scale according to the preprocessing method.
%       Lmodel.weight: [1xM] weight applied after the preprocessing method.
%       Lmodel.varl: {ncx1} label of each variable.
%
% test: [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% dummy: [Lx1] dummy variable containing weights for the observations to 
%   compare, and 0 for the rest of observations
%
% opt: (str or num) options for data plotting: binary code of the form 'abc' for:
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
%
% OUTPUTS:
%
% omedavec: [Mx1] oMEDA vector.
%
% lim: [Mx1] oMEDA limits.
%
%
% EXAMPLE OF USE: Anomaly on first observation and first 2 variables.
%
% nobs = 100;
% nvars = 10;
% nLVs = 10;
% X = simuleMV(nobs,nvars,'LevelCorr',6);
% Y = 0.1*randn(nobs,2) + X(:,1:2);
% Lmodel = iniLmodel(X,Y);
% Lmodel.multr = 100*rand(nobs,1); 
% Lmodel.lvs = 1:nLVs;
% 
% nobst = 10;
% test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(Lmodel.centr)*(nobst-1)/(Lmodel.N-1));
% test(1,1:2) = 10*max(abs(Lmodel.centr(:,1:2))); 
% dummy = zeros(10,1);
% dummy(1) = 1;
%
% omedavec = omedaLpls(Lmodel,test,dummy);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 19/Nov/2024
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

checkLmodel(Lmodel);

N = Lmodel.nc;
M = size(Lmodel.XX, 2);
L = size(test, 1);
if isempty(dummy), dummy = ones(L,1); end;
if nargin < 4 || isempty(opt), opt = '100'; end; 

A = length(Lmodel.lvs);

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'00'); end
if length(opt)<3, opt = strcat(opt,'0'); end

% Validate dimensions of input data
assert (A>0, 'Dimension Error: 1sr argument with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(test), [L M]), 'Dimension Error: 2nd argument must be L-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(dummy), [L 1]), 'Dimension Error: 3rd argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==3, 'Dimension Error: 4th argument must be a string or num of maximum 3 bits. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 4th argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

Lmodel = Lpls(Lmodel);
R = Lmodel.altweights;
P = Lmodel.loads;
    
testcs = preprocess2Dapp(test,Lmodel.av,'Scale',Lmodel.sc,'Weight',Lmodel.weight);
omedavec = omeda(testcs,dummy,R,P);

% heuristic: 95% limit for one-observation-dummy
xcs = Lmodel.centr;
xr = xcs*P*P';
omedax = abs((2*xcs-xr).*(xr));
lim = prctile(omedax,95)';
    

%% Show results

if opt(1) == '1'
    
    vec = omedavec;
 
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
    
    plotVec(vec,Lmodel.varl,[],{[],'d^2_A'},[limp -limp]);
    
end



        