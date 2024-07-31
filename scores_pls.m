
function [T,TT] = scores_pls(x,y,varargin)

% Compute and plot scores in PLS. This routine is deprecated and superseded 
% by scores.m (please, use the latter)
%
% T = scores_pls(x,y) % minimum call
% [T,TT] = scores_pls(x,y,'LatVars',lvs,'ObsTest',test,'PreprocessingX',prepx,'PreprocessingY',prepy,'Option',opt,'ObsLabel',label,'ObsClass',classes,'BlurIndex',blur) % complete call
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
% Optional Inputs (parameters):
%
% 'LatVars': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 1:rank(x)
%
% 'ObsTest': [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% 'PreprocessingX': [1x1] preprocesing of the x-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)  
%
% 'PreprocessingY': [1x1] preprocesing of the y-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)   
%
% 'Option': (str or num) options for data plotting: binary code of the form 'abcd' for:
%       a:
%           0: no plots
%           1: plot scores
%       b:
%           0: scatter plot of pairs of LVs 
%           1: bar plot of each single LV
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
% 'ObsLabel': [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% 'ObsClass': [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
%
% 'BlurIndex': [1x1] avoid blur when adding labels. The higher, the more labels 
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
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% T = scores_pls(X,Y,'LatVars',1:3);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% n_obs = 100;
% n_vars = 10;
% n_PCs = 10;
% X = simuleMV(n_obs,n_vars,'LevelCorr',8);
% Y = 0.1*randn(n_obs,2) + X(:,1:2);
% 
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,'LevelCorr',6,'Covar',corr(X)*(n_obst-1)/(n_obs-1));
% 
% scores_pls(X,Y,'LatVars',1,'ObsTest',test);
% scores_pls(X,Y,'LatVars',1:2,'ObsTest',test);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 23/Apr/2024
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
addParameter(p,'LatVars',1:rank(x));   
addParameter(p,'ObsTest',[]);   
addParameter(p,'Option','100');
addParameter(p,'PreprocessingX',2);
addParameter(p,'PreprocessingY',2);
addParameter(p,'ObsLabel',[]);
addParameter(p,'ObsClass',[]);
addParameter(p,'BlurIndex',1);
parse(p,varargin{:});


% Extract inputs from inputParser for code legibility
test = p.Results.ObsTest;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
lvs = p.Results.LatVars;
opt = p.Results.Option;
label = p.Results.ObsLabel;
classes = p.Results.ObsClass;
blur = p.Results.BlurIndex;

L = size(test, 1);
K = N+L;
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

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''LatVars'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LatVars'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: parameter ''ObsTest'' must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==4, 'Dimension Error: parameter ''Option''  must be a string or num of 4 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: parameter ''ObsLabel'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: parameter ''ObsClass'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LatVars'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[xcs,m,sd] = preprocess2D(x,'Preprocessing',prepx);
ycs = preprocess2D(y,'Preprocessing',prepy);

[beta,W,P,Q,R] = simpls(xcs,ycs,'LatVars',lvs); 
T = xcs*R;

if ~isempty(test)
    testcs = preprocess2Dapp(test,m,'SDivideTest',sd);
    TT = testcs*R;
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
    
    if length(lvs) == 1 || opt(2) == '1'
        for i=1:length(lvs)
            plot_vec(Tt(:,i), 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('Scores LV %d (%.0f%%)',lvs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2)))});
        end
    else
        for i=1:length(lvs)-1
            for j=i+1:length(lvs)
                plot_scatter([Tt(:,i),Tt(:,j)], 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{sprintf('Scores LV %d (%.0f%%)',lvs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2))),sprintf('Scores LV %d (%.0f%%)',lvs(j),100*sum(T(:,j).^2)/sum(sum(xcs.^2)))}','Option',opt(4),'BlurIndex',blur);
            end      
        end
    end
end
        