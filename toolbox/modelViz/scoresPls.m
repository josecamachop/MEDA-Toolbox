
function [T,TT] = scoresPls(x,y,varargin)

% Compute and plot scores in PLS. This routine is deprecated and superseded 
% by scores.m (please, use the latter)
%
% T = scoresPls(x,y) % minimum call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
%
% Optional Inputs (parameters):
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
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
% 'PlotType': str
%      'Scatter': scatterplot (by default for pairs of PCs)
%      'Bars': bar plot (by default for a single PC)
%
% 'PlotCal': bool
%      false: plot only test data
%      true: plot both calibration and test (by default)
%
% 'ObsLabel': [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% 'ObsClass': [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
%
% 'BlurIndex': [1x1] to avoid blur when adding labels. It reflects the
%   minimum distance with other points where a label is allowed to be 
%   visualized. For a value of 0, all labels are printed, while for a 
%   large value only uncluttered labels are printed. By default Inf is 
%   chosen, where only indices as visualized. 
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
% T = scoresPls(X,Y,'LVs',1:3);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% nobs = 100;
% nvars = 10;
% nPCs = 10;
% X = simuleMV(nobs,nvars,'LevelCorr',8);
% Y = 0.1*randn(nobs,2) + X(:,1:2);
% 
% nobst = 10;
% test = simuleMV(nobst,nvars,'LevelCorr',8,'Covar',corr(X)*(nobst-1)/(nobs-1));
% 
% scoresPls(X,Y,'LVs',1,'ObsTest',test);
% scoresPls(X,Y,'LVs',1:2,'ObsTest',test);
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 14/Aug/2025
% Dependencies: Matlab R2017b, MEDA v1.9
%
% Copyright (C) 2025  University of Granada, Granada
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
addParameter(p,'PlotType','Scatter');
addParameter(p,'PlotCal',true); 
addParameter(p,'PreprocessingX',2);
addParameter(p,'PreprocessingY',2);
addParameter(p,'ObsLabel',[]);
addParameter(p,'ObsClass',[]);
addParameter(p,'BlurIndex',Inf);
parse(p,varargin{:});


% Extract inputs from inputParser for code legibility
test = p.Results.ObsTest;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
lvs = p.Results.LVs;
plottype = p.Results.PlotType;
plotcal = p.Results.PlotCal;
label = p.Results.ObsLabel;
classes = p.Results.ObsClass;
blur = p.Results.BlurIndex;

L = size(test, 1);

if plotcal
    K = N+L;
else
    K = L;
end

if isempty(label) 
    if plotcal
        label = [1:N 1:L]; 
    else
        label = 1:L;
    end
end
if isempty(classes)
    if plotcal
        classes = [ones(N,1);2*ones(L,1)];  
    else 
        classes = ones(L,1); 
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
assert (A>0, 'Dimension Error: parameter ''LVs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: parameter ''ObsTest'' must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: parameter ''ObsLabel'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: parameter ''ObsClass'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[xcs,m,sd] = preprocess2D(x,'Preprocessing',prepx);
ycs = preprocess2D(y,'Preprocessing',prepy);

model = simpls(xcs,ycs,'LVs',lvs); 
R = model.altweights;
T = xcs*R;

if ~isempty(test)
    testcs = preprocess2Dapp(test,m,'Scale',sd);
    TT = testcs*R;
else
    TT = [];
end


%% Show results

if plotcal
    Tt = [T;TT];
else
    Tt = TT;
end

if strcmp(plottype,'Bars') || A == 1
    for i=1:length(lvs)
        plotVec(Tt(:,i), 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('Scores LV %d (%.0f%%)',lvs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2)))});
    end
elseif strcmp(plottype,'Scatter')
    for i=1:length(lvs)-1
        for j=i+1:length(lvs)
            plotScatter([Tt(:,i),Tt(:,j)], 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{sprintf('Scores LV %d (%.0f%%)',lvs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2))),sprintf('Scores LV %d (%.0f%%)',lvs(j),100*sum(T(:,j).^2)/sum(sum(xcs.^2)))}','BlurIndex',blur);
        end
    end
end
        