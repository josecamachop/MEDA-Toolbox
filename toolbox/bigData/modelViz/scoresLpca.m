function [T,TT,figH] = scoresLpca(Lmodel,varargin)

% Compute and plot compressed scores in PCA for large data. The original 
% paper is Camacho J. Visualizing Big data with Compressed Score Plots: 
% Approach and Research Challenges. Chemometrics and Intelligent Laboratory
% Systems, 2014, 135: 110-125.
%
% scoresLpca(Lmodel) % minimum call
% [T,TT] = scoresLpca(Lmodel,test,opt,label,classes) % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix
%       Lmodel.lvs: [1x1] number of PCs
%       Lmodel.sdT: [1xA] variance of PCs
%       Lmodel.centr: (LxM) centroids of the clusters of observations
%       Lmodel.multr: (Lx1) multiplicity of each cluster
%       Lmodel.class: (Lx1) class associated to each cluster
%       Lmodel.av: [1xM] sample average according to the preprocessing method
%       Lmodel.sc: [1xM] sample scale according to the preprocessing method
%       Lmodel.weight: [1xM] weight applied after the preprocessing method
%
% Optional INPUTS (parameters):
%
% 'Test': [LxM] data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data
%
% 'PlotType': str
%      'Scatter': scatterplot (by default for pairs of PCs)
%      'Bars': bar plot (by default for a single PC)
%
% 'PlotCal': bool
%      false: plot only test data
%      true: plot both calibration and test (by default)
%
% 'PlotMult': str
%      'none': do not plot multiplicity 
%      'size': plot multiplicity info in the size of the markers (by default)
%      'shape': plot multiplicity info in the shape of the markers
%      'zaxis': plot multiplicity information in the Z axis
%      'zsize': plot multiplicity info in the size of the markers and
%               classes in Z-axis
%           
%   By deafult, opt = '00000'. If less than 5 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0, c=0, d=00 
%
% 'ObsLabel': [Lx1] name of the test observations (numbers are used by default)
%
% 'ObsClass': [Lx1] groups in test for different visualization (a single group 
%   by default)
%
% 'BlurIndex': [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default).
%
%
%
% OUTPUTS:
%
% T: [LxA] calibration scores
%
% TT: [NxA] test scores
%
% figH: (Lx1) figure handles
%
%
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(20,10,'LevelCorr',8);
% Lmodel = iniLmodel(X);
% Lmodel.lvs = 1:3;
% T = scoresLpca(Lmodel);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% nobs = 100;
% nvars = 10;
% X = simuleMV(nobs,nvars,'LevelCorr',8);
% Lmodel = iniLmodel(X);
%
% nobst = 10;
% test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(X)*(nobst-1)/(nobs-1));
%
% Lmodel.multr = ceil(exp(10*randn(nobs,1)));
%
% Lmodel.lvs = 1;
% scoresLpca(Lmodel,'Test',test);
%
% Lmodel.lvs = 1:2;
% [T,TT] = scoresLpca(Lmodel,'Test',test);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 21/Jan/2025
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
  

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

checkLmodel(Lmodel);

N = size(Lmodel.centr,1);
M = size(Lmodel.XX, 2);
% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Test',[]);
addParameter(p,'ObsLabel',[]);
addParameter(p,'ObsClass',[]);
addParameter(p,'PlotMult','size'); 
addParameter(p,'PlotType','Scatter');
addParameter(p,'PlotCal',true);
addParameter(p,'BlurIndex',1);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
test = p.Results.Test;
label = p.Results.ObsLabel;
classes = p.Results.ObsClass;
plottype = p.Results.PlotType;
plotcal = p.Results.PlotCal;
plotmult = p.Results.PlotMult;
blur = p.Results.BlurIndex;

L = size(test, 1);
A = length(Lmodel.lvs);

if plotcal
    K = N+L;
else
    K = L;
end

if isempty(label) 
    if plotcal
        if isempty(Lmodel.obsl)
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
        label = cellstr(num2str((1:L)'));
    end
else
    if plotcal
        if isempty(Lmodel.obsl)
            lb1 = cellstr(num2str((1:N)'));
            label = {lb1{:} label{:}};
        else
            label = {Lmodel.obsl{:} label{:}};
        end
    end
end

label(find(ismember(label, 'mixed'))) = {''};

if isempty(classes)
    if plotcal
        classes = [Lmodel.class;2*ones(L,1)];
    else 
        classes = ones(L,1); 
    end
elseif plotcal && length(classes)==L
        classes = [Lmodel.class;2*classes];
end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
assert (A>0, 'Dimension Error: 1sr argument with non valid content. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: 2nd argument must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (isequal(size(label), [K 1]), 'Dimension Error: 4th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: 5th argument must be L-by-1. Type ''help %s'' for more info.', routine(1).name); 


%% Main code

Lmodel = Lpca(Lmodel);
P = Lmodel.loads;
T = Lmodel.centr*P;

if ~isempty(test)
    testcs = preprocess2Dapp(test,Lmodel.av,'Scale',Lmodel.sc,'Weight',Lmodel.weight);
    TT = testcs*P;
else
    TT = [];
end

%% Show results

figH = [];
     
if plotcal
    ttt = [T;TT];
    mult = [Lmodel.multr;ones(size(TT,1),1)];
else
    ttt = TT;
    mult = ones(size(TT,1));
end

tvar = varLpca(Lmodel,'Plot',false);

indM = floor(log10((max(Lmodel.multr)+1)));
indm = floor(log10((min(Lmodel.multr)+1)));
markers = 10.^[indm indm+(indM-indm)/2 indM];
if sum(markers)<100, markers = []; end

if length(Lmodel.lvs) == 1 || strcmp(plottype,'Bars')
    for i=1:length(Lmodel.lvs)
        figH(i) = plotVec(ttt(:,i),'EleLabel',label,'ObsClass',classes,'XYLabel',{'',sprintf('Compressed Scores PC %d (%.0f%%)',Lmodel.lvs(i),100*(tvar(Lmodel.lvs(i)) - tvar(Lmodel.lvs(i)+1)))},'Multiplicity',mult,'Markers',markers);
    end
elseif strcmp(plottype,'Scatter')
    h = 1;
    for i=1:length(Lmodel.lvs)-1
        for j=i+1:length(Lmodel.lvs)
            figH(h) = plotScatter([ttt(:,i),ttt(:,j)],'EleLabel',label,'ObsClass',classes,'XYLabel',{sprintf('Scores PC %d (%.0f%%)',Lmodel.lvs(i),100*(tvar(Lmodel.lvs(i)) - tvar(Lmodel.lvs(i)+1))),sprintf('Scores PC %d (%.0f%%)',Lmodel.lvs(j),100*(tvar(j) - tvar(j+1)))},'PlotMult',plotmult,'Multiplicity',mult,'Markers',markers,'BlurIndex',blur);
            h = h+1;
        end
    end
end



