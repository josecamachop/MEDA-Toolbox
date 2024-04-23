
function fig_h = scores(model,varargin)

% Compute and plot scores.
%
% fig_h = scores(model) % minimum call
% fig_h = scores(model,'ObsTest',test,'Option',opt,'Title',tit,'ObsLabel',label,'ObsClass',classes,'BlurIndex',blur) % complete call
%
% INPUTS:
%
% model (structure): structure with model parameters. 
%   var: [1x1] Total variance
%   lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%       first two LVs).
%   av: [1xM] centering parameters. 
%   sc: [1xM] scaling parameters. 
%   loads: [MxA] model parameters.
%   scores: [NxA] data scores.
%
% Optional INPUTS(parameters):
%
% 'ObsTest': [LxM] data set with the observations to be visualized in the model
%   space. By default, model.scores are plotted.
%
% 'Option': (str) options for data plotting: binary code of the form 'abcde' for:
%       a:
%           0: scatter plot of pairs of LVs 
%           1: bar plot of each single LV
%       b:
%           0: plot calibration and test data
%           1: plot only test data 
%       c:
%           0: plot for categorical classes (consistent with a legend)
%           1: plot for numerical classes (consistent with a colorbar) 
%       d:
%           0: do not plot multiplicity
%           1: plot multiplicity
%       e: (for d 0)
%           0: filled marks
%           1: empty marks
%       e: (for d 1)
%           00: plot multiplicity info in the size of the markers.
%           01: plot multiplicity info in the form of the markers.
%           10: plot multiplicity information in the Z axis.
%           11: plot multiplicity info in the size of the markers and
%               classes in Z-axis
%
%   By deafult, opt = '00000'. If less digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0, c=0, d=0
%   and e=0.
%
% 'Title': (str) title for the plots. Empty by default;
%
% 'ObsLabel': [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% 'ObsClass': [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
%
% 'BlurIndex': [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default)
%
%
% OUTPUTS:
%
% fig_h: set of figure handles
%
%
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(20,10,'LevelCorr',8);
% 
% model.lvs = 1:3;
% [Xcs,model.av,model.sc] = preprocess2D(X);
% model.var = trace(Xcs'*Xcs);
% [model.loads,model.scores] = pca_pp(Xcs,'Pcs',model.lvs);
% 
% scores(model);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% n_obs = 100;
% n_vars = 10;
% X = simuleMV(n_obs,n_vars,'LevelCorr',8);
% 
% model.lvs = 1:2;
% [Xcs,model.av,model.sc] = preprocess2D(X);
% model.var = trace(Xcs'*Xcs);
% [model.loads,model.scores] = pca_pp(Xcs,'Pcs',model.lvs);
% 
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,'LevelCorr',6,'Covar',corr(X)*(n_obst-1)/(n_obs-1));
% 
% scores(model,'ObsTest',test);

%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(model.scores, 1);
M = size(model.loads,1);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'ObsTest',[]);
addParameter(p,'Option','00000');
addParameter(p,'Title',' ');  
addParameter(p,'BlurIndex',1);
addParameter(p,'ObsLabel',[]);
addParameter(p,'ObsClass',[]);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility

test = p.Results.ObsTest;
opt = p.Results.Option;
tit = p.Results.Title;
blur = p.Results.BlurIndex;
label = p.Results.ObsLabel;
classes = p.Results.ObsClass;

L = size(test, 1);

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
while length(opt)<5, opt = strcat(opt,'0'); end
if length(opt)<6 && opt(4)==1, opt = strcat(opt,'0'); end

if opt(3) == '0', opt(3) = '1'; else,  opt(3) = '0'; end
if opt(2) == 1 || opt(2) == '1'
    K = L;
else
    K = N+L;
end


if isempty(label) 
    if opt(2) == 1 || opt(2) == '1'
        label = 1:L;
    else
        label = [1:N 1:L]; 
    end
end
if isempty(classes)
    if opt(2) == 1 || opt(2) == '1' 
        classes = ones(L,1); 
    else
        classes = [ones(N,1);2*ones(L,1)];  
    end
end



% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: parameter ''ObsTest'' must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (ischar(opt) && length(opt)>=5, 'Dimension Error: parameter ''Option'' must be a string or num of at least 5 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: parameter ''ObsLabel'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: parameter ''ObsClass'' must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

T = model.scores;
d = diag(T'*T);

if isfield(model,'scoresV')
    T = model.scoresV;
end

if ~isempty(test)
    testcs = preprocess2Dapp(test,model.av,'SDivideTest',model.sc);
    TT = testcs*model.loads;
else
    TT = [];
end


%% Show results

if opt(2) == '0'
    Tt = [T;TT];
else
    Tt = TT;
end

if ~isfield(model,'type') || strcmp(model.type,'PCA')
    dim = 'PC';
elseif strcmp(model.type,'PLS')
    dim = 'LV';
else
    dim = 'PC';
end

fig_h = [];
if length(model.lvs) == 1 || opt(1) == '1'
    for i=1:length(model.lvs)
        fig_h = [fig_h plot_vec(Tt(:,i), 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('Scores %s %d (%.0f%%)',dim,model.lvs(i),100*d(i)/model.var)})];
        title(tit);
    end
else
    for i=1:length(model.lvs)-1
        for j=i+1:length(model.lvs)
            fig_h = [fig_h plot_scatter([Tt(:,i),Tt(:,j)],'EleLabel',label,'ObsClass',classes,'XYLabel',{sprintf('Scores %s %d (%.0f%%)',dim,model.lvs(i),100*d(i)/model.var),sprintf('Scores %s %d (%.0f%%)',dim,model.lvs(j),100*d(j)/model.var)}','Option',opt(3:end),'BlurIndex',blur)];
            title(tit);
        end
    end
end
        