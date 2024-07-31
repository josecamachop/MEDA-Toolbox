
function fig_h =loadings(model,varargin)

% Compute and plot loadings.
%
% fig_h =loadings(model) % minimum call
% fig_h =loadings(model,'Option'opt,'Title',tit,'VarsLabel',label,'ObsClass,classes,'BlurIndex',blur) % complete call
%
% INPUTS:
%
% model (structure): structure with model parameters. 
%   lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%       first two LVs).
%   loads: [MxA] model parameters.
%   scores: [NxA] data scores. 
%
% Optional Inputs (parameter):
%
% 'Option': (str) options for data plotting: binary code of the form 'ab' for:
%       a:
%           0: scatter plot of pairs of LVs 
%           1: bar plot of each single LV
%       b:
%           0: plot for categorical classes (consistent with a legend)
%           1: plot for numerical classes (consistent with a colorbar) 
%
%   By deafult, opt = '00'. If less than 2 digits are specified, the least 
%   significant digit is set to 0, i.e. opt = 1 means a=1, b=0.
%
% 'Title': (str) title for the plots. Empty by default;
%
% 'VarsLabel': [Mx1] name of the variables (numbers are used by default)
%
% 'ObsClass': [Mx1] groups for different visualization (a single group 
%   by default)
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
% EXAMPLE OF USE: Random data
% 
% A = cell(1, 10);
% 
% for i = 1:10
%     A{i} = ['A_{', num2str(i), '}'];
% end
% 
% X = simuleMV(20,10,'LevelCorr',8);
% [~,~,model] = pca_pp(X,'Pcs',1:2);
% 
% P = loadings(model,'VarsLabel',A);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 22/Apr/2024
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
M = size(model.loads, 1);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Option',00);  
addParameter(p,'Title',' ');
addParameter(p,'VarsLabel',1:M);
addParameter(p,'ObsClass',ones(M,1));   
addParameter(p,'BlurIndex',1);     
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
opt = p.Results.Option;
tit = p.Results.Title;
label = p.Results.VarsLabel;
classes = p.Results.ObsClass;
blur = p.Results.BlurIndex;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
while length(opt)<2, opt = strcat(opt,'0'); end
if opt(2) == '0', opt(2) = '1'; else,  opt(2) = '0'; end


% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
assert (ischar(opt) && length(opt)==2, 'Dimension Error: parameter ''Option'' must be a string or num of 2 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: parameter''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: parameter ''ObsClass'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

P = model.loads;


%% Show results

if ~isfield(model,'type') || strcmp(model.type,'PCA')
    dim = 'PC';
elseif strcmp(model.type,'PLS')
    dim = 'LV'
else
    dim = 'PC';
end

fig_h = [];
if length(model.lvs) == 1 || opt(1) == '1'
    for i=1:length(model.lvs)
        fig_h = [fig_h plot_vec(P(:,i), 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('Loadings %s %d',dim,model.lvs(i))})];
        title(tit);
    end
else
    for i=1:length(model.lvs)-1
        for j=i+1:length(model.lvs)
            fig_h = [fig_h plot_scatter([P(:,i),P(:,j)], 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{sprintf('Loadings %s %d',dim,model.lvs(i)),sprintf('Loadings %s %d',dim,model.lvs(j))}','Option',opt(2),'BlurIndex',blur)];
            title(tit);
        end
    end
end
        