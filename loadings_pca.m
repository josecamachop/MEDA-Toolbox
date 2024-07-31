
function P = loadings_pca(x,varargin)


% Compute and plot loadings in PCA. This routine is deprecated and superseded 
% by loadings.m (please, use the latter)
%
% P = loadings_pca(x) % minimum call
% P = loadings_pca(x,'Pcs',pcs,'Preprocessing',prep,'Option',opt,'VarsLabel',label,'ObsClass',classes,'BlurIndex',blur) % complete call
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% Optional INPUTS (parameters):
%
% 'Pcs': [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 1:rank(xcs)
%
% 'Preprocessing': [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
%
% 'Option': (str or num) options for data plotting: binary code of the form 'ab' for:
%       a:
%           0: no plots
%           1: plot loadings
%       b:
%           0: scatter plot of pairs of PCs
%           1: bar plot of each single PC
%   By deafult, opt = '10'. If less than 2 digits are specified, least 
%   significant digit is set to 0, i.e. opt = 1 means a=1 and b=0. If a=0, 
%   then b is ignored.
%
% 'VarsLabel': [Mx1] name of the variables (numbers are used by default)
%
% 'ObsClass': [Mx1] groups for different visualization (a single group 
%   by default)
%
% 'BlurIndex': [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default).
%
%
% OUTPUTS:
%
% P: [MxA] loadings
%
%
% EXAMPLE OF USE: Scatter plot of random scores
%
% A = cell(1, 10);
% 
% for i = 1:10
%     A{i} = ['A_{', num2str(i), '}'];
% end
% 
% X = simuleMV(20,10,'LevelCorr',8);
% P = loadings_pca(X,'Pcs',1:3,'VarsLabel',A);
%
%
% EXAMPLE OF USE: Line plot of random scores
%
% X = real(ADICOV(randn(10,10).^19,randn(100,10),10));
% P = loadings_pca(X,'Pcs',1:3,'Option',11);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 22/ Apr/2024.
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
N = size(x, 1);
M = size(x, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
PCS = 1:rank(x);
addParameter(p,'Pcs',PCS);  
addParameter(p,'Preprocessing',2);
addParameter(p,'Option',10);  
addParameter(p,'VarsLabel',1:M);
addParameter(p,'ObsClass',ones(M,1));   
addParameter(p,'BlurIndex',1);     
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
pcs = p.Results.Pcs;
prep = p.Results.Preprocessing;
opt = p.Results.Option;
label = p.Results.VarsLabel;
classes = p.Results.ObsClass;
blur = p.Results.BlurIndex;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'0'); end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
A = length(pcs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''Pcs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''Pcs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==2, 'Dimension Error: parameter ''Option'' must be a string or num of 2 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: parameter ''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: parameter ''ObsClass'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''Pcs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,'Preprocessing',prep);
[P,T] = pca_pp(xcs,'Pcs',pcs);


%% Show results

if opt(1) == '1'
    
    if length(pcs) == 1 || opt(2) == '1'
        for i=1:length(pcs)
                plot_vec(P(:,i), 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('Loadings PC %d (%.0f%%)',pcs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2)))});
        end
    else
        for i=1:length(pcs)-1
            for j=i+1:length(pcs)
                plot_scatter([P(:,i),P(:,j)], 'EleLabel',label,'ObsClass' ,classes, 'XYLabel',{sprintf('Loadings PC %d (%.0f%%)',pcs(i),100*sum(T(:,i).^2)/sum(sum(xcs.^2))),sprintf('Loadings PC %d (%.0f%%)',pcs(j),100*sum(T(:,j).^2)/sum(sum(xcs.^2)))}','BlurIndex',blur);
            end      
        end
    end
end
        