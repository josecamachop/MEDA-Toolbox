
function [P,W,Q] = loadings_pls(x,y,varargin)

% Compute and plot loadings in PLS. This routine is deprecated and superseded 
% by loadings.m (please, use the latter)
%
% P = loadings_pls(x,y) % minimum call
% [P,W,Q] = loadings_pls(x,y,'LatVars',lvs,'PreprocessingX',prepx,'PreprocessingY',prepy,'Option',opt,'VarsLabel',label,'ObsClass',classes,'BlurIndex',blur) % complete call
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
% Optional INPUTS (parameters):
%
% 'LatVars': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 1:rank(x)
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
% 'Option': (str or num) options for data plotting: binary code of the form 'abc' for:
%       a:
%           0: no plots
%           1: plot loadings
%       b:
%           0: scatter plot of pairs of LVs
%           1: bar plot of each single LV
%       c:
%           0: plot weights
%           1: plot X-block loadings
%   By deafult, opt = '100'. If less than 3 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0 and c=0. 
%   If a=0, then b and c are ignored.
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
% P: [MxA] X-block loadings
%
% W: [MxA] X-block weights
%
% Q: [OxA] Y-block loadings
%
%
% EXAMPLE OF USE: Random loadings: bar and scatter plot of loadings
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% A = cell(1, 10);
% 
% for i = 1:10
%     A{i} = ['A_{', num2str(i), '}'];
% end
% 
% loadings_pls(X,Y,'LatVars',1);
% [P,W,Q] = loadings_pls(X,Y,'LatVars',1:3,'VarsLabel',A);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 22/Apr/24
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
O = size(y, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
LVS = 1:rank(x);
addParameter(p,'LatVars',LVS);  
addParameter(p,'PreprocessingX',2);
addParameter(p,'PreprocessingY',2);
addParameter(p,'Option','100');  
addParameter(p,'VarsLabel',1:M);
addParameter(p,'ObsClass',ones(M,1));   
addParameter(p,'BlurIndex',1);     
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LatVars;
prepx = p.Results.PreprocessingX;
prepy = p.Results.PreprocessingY;
opt = p.Results.Option;
label = p.Results.VarsLabel;
classes = p.Results.ObsClass;
blur = p.Results.BlurIndex;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Convert int arrays to str
if length(opt)<2, opt = strcat(opt,'00'); end
if length(opt)<3, opt = strcat(opt,'0'); end

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
assert (isequal(size(prepx), [1 1]), 'Dimension Error: parameter ''PreprocessingX'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: parameter ''PreprocessingY'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==3, 'Dimension Error: parameter ''Option'' must be a string or num of 3 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: parameter ''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: parameter ''ObsClass'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LatVars'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,'Preprocessing',prepx);
ycs = preprocess2D(y,'Preprocessing',prepy);

[beta,W,P,Q] = simpls(xcs,ycs,'LatVars',lvs); 

%% Show results

if opt(1) == '1'
    
    if opt(3) == '0'
        Pt = W;
        text = 'Weights';
    else
        Pt = P;
        text = 'X-block loadings';
    end
    
    if length(lvs) == 1 || opt(2) == '1'
        for i=1:length(lvs),
            plot_vec(Pt(:,i),  'EleLabel',label, 'ObsClass',classes, 'XYLabel',{'',sprintf('%s LV %d',text,lvs(i))});
        end
    else
        for i=1:length(lvs)-1,
            for j=i+1:length(lvs),
                plot_scatter([Pt(:,i),Pt(:,j)], 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{sprintf('%s LV %d',text,lvs(i)),sprintf('%s LV %d',text,lvs(j))}','BlurIndex',blur);
            end      
        end
    end
end
        