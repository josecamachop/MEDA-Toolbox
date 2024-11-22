function [P,figH] = loadingsLpca(Lmodel,varargin)

% Compute and plot loadings in PCA for large data.
%
% loadingsLpca(Lmodel) % minimum call
% loadingsLpca(Lmodel,opt,blur) % complete call
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.lvs: (1x1) number of PCs.
%       Lmodel.vclass: [Mx1] class associated to each variable.
%       Lmodel.varl: {ncx1} label of each variable.
%
% Optional INPUTS (parameters):
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
% 'BlurIndex': [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default)
%
%
% OUTPUTS:
%
% P: [MxA] scores
%
% figH: (Lx1) figure handles
%
%
% EXAMPLE OF USE: Scatter plot of random scores
%
% X = simuleMV(20,10,'LevelCorr',8);
% Lmodel = iniLmodel(X);
% Lmodel.lvs = 1:3;
% P = loadingsLpca(Lmodel);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 22/Nov/2024
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
[ok, Lmodel] = checkLmodel(Lmodel);
% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Option','10');
addParameter(p,'BlurIndex','1');  
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
opt = p.Results.Option;
blur = p.Results.BlurIndex;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'0'); end


% Validate dimensions of input data
assert (ischar(opt) && length(opt)==2, 'Dimension Error: 2nd argument must be a string or num of maximum 2 bits. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3nd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

Lmodel = Lpca(Lmodel);
P = Lmodel.loads;


%% Show results

figH = [];
if opt(1) == '1'
    
    tvar = varLpca(Lmodel,'Option',0);
    
    if length(Lmodel.lvs) == 1 || opt(2) == '1'
        for i=1:length(Lmodel.lvs)
                figH(i) = plotVec(P(:,i),'EleLabel',Lmodel.varl,'ObsClass',Lmodel.vclass,'XYLabel',{'',sprintf('Loadings PC %d (%.0f%%)',Lmodel.lvs(i),100*(tvar(i) - tvar(i+1)))});
        end
    else
        h = 1;
        for i=1:length(Lmodel.lvs)-1
            for j=i+1:length(Lmodel.lvs)
                figH(h) = plotScatter([P(:,i),P(:,j)],'EleLabel',Lmodel.varl,'ObsClass',Lmodel.vclass,'XYLabel',{sprintf('Loadings PC %d (%.0f%%)',Lmodel.lvs(i),100*(tvar(i) - tvar(i+1))),sprintf('Loadings PC %d (%.0f%%)',Lmodel.lvs(j),100*(tvar(j) - tvar(j+1)))}','BlurIndex',blur);
                h = h+1;
            end      
        end
    end
end
        