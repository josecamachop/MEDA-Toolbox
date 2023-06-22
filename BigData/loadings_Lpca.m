function [P,fig_h] = loadings_Lpca(Lmodel,opt,blur)

% Compute and plot loadings in PCA for large data.
%
% loadings_Lpca(Lmodel) % minimum call
% loadings_Lpca(Lmodel,opt,blur) % complete call
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.lvs: (1x1) number of PCs.
%       Lmodel.vclass: [Mx1] class associated to each variable.
%       Lmodel.var_l: {ncx1} label of each variable.
%
% opt: (str or num) options for data plotting: binary code of the form 'ab' for:
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
% blur: [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default)
%
%
% OUTPUTS:
%
% P: [MxA] scores
%
% fig_h: (Lx1) figure handles
%
%
% EXAMPLE OF USE: Scatter plot of random scores
%
% X = simuleMV(20,10,8);
% Lmodel = Lmodel_ini(X);
% Lmodel.lvs = 1:3;
% P = loadings_Lpca(Lmodel);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 21/Jun/2023
%
% Copyright (C) 2023  University of Granada, Granada
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
[ok, Lmodel] = check_Lmodel(Lmodel);
if nargin < 2 || isempty(opt), opt = '10'; end; 

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'0'); end

if nargin < 3 || isempty(blur),    blur    = 1;       end;

% Validate dimensions of input data
assert (ischar(opt) && length(opt)==2, 'Dimension Error: 2nd argument must be a string or num of maximum 2 bits. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3nd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

P = Lpca(Lmodel);


%% Show results

fig_h = [];
if opt(1) == '1'
    
    t_var = var_Lpca(Lmodel,0);
    
    if length(Lmodel.lvs) == 1 || opt(2) == '1'
        for i=1:length(Lmodel.lvs)
                fig_h(i) = plot_vec(P(:,i), Lmodel.var_l, Lmodel.vclass, {'',sprintf('Loadings PC %d (%.0f%%)',Lmodel.lvs(i),100*(t_var(i) - t_var(i+1)))});
        end
    else
        h = 1;
        for i=1:length(Lmodel.lvs)-1
            for j=i+1:length(Lmodel.lvs)
                fig_h(h) = plot_scatter([P(:,i),P(:,j)], Lmodel.var_l, Lmodel.vclass, {sprintf('Loadings PC %d (%.0f%%)',Lmodel.lvs(i),100*(t_var(i) - t_var(i+1))),sprintf('Loadings PC %d (%.0f%%)',Lmodel.lvs(j),100*(t_var(j) - t_var(j+1)))}',[],[],[],[],blur);
                h = h+1;
            end      
        end
    end
end
        