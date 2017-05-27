
function L = leverages_Lpls(Lmodel,opt)

% Compute and plot the leverages of variables in PLS for large data.
%
% L = leverages_Lpls(Lmodel) % minimum call
% L = leverages_Lpls(Lmodel,opt) % complete call
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PLS
%   model:
%       Lmodel.XX: [MxM] X-block cross-product matrix.
%       Lmodel.XY: [MxO] cross-product matrix between the x-block and the
%           y-block.
%       Lmodel.lvs: [1x1] number of Latent Variables.
%       Lmodel.vclass: [Mx1] class associated to each variable.
%       Lmodel.var_l: {ncx1} label of each variable.
%
% opt: (str or num) options for data plotting
%       0: no plots.
%       1: plot bar plot of leverages (default) 
%
%
% OUTPUTS:
%
% L: [Mx1] leverages of the variables
%
%
% EXAMPLE OF USE: Random leverages
%
% X = simuleMV(20,10,8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Lmodel = Lmodel_ini(X,Y);
% Lmodel.lvs = 1:3;
% L = leverages_Lpls(Lmodel);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/May/17.
%
% Copyright (C) 2017  University of Granada, Granada
% Copyright (C) 2017  Jose Camacho Paez
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
if nargin < 2 || isempty(opt), opt = 1; end; 

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Validate dimensions of input data
assert (~isempty(Lmodel.XY), 'Dimension Error: Empty XY. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 2nd argument must be a string or num of maximum 3 bits. Type ''help %s'' for more info.', routine(1).name);
  
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 2nd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

[beta,W,P,Q,R] = Lpls(Lmodel);

L = diag(W*W');


%% Show results

if opt == '1', 
    plot_vec(L, Lmodel.var_l, Lmodel.vclass, {'Variables','Leverages'});
end        