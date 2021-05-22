function x_var = var_Lpca(Lmodel,opt)

% Variability captured in terms of the number of PCs.
%
% var_Lpca(Lmodel) % minimum call
% x_var = var_Lpca(Lmodel,opt) %complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.lvs: (1x1) number of PCs.
%
% opt: (str or num) options for data plotting.
%       0: no plots.
%       1: bar plot (default)
%
%
% OUTPUTS:
%
% x_var: [Ax1] Percentage of captured variance of X.
%
%
% EXAMPLE OF USE: Random data
%
% Lmodel = Lmodel_ini(simuleMV(20,10,8));
% Lmodel.lvs = 0:10;
% x_var = var_Lpca(Lmodel);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 30/Oct/2016
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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

check_Lmodel(Lmodel);

% Preprocessing
Lmodel.lvs = unique([0 Lmodel.lvs]);

if nargin < 2 || isempty(opt), opt = '1'; end;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Validate dimensions of input data
assert (ischar(opt) && length(opt)==1, 'Dimension Error: 2nd argument must be a string or num of 1 bit. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 2nd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

Lmodel.lvs = 0:max(Lmodel.lvs);
P = Lpca(Lmodel);

totalVx = sum(eig(Lmodel.XX));
x_var = ones(max(Lmodel.lvs)+1,1);
for i=1:max(Lmodel.lvs),
    x_var(i+1) = x_var(i+1) - sum(eig(P(:,1:i)'*Lmodel.XX*P(:,1:i)))/totalVx;
end
    
%% Show results

if opt == '1',
    plot_vec(x_var,Lmodel.lvs,[],{'#PCs','% Residual Variance'},[],0);
end
