function lim = spe_lim(res,p_value)

% Control limit for SPE statistic based on J.E. Jackson and G.S. Mudholkar. 
% Control procedures for residuals associated with principal component 
% analysis. Technometrics, 21:331-349, 1979.
%
% lim = spe_lim(res,p_value)        % complete call
%
% INPUTS:
%
% res: [NxM] Two-way residuals data matrix
%
% p_value: [1x1] p-value of the test, in (0,1]
%
%
% OUTPUTS:
%
% lim: [1x1] control limit at a 1-p_value confidence level.
%
%
% EXAMPLE OF USE: Compute the 99% confidence limit for 2 PCs and 100 
%   observations:
%
% X = simuleMV(100,10,'LevelCorr',8);
% Xcs = preprocess2D(X,'preprocessing',2);
% pcs = 1:2;
% [p,t] = pca_pp(Xcs,'Pcs',pcs);
% res = Xcs - t*p'; 
% lim = spe_lim(res,0.01)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Apr/2024
%
% Copyright (C) 2024 University of Granada, Granada
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(res, 1);
M = size(res, 2);

% Validate dimensions of input data
assert (isequal(size(p_value), [1 1]), 'Dimension Error: paramter ''p_value'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (p_value>=0 && p_value<1, 'Value Error: paramter ''p_value'' must be in (0,1]. Type ''help %s'' for more info.', routine(1).name);


%% Main code

pcs_left = rank(res);

if N>M
    lambda = eig(1/(N-1)*res'*res);
else
    lambda = eig(1/(N-1)*res*res');
end
[kk,ord]=sort(abs(lambda),'descend');
lambda = lambda(ord);

theta1 = sum(lambda(1:pcs_left));
theta2 = sum(lambda(1:pcs_left).^2);
theta3 = sum(lambda(1:pcs_left).^3);

h0 = 1-2*theta1*theta3/(3*theta2^2);

z = norminv(1-p_value);

lim = theta1*(z*sqrt(2*theta2*h0^2)/theta1 + 1 + theta2*h0*(h0-1)/(theta1^2))^(1/h0);


