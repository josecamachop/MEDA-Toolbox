function X = simuleMV(obs,vars,lcorr,corM)

% Simulation of MV data with ADICOV, submitted to Chemolab
%
% X = simuleMV(obs,vars) % minimum call
% X = simuleMV(obs,vars,lcorr,corM)% complete call
%
%
% INPUTS:
%
% obs: [1x1] number of observations (rows) in the output.
%
% vars: [1x1] number of variables (columns) in the output.
%
% lcorr: [1x1] level of correlation among variables in [0,10] (5 by default) 
%
% corM: [vars x vars] covariance for simulation (empty by default) 
%
%
% OUTPUTS:
%
% X: [obs x vars] data matrix generated.
%
%
% EXAMPLE OF USE: To obtain a matrix 100x10 with random covariance matrix, 
% use the following call:
%
% X = simuleMV(100,10,6);
% plot_map(corr(X)); % visualization 
% var_pca(X)
%
%
% EXAMPLE OF USE: To obtain a matrix 100x10 with random covariance matrix 
%   of complete rank:
%
% X = simuleMV(100,10,8) + 0.1*randn(100,10);
% plot_map(corr(X)); % visualization 
% var_pca(X)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 21/Sep/16.
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if nargin < 3 || isempty(lcorr), lcorr = 5; end;
if nargin < 4 || isempty(corM), 
    uselevel = true; 
    corM = eye(vars); 
else
    uselevel = false;
end;
    
% Validate dimensions of input data
assert (isequal(size(obs), [1 1]), 'Dimension Error: 1st argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(vars), [1 1]), 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lcorr), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(corM), [vars vars]), 'Dimension Error: 4th argument must be vars-by-vars. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (obs>0, 'Value Error: 1st argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(obs), obs), 'Value Error: 1st argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (vars>0, 'Value Error: 2nd argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(vars), vars), 'Value Error: 2nd argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (lcorr >= 0, 'Value Error: 3rd argument must be above or equal to 0. Type ''help %s'' for more info.', routine(1).name);
assert (lcorr<=10, 'Value Error: 3rd argument must be equal to or below 10. Type ''help %s'' for more info.', routine(1).name);


%% Main code

if uselevel,
    if lcorr>0,
        x = [0.2837   17.9998   19.3749    2.0605    0.3234 0.4552   12.0737   16.6831    5.2423    0.5610 0.3000   14.7440 4.4637e+04 7.1838    0.8429];
        obs2 = round(Floc(x,[lcorr,vars]).^2);
        obs2 = max(2,obs2);
        if obs < obs2,
            disp('Warning: correlation level too low. Resulting matrix may show a higher correlation due to structural constraints.')
        end
        X = real(ADICOV(eye(vars),randn(obs2,vars),vars));
        COV = corr(X);
        corM = COV + 0.01*eye(vars);
    else
        if obs < vars,
            disp('Warning: correlation level too low. Resulting matrix may show a higher correlation due to structural constraints.')
        end
    end
end

X = real(ADICOV(corM,randn(obs,vars),vars));


function y = Floc(x,xdata)

y = zeros(size(xdata,1),1);

for i=1:size(xdata,1),
    
    switch xdata(i,1),
    
        case {1, 2, 3},
            y(i) = (x(1)*(x(2)-xdata(i,1)).*((log(xdata(i,2))/log(x(3))).^(x(4)*exp(-x(5)*xdata(i,1)))));
        case {4, 5, 6, 7},
            y(i) = (x(6)*(x(7)-xdata(i,1)).*((log(xdata(i,2))/log(x(8))).^(x(9)*exp(-x(10)*xdata(i,1)))));
        case {8, 9, 10},
            y(i) = (x(11)*(x(12)-xdata(i,1)).*((log(xdata(i,2))/log(x(13))).^(x(14)*exp(-x(15)*xdata(i,1)))));
            
    end

end;
