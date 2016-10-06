function [cumpress,press,term1,term2,term3] = ckf(xcs,T,P,opt)

% CKF Algorithm: Journal of Chemometrics, 29(8): 467-478, 2015
%
% cumpress = ckf(xcs,T,P) % minimum call
% [cumpress,press,term1,term2,term3] = ckf(xcs,T,P,opt) % complete call
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% T: [NxA] scores.
%
% P: [MxA] loadings.
%
% opt: (str or num) options for data plotting.
%       0: no plots.
%       1: plot (default)
%
%
% OUTPUTS:
%
% cumpress: [Ax1] ckf curve.
%
% press: [AxM] PRESS per variable.
%
% term1: [NxM] first error term, according to referred paper.
%
% term2: [NxM] second error term, according to referred paper.
%
% term3: [NxM] third error term, according to referred paper.
%
%
% EXAMPLE OF USE: Random curve
%
% X = simuleMV(20,10,8);
% [P,T] = pca_pp(X);
% cumpress = ckf(X,T,P);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Apr/2016
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
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(xcs, 1);
M = size(xcs, 2);
A = size(T, 2);
if nargin < 4 || isempty(opt), opt=1; end;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Validate dimensions of input data
assert (isequal(size(T), [N A]), 'Dimension Error: 2nd argument must be N-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(P), [M A]), 'Dimension Error: 3rd argument must be M-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(opt), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3rd argument must contain a binary value. Type ''help %s'' for more info.', routine(1).name);


%% Main code

cumpress = zeros(A+1,1);
press = zeros(A+1,size(P,1));

s = size(xcs);

for i=0:A,
    
    if i > 0, % PCA Modelling
        p2 = P(:,1:i);
        srec = T(:,1:i)*p2';
        erec = xcs - srec;
        term3_p = erec;
        term1_p = xcs.*(ones(s(1),1)*(sum(p2.*p2,2))');
    else % Modelling with the average
        term1_p = zeros(size(xcs));
        term3_p = xcs;
    end
    
    term1 = term1_p.^2;
    term2 = 2*term1_p.*term3_p;
    term3 = term3_p.^2;
    
    press(i+1,:) = sum([sum(term1,1);sum(term2,1);sum(term3,1)]);
    
    cumpress(i+1) = sum(press(i+1,:));
end
    
    
%% Show results

if opt == '1',    
    fig_h = plot_vec(cumpress/cumpress(1),0:A,[],{'#PCs','ckf'},[],0); 
end

        