
function cumpress = ckf(xcs,T,P,opt)

% CKF Algorithm: Journal of Chemometrics, 29(8): 467-478, 2015
%
% ckf(xcs,T,P) % minimum call
% ckf(xcs,T,P,opt) % complete call
%
%
% INPUTS:
%
% xcs: (LxM) billinear data set preprocessed
%
% T: (LxA) scores.
%
% P: (MxA) loadings.
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: plot (default)
%
%
% OUTPUTS:
%
% cumpress: ((maxpcs+1)x1) ckf curve.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 07/Sep/15.
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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

if nargin < 3, error('Error in the number of arguments.'); end;
if nargin < 4, opt = 1; end;

maxpcs = size(P,2);

%% Main code

cumpress = zeros(maxpcs+1,1);
press = zeros(maxpcs+1,size(P,1));

s = size(xcs);

for i=0:maxpcs,
    
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
    
    term1 = sum(term1_p.^2,1);
    term2 = sum(2*term1_p.*term3_p,1);
    term3 = sum(term3_p.^2,1);
    
    press(i+1,:) = sum([term1;term2;term3]);
    
    cumpress(i+1) = sum(press(i+1,:));
    end
    
%% Show results

if opt == 1,
        fig_h = plot_vec(cumpress/cumpress(1),num2str((0:maxpcs)')','ckf',[],1,'r--');
end

        