function [p,t,bel] = gpca(x,states,pc,opt)

% Group-wise Principal Component Analysis.
%
% [p,t,bel] = gpca(x,states,pc)     % minimum call
% [p,t,bel] = gpca(x,states,pc,opt)     % complete call
%
% INPUTS:
%
% x: (NxM) Two-way batch data matrix, N(observations) x M(variables)
%
% states: {Sx1} Cell with the groups of variables.
%
% pc: number of principal components. 
%
% opt: options
%   - 0: set the number of pcs to pc (by default)
%   - 1: extract at least 1 PC per state.
%
%
% OUTPUTS:
%
% p: (M x A) matrix of loadings.
%
% t: (N x A) matrix of scores.
%
% bel: (A x 1) correspondence between PCs and States.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 28/Aug/15.
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

% Parameters checking

if nargin < 3, error('Error in the number of arguments.'); end;
if nargin < 4, opt=0; end;

% Main code

map = x'*x;
I =  eye(size(map));
B = I;
j=1;
if opt,
    finish = false;
else
    finish = j <= pc,
end
while ~finish, 
    
    for i=1:length(states), % construct eigenvectors according to states
        map_aux = zeros(size(map));
        map_aux(states{i},states{i})= map(states{i},states{i});
        [V,D] = eig(map_aux);
        ind = find(diag(D)==max(diag(D)),1);
        R(:,i) = V(:,ind);
        S(:,i) = x*R(:,i);       
    end

    sS = sum(S.^2,1); % select pseudo-eigenvector with the highest variance
    ind = find(sS==max(sS),1);
    p(:,j) = R(:,ind);
    t(:,j) = S(:,ind);
    bel(j) = ind;
    
    q = B*R(:,ind); % deflate (Mackey'09)
    map = (I-q*q')*map*(I-q*q');
    x = x*(I-q*q');
    B = B*(I-q*q');
    
    j = j+1;
    
    if opt,
        if length(unique(bel))==length(states),
            finish = true;
        end
    else
        finish = j <= pc,
    end
    
end