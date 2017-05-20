function [p,t,bel] = Lgpca(Lmodel,states,opt)

% Group-wise Principal Component Analysis for large data.
%
% [p,t,bel] = Lgpca(Lmodel,states)     % minimum call
% [p,t,bel] = Lgpca(Lmodel,states,opt)     % complete call
%
% INPUTS:
%

% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.lv: (1x1) number of PCs A.
%
% states: {Sx1} Cell with the groups of variables.
%
% opt: options
%   - 0: set the number of pcs to Lmodel.lv (by default)
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
% last modification: 05/Sep/15.
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

% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if nargin < 3, opt=0; end;

% Main code

map = Lmodel.XX;
x = Lmodel.centr;
I =  eye(size(map));
B = I;
j=1;
if opt,
    finish = false;
else
    finish = (j > max(Lmodel.lvs));
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
         finish = (j > max(Lmodel.lvs));
    end
    
end