function [centr,multn,labn,obslist] = psc(x,n_min,mult,lab,mat,obslist)

% Projection sequential clustering. 
%
% [centr,multn,labn] = psc(x,n_min)          % minimum call
% [centr,multn,labn] = psc(x,n_min,mult,lab,mat,obslist) % complete call
%
%
% INPUTS:
%
% x: (NxM) original matrix with centroids.
%
% n_min: (1x1) number of output clusters.
%
% mult: (Nx1) multiplicity of each cluster.
%
% lab: (Nx1) class associated to each cluster.
%
% mat: (MxA) projection matrix for distance computation.
%
% obslist: (Nx1) list of observations for the clustering file system.
%
%
% OUTPUTS:
%
% centr: (n_minxM) output centroids.
%
% multn: (n_minx1) output multiplicity.
%
% labn: (n_minx1) output classes.
%
% obslist: (n_minx1) output list of observations for the clustering file
%   system.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 02/Jul/13.
%
% Copyright (C) 2014  José Camacho Páez
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

if nargin < 2, error('Error in the number of arguments.'); end;
if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if nargin < 3, mult = ones(s(1),1); end;
if nargin < 4, lab = ones(s(1),1); end;
if nargin < 5, mat = eye(s(2)); end;
if nargin < 6, obslist = []; end;


% Main code

u = x*mat;
D = Inf*ones(s(1)); % initialization of the (upper triangular) matrix of Mahalanobis distances 
for i=1:s(1),
    for j=i+1:s(1),
        if lab(i)==lab(j), % if belong to the same class, compute the distance between observations i and j
            r = (u(i,:)-u(j,:))';
            D(i,j) = r'*r; 
        else % of different classes, not comparable
            D(i,j) = Inf;
        end
        
    end
end

centr = x;
multn = mult;
labn = lab;
for i=s(1)-1:-1:n_min, % reduction to n_min clusters
    
    % Computation of the minimum distance between observations or clusters
    min_dist = find(min(min(D))==D,1);
    
    % Location of the element with minimum distance
    row = mod(min_dist,i+1);
    if ~row, row = i+1; end
    column = ceil(min_dist/(i+1)); 
    if multn(column)>multn(row),
        aux_v = row;
        row = column;
        column = aux_v;
    end     
        
    % Actualization of the list
    if ~isempty(obslist),
        obslist{row} = [obslist{row} obslist{column}];
    end
    
    % Actualization of centroids and multiplicity 
    centr(row,:) = (multn(row)*centr(row,:)+multn(column)*centr(column,:))/sum(multn([row column])); % centroids 
    multn(row) = sum(multn([row column])); 
   
    % Actualization of the distance
    u(row,:) = centr(row,:)*mat;
    r = u-ones((i+1),1)*u(row,:);
    new_dist = sum(r.^2,2);
    labdiff = find(labn~=labn(row)); 
    new_dist(labdiff) = Inf; 
    D(1:(row-1),row) = new_dist(1:(row-1));
    D(row,(row+1):end) = new_dist((row+1):end);
   
    % Delete old cluster references
    if ~isempty(obslist),
        obslist = obslist([1:(column-1) (column+1):end]);
    end
    multn = multn([1:(column-1) (column+1):end]);
    labn = labn([1:(column-1) (column+1):end]); 
    centr = centr([1:(column-1) (column+1):end],:);
    u = u([1:(column-1) (column+1):end],:);
    D = D([1:(column-1) (column+1):end],[1:(column-1) (column+1):end]);
    
    if isempty(find(D<Inf)),
        indmin = i;
        break;
    end
        
end