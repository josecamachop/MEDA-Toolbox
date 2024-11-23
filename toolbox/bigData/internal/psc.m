function [centr,multn,classn,olabn,updatedn,obslist] = psc(x,nmin,mult,class,olab,updated,mat,obslist)

% Projected Sequential Clustering. 
%
% [centr,multn,classn] = psc(x,nmin)          % minimum call
% [centr,multn,classn,olabn,updatedn,obslist] = psc(x,nmin,mult,class,olab,updated,mat,obslist) % complete call
%
%
% INPUTS:
%
% x: [NxM] original matrix with centroids.
%
% nmin: [1x1] number of output clusters.
%
% mult: [Nx1] multiplicity of each cluster.
%
% class: [Nx1] class associated to each cluster.
%
% olab: {Nx1} label of each cluster.
%
% updated: [Nx1] specifies if the data is a new point
%   0: old point
%   1: new point
%
% mat: [MxA] projection matrix for distance computation.
%
% obslist: [Nx1] list of observations for the clustering file system.
%
%
% OUTPUTS:
%
% centr: [nminxM] output centroids.
%
% multn: [nminx1] output multiplicity.
%
% classn: [nminx1] output classes.
%
% olabn: {nminx1} output labels.
%
% updatedn: {nminx1} output updated values
%   0: old point
%   1: new point or combination with a new point
%
% obslist: [nminx1] output list of observations for the clustering file
%   system.
%
%
% EXAMPLE OF USE: Random values from 1000 to 20 clusters
%
% X = simuleMV(1000,2,'LevelCorr',8);
% [centr,multn] = psc(X,20);
% plotScatter(centr,[],[],{'Var 1', 'Var 2'},[],[],multn);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 20/Nov/2024
%
% Copyright (C) 2024  University of Granada, Granada
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

routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);

if nargin < 3 || isempty(mult), mult = ones(N,1); end;
if nargin < 4 || isempty(class), class = ones(N,1); end;
if nargin < 5, olab = {}; end;
if nargin < 6 || isempty(updated), updated = ones(N,1); end;
if nargin < 7 || isempty(mat), mat = eye(M); end;
if nargin < 8, obslist = {}; end;

% Validate dimensions of input data
assert (isequal(size(nmin), [1 1]), 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(mult), [N 1]), 'Dimension Error: 3rd argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(class), [N 1]), 'Dimension Error: 4th argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(updated), [N 1]), 'Dimension Error: 6th argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(mat,1), M), 'Dimension Error: 7th argument must be M-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (size(mat,2)<=M, 'Dimension Error: 7th argument must be M-by-A. Type ''help %s'' for more info.', routine(1).name);


% Validate values of input data
assert (nmin>0, 'Value Error: 2nd argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(mult<=0)), 'Value Error: 3rd argument must be above 0. Type ''help %s'' for more info.', routine(1).name);


%% Main code

u = x*mat;
D = Inf*ones(N); % initialization of the (upper triangular) matrix of Mahalanobis distances 
uc = unique(class);
for i=1:length(uc)
    ind = find(ismember(class,uc(i)));
    if length(ind) > 1
        indM = reshape(1:length(ind)^2,length(ind),length(ind));
        indM = triu(indM) - diag(diag(indM));
        indM = indM';
        indM(indM == 0) = [];
        D2 = Inf(length(ind));
        D2(indM) = pdist(u(ind,:),'squaredeuclidean');
        D(ind,ind) = D2;
    end
end


centr = x;
multn = mult;
classn = class;
olabn = olab;
updatedn = updated;
for i=N-1:-1:nmin % reduction to nmin clusters
    
    % Computation of the minimum distance between observations or clusters
    mindist = find(min(min(D))==D,1);
    
    % Location of the element with minimum distance
    row = mod(mindist,i+1);
    if ~row, row = i+1; end
    column = ceil(mindist/(i+1)); 
    if multn(column)>multn(row)
        auxv = row;
        row = column;
        column = auxv;
    end
    
    % Actualization of the list
    if ~isempty(obslist)
        obslist{row} = [obslist{row} obslist{column}];
    end
    
    % Actualization of labels
    if ~isempty(olabn)
        if ~isequal(olabn{row},olabn{column})
            olabn{row} = 'mixed';
        end      
    end
    
    % Actualization of centroids, multiplicity and updated values 
    centr(row,:) = (multn(row)*centr(row,:)+multn(column)*centr(column,:))/sum(multn([row column])); % centroids 
    multn(row) = sum(multn([row column])); 
    updatedn(row) = max(updatedn([row column])); 
   
    % Actualization of the distance
    u(row,:) = centr(row,:)*mat;
    r = u-ones((i+1),1)*u(row,:);
    newdist = sum(r.^2,2); 
    classdiff = find(~ismember(classn,classn(row))); 
    newdist(classdiff) = Inf; 
    D(1:(row-1),row) = newdist(1:(row-1));
    D(row,(row+1):end) = newdist((row+1):end);
   
    % Delete old cluster references
    if ~isempty(obslist)
        obslist = obslist([1:(column-1) (column+1):end]);
    end
    if ~isempty(olabn)
        olabn = olabn([1:(column-1) (column+1):end]);
    end
    multn = multn([1:(column-1) (column+1):end]);
    classn = classn([1:(column-1) (column+1):end]);
    updatedn = updatedn([1:(column-1) (column+1):end]); 
    centr = centr([1:(column-1) (column+1):end],:);
    u = u([1:(column-1) (column+1):end],:);
    D = D([1:(column-1) (column+1):end],[1:(column-1) (column+1):end]);
    
    if isempty(find(~isinf(D)))
        indmin = i;
        break;
    end
        
end