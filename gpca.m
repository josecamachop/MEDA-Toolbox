function [p,t,bel] = gpca(xcs,states,pcs)

% Group-wise Principal Component Analysis.
%
% p = gpca(xcs,states)     % minimum call
% [p,t,bel] = gpca(xcs,states,pcs)     % complete call
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% states: {Sx1} Cell with the groups of variables.
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:rank(xcs)
%
%
% OUTPUTS:
%
% p: [MxA] matrix of loadings.
%
% t: [NxA] matrix of scores.
%
% bel: [Ax1] correspondence between PCs and States.
%
%
% EXAMPLE OF USE: Random data:
%
% X = real(ADICOV(randn(10,10).^19,randn(100,10),10));
% pcs = 1:rank(X);
% map = meda_pca(X,pcs);
% [bel,states] = gia(map,0.3);
% Xcs = preprocess2D(X,2);
% pcs = 1:length(states);
% [p,t,bel] = gpca(Xcs,states,pcs);
% 
% f = figure;
% for i=pcs,
%   subplot(max(pcs),1,i); 
%   plot_vec(p(:,i),[],[],sprintf('PC %d',i),[],[],f);
% end
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 05/Apr/16.
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(xcs, 1);
M = size(xcs, 2);
if nargin < 3 || isempty(pcs), pcs = 0:rank(xcs); end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
pcs(find(pcs>rank(xcs))) = [];
A = length(pcs);

% Validate dimensions of input data
assert (length(states)>0, 'Dimension Error: 2nd argument must contain at least one state. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);


% Validate values of input data
assert (iscell(states), 'Value Error: 2nd argument must be a cell of positive integers. Type ''help %s'' for more info.', routine(1).name);
for i=1:length(states),
    assert (isempty(find(states{i}<1)) && isequal(fix(states{i}), states{i}), 'Value Error: 2nd argument must be a cell of positive integers. Type ''help %s'' for more info.', routine(1).name);
    assert (isempty(find(states{i}>M)), 'Value Error: 2nd argument must contain values not higher than M. Type ''help %s'' for more info.', routine(1).name);
end
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: 3rd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);



%% Main code

map = xcs'*xcs;
I =  eye(size(map));
B = I;
for j = 1:max(pcs),  
    
    for i=1:length(states), % construct eigenvectors according to states
        map_aux = zeros(size(map));
        map_aux(states{i},states{i})= map(states{i},states{i});
        [V,D] = eig(map_aux);
        ind = find(diag(D)==max(diag(D)),1);
        R(:,i) = V(:,ind);
        S(:,i) = xcs*R(:,i);       
    end

    sS = sum(S.^2,1); % select pseudo-eigenvector with the highest variance
    ind = find(sS==max(sS),1);
    p(:,j) = R(:,ind);
    t(:,j) = S(:,ind);
    bel(j) = ind;
    
    q = B*R(:,ind); % deflate (Mackey'09)
    map = (I-q*q')*map*(I-q*q');
    xcs = xcs*(I-q*q');
    B = B*(I-q*q');
    
end

% Postprocessing
p = p(:,pcs);
t = t(:,pcs);
bel = bel(pcs);
