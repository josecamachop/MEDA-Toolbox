function [p,t,bel,e] = gpca(xcs,states,varargin)

% Group-wise Principal Component Analysis. The original paper is Camacho, J., 
% Rodr�guez-G�mez, R., Saccenti, E. Group-wise Principal Component Analysis 
% for Exploratory Data Analysis. Journal of Computational and  Graphical 
% Statistics, 2017.
%
% p = gpca(xcs,states)     % minimum call
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% states: {Sx1} Cell with the groups of variables.
%
%
% Optional INPUTS (parameters):
%
% 'PCs': [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:rank(xcs)
%
% 'Tolerance': [1x1] tolerance value. By default, 1e-15.
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
% e: [NxM] matrix of residuals.
%
%
% EXAMPLE OF USE: Random data:
%
% x = simuleMV(20,10,'LevelCorr',8);
% pcs = 1:2;
% [map,~,ord] = medaPca(x,'PCs',pcs,'Seriated',true);
% x = x(:,ord);
% [bel,states] = gia(map,'Gamma',0.3);
% Xcs = preprocess2D(x,'Preprocessing',2);
% [p,t,bel] = gpca(Xcs,states,'PCs',pcs);
% for i=pcs,
%   plotVec(p(:,i),'XYLabel',{'',sprintf('Loadings PC %d',i)});
%   plotVec(t(:,i),'XYLabel',{'',sprintf('Scores PC %d',i)});
% end
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 03/Feb/2025
%
% Copyright (C) 2025  University of Granada, Granada
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
if length(states)==0
    states{1} = 1:M;
end


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
PCS = 0:rank(xcs);
Tol = 1e-15;
addParameter(p,'PCs',PCS);  
addParameter(p,'Tolerance',Tol);          
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
pcs = p.Results.PCs;
tol = p.Results.Tolerance;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
%pcs(find(pcs>rank(xcs))) = [];
A = length(pcs);

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''PCs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);


% Validate values of input data
assert (iscell(states), 'Value Error: parameter ''states'' must be a cell of positive integers. Type ''help %s'' for more info.', routine(1).name);
for i=1:length(states)
    assert (isempty(find(states{i}<1)) && isequal(fix(states{i}), states{i}), 'Value Error: 2nd argument must be a cell of positive integers. Type ''help %s'' for more info.', routine(1).name);
    assert (isempty(find(states{i}>M)), 'Value Error: parameter ''states'' must contain values not higher than M. Type ''help %s'' for more info.', routine(1).name);
end
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''PCs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);



%% Main code

map = xcs'*xcs;
I =  eye(size(map));
B = I;

p = [];
t = [];
bel = [];
for j = 1:max(pcs)  
    
    R = zeros(M,length(states));
    S = zeros(N,length(states));
    for i=1:length(states) % construct eigenvectors according to states
        mapaux = zeros(size(map));
        mapaux(states{i},states{i})= map(states{i},states{i});
       
        if find(mapaux>tol)
            [V,D] = eig(mapaux);
            ind = find(diag(D)==max(diag(D)),1);
            R(:,i) = V(:,ind)/sqrt((V(:,ind)'*B*V(:,ind)));
            S(:,i) = xcs*R(:,i);   
        end
    end

    sS = sum(S.^2,1); % select pseudo-eigenvector with the highest variance
    ind = find(sS==max(sS),1);
    
    q = B*R(:,ind); % deflate (Mackey'09)
    p(:,j) = R(:,ind)/norm(R(:,ind));
    t(:,j) = xcs*p(:,j);
    bel(j) = ind;
    map = (I-q*q')*map*(I-q*q');
    xcs = xcs*(I-q*q');
    B = B*(I-q*q');
    
end

% Postprocessing
p = p(:,pcs);
t = t(:,pcs);
bel = bel(pcs);

e = xcs;
