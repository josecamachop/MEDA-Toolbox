function Lmodel = Lgpca(Lmodel,states)

% Group-wise Principal Component Analysis for large data. The original 
% paper is Camacho, J., Rodríguez-Gómez, R., Saccenti, E. Group-wise 
% Principal Component Analysis for Exploratory Data Analysis. Journal of 
% Computational and  Graphical Statistics, 2017.
%
% Lmodel = Lgpca(Lmodel,states)     % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the GPCA
%   model:
%       Lmodel.XX: [MxM] X-block cross-product matrix.
%       Lmodel.lvs: [1x1] number of PCs A.
%
% states: {Sx1} Cell with the groups of variables.
%
%
% OUTPUTS:
%
% Lmodel: (struct Lmodel) model after integrity checking.
%
%
% EXAMPLE OF USE: Random data:
%
% X = simuleMV(20,10,'LevelCorr',6);
% Lmodel = iniLmodel(X);
% Lmodel.lvs = 0:10;
% map = medaLpca(Lmodel);
% [bel,states] = gia(map,'Gamma',0.3);
% Lmodel.lvs = 1:length(states);
% Lmodel = Lgpca(Lmodel,states);
% 
% for i=Lmodel.lvs,
%   plotVec(Lmodel.loads(:,i),'XYLabel',{'',sprintf('PC %d',i)});
% end
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 19/Nov/2024
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

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
[ok, Lmodel] = checkLmodel(Lmodel);


%% Main code

map = Lmodel.XX;
xcs = Lmodel.centr;
N = size(xcs, 1);
M = size(xcs, 2);
I =  eye(size(map));
B = I;

P = [];
T = [];
bel = [];
for j = 1:max(Lmodel.lvs)  
    
    R = zeros(M,length(states));
    S = zeros(N,length(states));
    
    for i=1:length(states) % construct eigenvectors according to states
        mapaux = zeros(size(map));
        mapaux(states{i},states{i})= map(states{i},states{i});
        if rank(mapaux)
            [V,D] = eig(mapaux);
            ind = find(diag(D)==max(diag(D)),1);
            R(:,i) = V(:,ind);
            S(:,i) = xcs*R(:,i);    
        end
    end

    sS = sum(S.^2,1); % select pseudo-eigenvector with the highest variance
    ind = find(sS==max(sS),1);
    P(:,j) = R(:,ind);
    T(:,j) = S(:,ind);
    bel(j) = ind;
    
    q = B*R(:,ind); % deflate (Mackey'09)
    map = (I-q*q')*map*(I-q*q');
    xcs = xcs*(I-q*q');
    B = B*(I-q*q');
    
end

Lmodel.loads = P;
Lmodel.scores = T;
Lmodel.bel = bel;
Lmodel.residuals = xcs;