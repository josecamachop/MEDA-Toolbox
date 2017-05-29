function [P,T,bel,E,Lmodel] = Lgpca(Lmodel,states)

% Group-wise Principal Component Analysis for large data. The original 
% paper is Camacho, J., Rodríguez-Gómez, R., Saccenti, E. Group-wise 
% Principal Component Analysis for Exploratory Data Analysis. Journal of 
% Computational and  Graphical Statistics, 2017.
%
% [P,T,bel,E,Lmodel] = Lgpca(Lmodel,states)     % complete call
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
% P: [MxA] matrix of loadings.
%
% T: [NxA] matrix of scores.
%
% bel: [Ax1] correspondence between PCs and States.
%
% E: [NxM] matrix of residuals.
%
% Lmodel: (struct Lmodel) model after integrity checking.
%
%
% EXAMPLE OF USE: Random data:
%
% X = simuleMV(20,10,8);
% Lmodel = Lmodel_ini(X);
% Lmodel.lvs = 0:10;
% map = meda_Lpca(Lmodel);
% [bel,states] = gia(map,0.3);
% Lmodel.lvs = 1:length(states);
% [P,T,bel,E] = Lgpca(Lmodel,states);
% 
% for i=Lmodel.lvs,
%   plot_vec(P(:,i),[],[],{'',sprintf('PC %d',i)});
% end
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 21/May/2017
%
% Copyright (C) 2017  University of Granada, Granada
% Copyright (C) 2017  Jose Camacho Paez
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
[ok, Lmodel] = check_Lmodel(Lmodel);


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
for j = 1:max(Lmodel.lvs),  
    
    R = zeros(M,length(states));
    S = zeros(N,length(states));
    
    for i=1:length(states), % construct eigenvectors according to states
        map_aux = zeros(size(map));
        map_aux(states{i},states{i})= map(states{i},states{i});
        if rank(map_aux),
            [V,D] = eig(map_aux);
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

E = xcs;