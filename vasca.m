function vascao = vasca(parglmoVS,siglev)

% Variable-selection ASCA is a data analysis algorithm for designed 
% experiments. Camacho J, Vitale R, Morales-Jiménez D, Gómez-Llorente C. 
% Variable-selection ANOVA Simultaneous Component Analysis (VASCA). 
% Bioinformatics. 2023 Jan 1;39(1):btac795. 
%
% Related routines: parglmVS, parglm, asca
%
% vascao = vasca(parglmVS,singlev)   % complete call
%
%
% INPUTS:
%
% parglmoVS (structure): structure with the factor and interaction matrices, 
% p-values and explained variance. Obtained with parallel general linear model
% with variable selection (parglmVS).
%
% siglev: [1x1] significance level (0.01 by default)
%
%
% OUTPUTS:
%
% vascao (structure): structure that contains scores, loadings, singular
% values and projections of the factors and interactions.
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Random data, one factor and two levels, three variables with information 
% on the factor.
%
% n_obs = 40;
% n_vars = 400;
% 
% class = (randn(n_obs,1)>0)+1;
% X = simuleMV(n_obs,n_vars,'LevelCorr',8);
% X(class==2,1:3) = X(class==2,1:3) + 10;
% 
% [TVS, parglmoVS] = parglmVS(X, class); % With variable selection
% TVS
% 
% vascao = vasca(parglmoVS);
% 
% if vascao.factors{1}.stasig
%    scores(vascao.factors{1},'Title',sprintf('Factor %d',1),'ObsClass',vascao.design(:,1));
%    loadings(vascao.factors{1},'Title',sprintf('Factor %d',1));
% end
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 23/Apr/2024
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if nargin < 2 || isempty(siglev), siglev = 0.01; end;

% Validate dimensions of input data
assert (isequal(size(siglev), [1 1]), 'Dimension Error: parameter ''singlev'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

%% Main code

vascao = parglmoVS;

%Do PCA on level averages for each factor
for factor = 1 : vascao.n_factors
    
    pvals = parglmoVS.p(parglmoVS.ord_factors(factor,:),factor); 
    M = find(pvals==min(pvals)); M = M(end);
    if ~isempty(pvals(M) <= siglev)
        vascao.factors{factor}.stasig = true;
        ind = parglmoVS.ord_factors(factor,1:M);
        xf = vascao.factors{factor}.matrix(:,ind);
        p = pca_pp(xf,'Pcs',1:rank(xf));
    
        vascao.factors{factor}.ind = ind;
        vascao.factors{factor}.var = trace(xf'*xf);
        vascao.factors{factor}.lvs = 1:size(p,2);
        vascao.factors{factor}.loads = p;
        vascao.factors{factor}.scores = xf*p;
        vascao.factors{factor}.scoresV = (xf+vascao.residuals(:,ind))*p;
    else
        vascao.factors{factor}.stasig = false;
    end
end

%Do PCA on interactions
for interaction = 1 : vascao.n_interactions
    
    pvals = parglmoVS.p(parglmoVS.ord_interactions(interaction,:),interaction+vascao.n_factors); 
    M = find(pvals==min(pvals)); M = M(end);
    if ~isempty(pvals(M) <= siglev)
        vascao.interactions{interaction}.stasig = true;
        ind = parglmoVS.ord_interactions(interaction,1:M);
        xf = vascao.interactions{interaction}.matrix(:,ind);
        for factor = 1 : vascao.interactions{1}.factors
            xf = xf + vascao.factors{factor}.matrix(:,ind);
        end
        p = pca_pp(xf,1:rank(xf));
    
        vascao.interactions{interaction}.ind = ind;
        vascao.interactions{interaction}.var = trace(xf'*xf);
        vascao.interactions{interaction}.lvs = 1:size(p,2);
        vascao.interactions{interaction}.loads = p;
        vascao.interactions{interaction}.scores = xf*p;
        vascao.interactions{interaction}.scoresV = (xf+vascao.residuals(:,ind))*p;
    else
        vascao.interactions{interaction}.stasig = false;
    end
end

vascao.type = 'VASCA';
