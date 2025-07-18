function gascao = gasca(paranovaost,c)

% GASCA is a data analysis algorithm for designed experiments. It does a
% group-wise principal component analysis on the level averages of each
% experimental factor in a designed experiment with balanced data.
% Interactions between two factors can also be calculated. The original
% paper is Saccenti, E., Smilde, A.K. and Camacho, J. Group-wise ANOVA
% simultaneous component analysis for designed omics experiments.
% Submitted to Metabolomics, 2018.
%
% gascao = gasca(paranovaost,c)   % complete call
%
% See also: parglm, paranova, asca, apca, createDesign 
%
%
% INPUTS:
%
% paranovaost (structure): structure with the factor and interaction
% matrices, p-values and explained variance. Obtained with parallel anova
% and where the field 'states' contains cells with the groups of variables 
% per factor and interaction.
%
%
% OUTPUTS:
%
% gascao (structure): structure that contains scores, loadings, singular
% values and projections of the factors and interactions.
%
%
% EXAMPLE OF USE: Random data, two significative factors, with 4 and 3 
%   levels, and 4 replicates, sparse relevant loadings:
%
% reps = 4;
% vars = 50;
% levels = {[1,2,3,4],[1,2,3]};
% int1 = 10:15;
% int2 = 30:37;
% 
% F = createDesign(levels,'Replicates',reps);
% 
% X = 0.1*randn(size(F,1),vars);
% for i = 1:length(levels{1}),
%   X(find(F(:,1) == levels{1}(i)),int1) = X(find(F(:,1) == levels{1}(i)),int1) + simuleMV(reps*length(levels{2}),length(int1),'LevelCorr',8) + repmat(randn(1,length(int1)),reps*length(levels{2}),1);
% end  
% for i = 1:length(levels{2}),
%   X(find(F(:,2) == levels{2}(i)),int2) = X(find(F(:,2) == levels{2}(i)),int2) + simuleMV(reps*length(levels{1}),length(int2),'LevelCorr',8) + repmat(randn(1,length(int2)),reps*length(levels{1}),1);
% end
% 
% [table, paranovaost] = parglm(X, F);
% 
% for i=1:length(paranovaost.factors)
%   map = corr(paranovaost.factors{i}.matrix);
%   plotMap(map);
%   c = input('Introduce threshold for correlation in interval (0,1): ');
%   [bel,paranovaost.factors{i}.states] = gia(map,'Gamma',c);
% end
% 
% gascao = gasca(paranovaost,c);
% 
% for i=1:2,
%   scores(gascao.factors{i},'Title',sprintf('Factor %d',i),'ObsClass',gascao.design(:,i));
%   loadings(gascao.factors{i},'Title',sprintf('Factor %d',i));
% end
%
% EXAMPLE OF USE: Same example with MEDA:
%
% reps = 4;
% vars = 50;
% levels = {[1,2,3,4],[1,2,3]};
% int1 = 10:15;
% int2 = 30:37;
% 
% F = createDesign(levels,'Replicates',reps);
% 
% X = 0.1*randn(size(F,1),vars);
% for i = 1:length(levels{1}),
%   X(find(F(:,1) == levels{1}(i)),int1) = X(find(F(:,1) == levels{1}(i)),int1) + simuleMV(reps*length(levels{2}),length(int1),'LevelCorr',8) + repmat(randn(1,length(int1)),reps*length(levels{2}),1);
% end  
% for i = 1:length(levels{2}),
%   X(find(F(:,2) == levels{2}(i)),int2) = X(find(F(:,2) == levels{2}(i)),int2) + simuleMV(reps*length(levels{1}),length(int2),'LevelCorr',8) + repmat(randn(1,length(int2)),reps*length(levels{1}),1);
% end
% 
% [table, paranovaost] = parglm(X, F);
% 
% for i=1:length(paranovaost.factors),
%   map = medaPca(paranovaost.factors{i}.matrix+paranovaost.residuals,'Preprocessing',0,'Threshold',0.3,'Option','100');
%   c = input('Introduce threshold for correlation in interval (0,1): ');
%   [bel,paranovaost.factors{i}.states] = gia(map,'Gamma',c);
% end
% 
% gascao = gasca(paranovaost,c);
% 
% for i=1:2,
%   scores(gascao.factors{i},'Title',sprintf('Factor %d',i),'ObsClass',gascao.design(:,i));
%   loadings(gascao.factors{i},'Title',sprintf('Factor %d',i));
% end
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 1/Jul/2025
% Dependencies: Matlab R2017b, MEDA v1.9
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


%% Main code

gascao = paranovaost;

%Do GPCA on level averages for each factor
for factor = 1 : gascao.nFactors
    
    xf = gascao.factors{factor}.matrix+gascao.residuals;
    map = medaPca(xf,'Preprocessing',0,'Seriated',true);
    
    gascao.factors{factor}.states = transformCrit(map,c(factor));
    
    p = gpca(xf,gascao.factors{factor}.states,'PCs',1:rank(xf));
    
    gascao.factors{factor}.var = trace(xf'*xf);
    gascao.factors{factor}.lvs = 1:size(p,2);
    gascao.factors{factor}.loads = p;
    gascao.factors{factor}.scores = xf*p;
    gascao.factors{factor}.scoresV = (xf+gascao.residuals)*p;
end

%Do GPCA on interactions
for interaction = 1 : gascao.nInteractions
    
    xf = gascao.interactions{interaction}.matrix+gascao.residuals;
    map = medaPca(xf,'Preprocessing',0,'Seriated',true);
    gascao.interactions{interaction}.states = transformCrit(map,c(length(gascao.factors)+interaction));
    
    p = gpca(xf,gascao.interactions{interaction}.states,'PCs',1:rank(xf));
    
    gascao.factors{factor}.var = trace(xf'*xf);
    gascao.interactions{interaction}.lvs = 1:size(p,2);
    gascao.interactions{interaction}.loads = p;
    gascao.interactions{interaction}.scores = xf*p;
    gascao.interactions{interaction}.scoresV = (xf+gascao.residuals)*p;
end

gascao.type = 'GASCA';

end

%% Auxiliary

function states = transformCrit(map,c)

lim = 1e-5;

if c<0
    [bel,states] = gia(map,'Gamma',-c);
else
    c2 = 0.99;
    [bel,states] = gia(map,'Gamma',c2);
    len = max(cellfun('length',states));
    if isempty(len), len=0; end
    while len~=c && c2 < 1-lim && c2 > lim
        diff = c - len;
        c2 = (-diff/(abs(diff)+10))*(1-c2) + c2;        
        if c2 < 1-lim && c2 > lim
            [bel,states] = gia(map,'Gamma',c2);
            len = max(cellfun('length',states));
            if isempty(len), len=0; end
        end
    end
    
end

end