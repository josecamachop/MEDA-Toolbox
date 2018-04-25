function gascao = gasca(paranovao_st)

% GASCA is a data analysis algorithm for designed experiments. It does a
% group-wise principal component analysis on the level averages of each
% experimental factor in a designed experiment with balanced data.
% Interactions between two factors can also be calculated. The original
% paper is Saccenti, E., Smilde, A.K. and Camacho, J. Group-wise ANOVA
% simultaneous component analysis for designed omics experiments.
% Submitted to Metabolomics, 2018.
%
% ggascao = gasca(paranovao_st)   % complete call
%
%
% INPUTS:
%
% paranovao_st (structure): structure with the factor and interaction
% matrices, p-values and explained variance. Obtained with parallel anova
% (paranovao_st) and where the field 'states' contains cells with the groups
% of variables per factor and interaction.
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
% F = create_design(levels,reps);
%
% X = 0.1*randn(size(F,1),vars);
% for i = 1:length(levels{1}),
%   X(find(F(:,1) == levels{1}(i)),int1) = X(find(F(:,1) == levels{1}(i)),int1) + simuleMV(reps*length(levels{2}),length(int1),8) + repmat(randn(1,length(int1)),reps*length(levels{2}),1);
% end  
% for i = 1:length(levels{2}),
%   X(find(F(:,2) == levels{2}(i)),int2) = X(find(F(:,2) == levels{2}(i)),int2) + simuleMV(reps*length(levels{1}),length(int2),8) + repmat(randn(1,length(int2)),reps*length(levels{1}),1);
% end
%
% paranovao_st = paranova(X, F);
% 
% for i=1:length(paranovao.factors),
%   map = corr(paranovao_st.factors{i}.means);
%   plot_map(map);
%   c = input('Introduce threshold for correlation in interval (0,1): ');
%   [bel,paranovao_st.factors{i}.states] = gia(map,c);
% end
%         
% gascao = gasca(paranovao_st);
%
% for i=1:2,
%   scores(gascao.factors{i},[],[],sprintf('Factor %d',i),[],gascao.design(:,i));
%   loadings(gascao.factors{i},11,sprintf('Factor %d',i));
% end
%
%
% EXAMPLE OF USE: Same example with MEDA:
%
% reps = 4;
% vars = 50;
% levels = {[1,2,3,4],[1,2,3]};
% int1 = 10:15;
% int2 = 30:37;
%
% F = create_design(levels,reps);
%
% X = 0.1*randn(size(F,1),vars);
% for i = 1:length(levels{1}),
%   X(find(F(:,1) == levels{1}(i)),int1) = X(find(F(:,1) == levels{1}(i)),int1) + simuleMV(reps*length(levels{2}),length(int1),8) + repmat(randn(1,length(int1)),reps*length(levels{2}),1);
% end  
% for i = 1:length(levels{2}),
%   X(find(F(:,2) == levels{2}(i)),int2) = X(find(F(:,2) == levels{2}(i)),int2) + simuleMV(reps*length(levels{1}),length(int2),8) + repmat(randn(1,length(int2)),reps*length(levels{1}),1);
% end
%
% paranovao_st = paranova(X, F);
% 
% for i=1:length(paranovao.factors),
%   map = meda_pca(paranovao_st.factors{i}.means+paranovao_st.residuals,[],0,0.3,'100');
%   c = input('Introduce threshold for correlation in interval (0,1): ');
%   [bel,paranovao_st.factors{i}.states] = gia(map,c);
% end
%         
% gascao = gasca(paranovao_st);
%
% for i=1:2,
%   scores(gascao.factors{i},[],[],sprintf('Factor %d',i),[],gascao.design(:,i));
%   loadings(gascao.factors{i},11,sprintf('Factor %d',i));
% end
%
%
% Related routines: paranova, asca, apca, create_design 
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 25/Apr/18
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez
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

gascao = paranovao_st;

%Do GPCA on level averages for each factor
for factor = 1 : gascao.n_factors
    
    xf = gascao.factors{factor}.means;
    p = gpca(xf,gascao.factors{factor}.states,1:length(gascao.factors{factor}.states));
    
    gascao.factors{factor}.lvs = 1:size(p,2);
    gascao.factors{factor}.loads = p;
    gascao.factors{factor}.scores = (xf+gascao.residuals)*p;
end

%Do GPCA on interactions
for interaction = 1 : gascao.n_interactions
    
    xf = gascao.interactions{interaction}.means;
    p = gpca(xf,gascao.interactions{interaction}.states,1:length(gascao.interactions{interaction}.states));
    
    gascao.interactions{interaction}.lvs = 1:size(p,2);
    gascao.interactions{interaction}.loads = p;
    gascao.interactions{interaction}.scores = (xf+gascao.residuals)*p;
end

