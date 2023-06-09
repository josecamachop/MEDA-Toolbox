function ascao = asca(parglmo)

% ASCA is a data analysis algorithm for designed experiments. It does a 
% principal component analysis on the level averages of each experimental 
% factor in a designed experiment with balanced data. Interactions between 
% two factors can also be calculated. The original paper for this software 
% is Zwanenburg, G, Hoefsloot, HCJ, Westerhuis, JA, Jansen, JJ, Smilde, AK.
% ANOVA–principal component analysis and ANOVA–simultaneous component 
% analysis: a comparison. Journal of Chemometrics, 2018, 25:561-567.
%
% Related routines: parglm, paranova, apca, gasca, create_design
%
% ascao = asca(parglm)   % complete call
%
%
% INPUTS:
%
% parglmo (structure): structure with the factor and interaction matrices, 
% p-values and explained variance. Obtained with parallel anova (deprecated) 
% or parallel general linear model.
%
%
% OUTPUTS:
%
% ascao (structure): structure that contains scores, loadings, singular
% values and projections of the factors and interactions.
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two factors, with 4 and 3 levels, but only the first one 
%   is significative, and 4 replicates:
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1);
% end
%
% [table, parglmo] = parglm(X, F);
% table
% 
% ascao = asca(parglmo);
%
% for i=1:2, % Note, the second factor is shown for the sake of illustration
%   scores(ascao.factors{i},[],[],sprintf('Factor %d',i),[],ascao.design(:,i));
%   loadings(ascao.factors{i},[],sprintf('Factor %d',i));
% end
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two significative factors, with 4 and 3 levels, and 4 replicates:
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     fi{i} = randn(1,vars);
% end
% for j = 1:length(levels{2}),
%     fj{j} = randn(1,vars);
% end
% for i = 1:length(levels{1}),
%     for j = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,8) + repmat(fi{i} + fj{j},reps,1);
%     end
% end
%
% [table, parglmo] = parglm(X, F, {[1 2]});
% table
% 
% ascao = asca(parglmo);
%
% for i=1:2,
%   scores(ascao.factors{i},[],[],sprintf('Factor %d',i),[],ascao.design(:,i));
%   loadings(ascao.factors{i},[],sprintf('Factor %d',i));
% end
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two factors with 4 and 3 levels, and 4 replicates, with 
%   significant interaction:
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     for j = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,8) + repmat(randn(1,vars),reps,1);
%     end
% end
%
% [table, parglmo] = parglm(X, F, {[1 2]});
% table
% 
% ascao = asca(parglmo);
%
% M = ascao.factors{1}.matrix + ascao.factors{2}.matrix + ascao.interactions{1}.matrix;
% code_levels = F(:,1)*10+F(:,2);
% scores_pca(M,1:2,X,0,101,[],code_levels);
% legend(num2str(unique(code_levels)))
%
% loadings_pca(M,1:2,0);
%
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 19/May/23
%
% Copyright (C) 2023  University of Granada, Granada
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

ascao = parglmo;

%Do PCA on level averages for each factor
for factor = 1 : ascao.n_factors
    
    xf = ascao.factors{factor}.matrix;
    p = pca_pp(xf,1:rank(xf));
    
    ascao.factors{factor}.var = trace(xf'*xf);
    ascao.factors{factor}.lvs = 1:size(p,2);
    ascao.factors{factor}.loads = p;
    ascao.factors{factor}.scores = xf*p;
    ascao.factors{factor}.scoresV = (xf+ascao.residuals)*p;
end

%Do PCA on interactions
for interaction = 1 : ascao.n_interactions
    
    xf = ascao.interactions{interaction}.matrix;
    p = pca_pp(xf,1:rank(xf));
    
    ascao.interactions{interaction}.var = trace(xf'*xf);
    ascao.interactions{interaction}.lvs = 1:size(p,2);
    ascao.interactions{interaction}.loads = p;
    ascao.interactions{interaction}.scores = xf*p;
    ascao.interactions{interaction}.scoresV = (xf+ascao.residuals)*p;
end

ascao.type = 'ASCA';

