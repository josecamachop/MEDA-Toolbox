function apcao = apca(parglmo)

% ANOVA-PCA (APCA) is a data analysis algorithm for the analysis of designed 
% experiments. The input is a General Linear Models (GLM) factorization of
% the data (done with parglm and stored in parglmo) and the code applies 
% Principal Component Analysis to the factor/interaction matrices plus the
% residuals. 
%
% apcao = apca(parglmo)   % minimum call
%
% See also: parglm, asca, gasca, createDesign
%
%
% INPUTS:
%
% parglmo (structure): structure with the GLM decomposition with factor and 
% interaction matrices, p-values and explained variance. 
%
%
% OUTPUTS:
%
% apcao (structure): structure that contains scores, loadings, singular
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
% F = createDesign(levels,'Replicates',reps);
% 
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,'LevelCorr',8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1);
% end
% 
% [table, parglmo] = parglm(X, F);
% table
% 
% apcao = apca(parglmo);
% 
% for i=1:2,
%   scores(apcao.factors{i},'Title',sprintf('Factor %d',i),'ObsClass',apcao.design(:,i));
%   loadings(apcao.factors{i},'Title',sprintf('Factor %d',i));
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
% F = createDesign(levels,'Replicates',reps);
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
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,'LevelCorr',8) + repmat(fi{i} + fj{j},reps,1);
%     end
% end
% 
% [table, parglmo] = parglm(X, F, 'Model',[1 2]);
% table
% 
% apcao = apca(parglmo);
% 
% for i=1:2,
%   scores(apcao.factors{i},'Title',sprintf('Factor %d',i),'ObsClass',apcao.design(:,i));
%   loadings(apcao.factors{i},'Title',sprintf('Factor %d',i));
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
% F = createDesign(levels,'Replicates',reps);
% 
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     for j = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,'LevelCorr',8) + repmat(randn(1,vars),reps,1);
%     end
% end
% 
% [table, parglmo] = parglm(X, F, 'Model',[1 2]);
% table
% 
% apcao = apca(parglmo);
% 
% M = apcao.factors{1}.matrix + apcao.factors{2}.matrix + apcao.interactions{1}.matrix;
% codeLevels = {};
% for i=1:size(F,1), codeLevels{i} = sprintf('F1:%d,F2:%d',F(i,1),F(i,2));end;
% scoresPca(M,'PCs',1:2,'ObsTest',X,'Preprocessing',0,'PlotCal',false,'ObsClass',codeLevels);
% legend(unique(codeLevels))
% 
% loadingsPca(M,'PCs',1:2,'Preprocessing',0);
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);


%% Main code

apcao = parglmo;

%Do PCA on level averages for each factor
for factor = 1 : apcao.nFactors
    
    xf = apcao.factors{factor}.matrix+apcao.residuals;
    model = pcaEig(xf,'PCs',1:rank(apcao.factors{factor}.matrix));
    
    fnames = fieldnames(model);
    for n = 1:length(fnames)
        apcao.factors{factor} = setfield(apcao.factors{factor},fnames{n},getfield(model,fnames{n}));
    end
end

%Do PCA on interactions
for interaction = 1 : apcao.nInteractions
    
    xf = apcao.interactions{interaction}.matrix+apcao.residuals;
    model = pcaEig(xf,'PCs',1:rank(apcao.interactions{interaction}.matrix));

    fnames = fieldnames(model);
    for n = 1:length(fnames)
        apcao.interactions{interaction} = setfield(apcao.interactions{interaction},fnames{n},getfield(model,fnames{n}));
    end
end

apcao.type = 'APCA'
