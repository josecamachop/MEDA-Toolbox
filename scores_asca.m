function T = scores_asca(paranovao, alpha, labels)

% ASCA is a data analysis algorithm for designed experiments. It does a 
% principal component analysis on the level averages of each experimental 
% factor in a designed experiment with balanced data. Interactions between 
% two factors can also be calculated. The original paper for this software 
% is Zwanenburg, G, Hoefsloot, HCJ, Westerhuis, JA, Jansen, JJ, Smilde, AK.
% ANOVA–principal component analysis and ANOVA–simultaneous component 
% analysis: a comparison. Journal of Chemometrics, 2018, 25:561-567.
%
% T = scores_asca(paranovao)   % minimum call
% T = scores_asca(paranovao, alpha, labels)   % complete call
%
%
% INPUTS:
%
% paranovao (structure): structure with the factor and interaction
% matrices, p-values and explained variance. Obtained with parallel anova
% (paranova)
%
% alpha: [1x1] significance level (0.05 by default)
%
% labels: [Kx1] name of the factors and interactions (empty by default)
%
%
% OUTPUTS:
%
% T: [NxA] calibration scores
%
%
% EXAMPLE OF USE: Random data, two significative factors, with 4 and 3 
%   levels, and 4 replicates:
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
% paranovao = paranova(X, F);
%
% T = scores_asca(paranovao,0.05,{'Factor 1','Factor 2'});
%
%
% EXAMPLE OF USE: Random data, two factors, with 4 and 3 levels, but only
%   the first one is significative, and 4 replicates:
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
% paranovao = paranova(X, F);
%
% T = scores_asca(paranovao,0.05,{'Factor 1','Factor 2'});
%
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 21/Mar/18
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
if nargin < 2 || isempty(alpha), alpha = 0.05; end;
if nargin < 3, labels = []; end
    
% Validate dimensions of input data
assert (isequal(size(alpha), [1 1]), 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

T = {};
ind = find(paranovao.p<alpha);
for i=1:length(ind)
    xf = paranovao.factors.means{ind(i)};
    e = paranovao.residuals;
    T{ind(i)} = scores_pca(xf,1:rank(xf),xf+e,0,101,[],paranovao.design(:,ind(i)));
    if ~isempty(labels),
        title(labels{ind(i)})
    end
end
