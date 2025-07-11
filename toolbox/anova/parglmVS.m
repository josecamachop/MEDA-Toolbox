function [T, parglmo] = parglmVS(X, F, varargin)

% Parallel General Linear Model with variable-selection to factorize in 
% multivariate factor and interaction matrices in an experimental design 
% and permutation test for multivariate statistical significance. These 
% represent the two first steps in an variable-selection ANOVA Simultaneous 
% Component Analsysis (vASCA) 
%
% T = parglmVS(X, F)   % minimum call
%
% See also: parglm, parglmMC, vasca, createDesign
%
%
% INPUTS:
%
% X: [NxM] billinear data set for model fitting, where each row is a
% measurement, each column a variable
%
% F: [NxF] design matrix, cell or array, where columns correspond to 
% factors and rows to levels.
%
%
% Optional INPUTS (parameters):
%
% 'Model': This paremeter is similar to 'model' of anovan. It could be:
%       'linear': only main effects are provided (by default)
%       'interaction': two order interactions are provided
%       'full': all potential interactions are provided
%       [1x1]: maximum order of interactions considered
%       [ix2]: array with two order interactions
%       cell: with each element a vector of factors
%
% 'Preprocessing': [1x1] preprocesing:
%       0: no preprocessing 
%       1: mean-centering 
%       2: auto-scaling (default)  
%
% 'Permutations': [1x1] number of permutations (1000 by default)
%
% 'Ts': [1x1] Use SSQ (0) or the F-value (otherwise, by default) as test statistic  
%       0: Sum-of-squares of the factor/interaction
%       1: F-ratio of the SS of the factor/interaction divided by the SS of 
%       the residuals
%       2: F-ratio following the factors/interactions hierarchy (see
%       Montogomery, by default)
%
% 'Ordinal': [1xF] whether factors are nominal or ordinal
%       0: nominal (default)
%       1: ordinal
%
% 'Random': [1xF] whether factors are fixed or random
%       0: fixed (default)
%       1: random
%
% 'Mtc': [1x1] Multiple test correction
%       1: Bonferroni 
%       2: Holm step-up or Hochberg step-down
%       3: Benjamini-Hochberg step-down (FDR, by default)
%       -pi0: Q-value assuming 0<pi0<=1 the ratio of null variables   
% 
% 'Fmtc': [1x1] correct for multiple-tesis when multifactorial (multi-way)
% analysis.
%       0: do not correct (by default)
%       1: same as mtc 
%
% 'Coding': [1xF] type of coding of factors
%       0: sum/deviation coding (default)
%       1: reference coding (reference is the last level)
%
% 'Nested': [nx2] pairs of neted factors, e.g., if factor 2 is nested in 1,
%   and 3 in 2, then nested = [1 2; 2 3]
%
%
% OUTPUTS:
%
% T (table): ANOVA-like output table
%
% parglmoMV (structure): structure with the factor and interaction
% matrices, univariate p-values (corrected, depending on mtc and on fmtc), 
% multivariate p-values (corrected, depending on fmtc) and explained variance 
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Random data, one factor and two levels, three variables with information 
% on the factor, high background correlation. 
%
% nObs = 40;
% nVars = 100;
% 
% class = (randn(nObs,1)>0)+1;
% X = simuleMV(nObs,nVars,'LevelCorr',8);
% X(class==2,1:3) = X(class==2,1:3) + 10;
% 
% S = rng; % Use same seed for random generators to improve comparability of results
% [T, parglmo] = parglm(X,class); % No variables selection 
% rng(S);
% [TVS, parglmoVS] = parglmVS(X, class); % With variable selection
% 
% h = figure; hold on
% plot([1 nVars],-log10([parglmo.p parglmo.p]),'b-.')
% plot(-log10(parglmoVS.p(parglmoVS.ordFactors)),'g-o')
% plot([0,size(X,2)],-log10([0.01 0.01]),'r--')
% legend('ASCA','VASCA','alpha=0.01','Location','southeast')
% a=get(h,'CurrentAxes');
% set(a,'FontSize',14)
% ylabel('-log10(p-values)','FontSize',18)
% xlabel('Responses in selected order','FontSize',18)
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Random data, one factor and two levels, no variable with information 
% on the factor. 
%
% nObs = 40;
% nVars = 100;
% 
% class = (randn(nObs,1)>0)+1;
% X = simuleMV(nObs,nVars,'LevelCorr',8);
% 
% S = rng; % Use same seed for random generators to improve comparability of results
% [T, parglmo] = parglm(X,class); % No variables selection 
% rng(S);
% [TVS, parglmoVS] = parglmVS(X, class); % With multivariate variable selection
% 
% h = figure; hold on
% plot([1 nVars],-log10([parglmo.p parglmo.p]),'b-.')
% plot(-log10(parglmoVS.p(parglmoVS.ordFactors)),'g-o')
% plot([0,size(X,2)],-log10([0.01 0.01]),'r--')
% legend('ASCA','VASCA','alpha=0.01','Location','southeast')
% a=get(h,'CurrentAxes');
% set(a,'FontSize',14)
% ylabel('-log10(p-values)','FontSize',18)
% xlabel('Responses in selected order','FontSize',18)
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two factors with 4 and 3 levels, and 4 replicates, with 
%   significant interaction but only in the first 3 variables
%
% reps = 4;
% vars = 3;
% levels = {[1,2,3,4],[1,2,3]};
% 
% F = createDesign(levels,'Replicates',reps);
% 
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1})
%     for j = 1:length(levels{2})
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,'LevelCorr',8) + repmat(randn(1,vars),reps,1);
%     end
% end
% 
% X = [X simuleMV(length(F),397,'LevelCorr',8)];
% 
% table = parglm(X, F, 'Model',{[1 2]}, 'Random', [1 1])
% 
% [tableVs, parglmoVs] = parglmVS(X, F, 'Model',{[1 2]}, 'Random', [1 1]);
% tableVs
% 
% parglmoVs.p(1:10,:) % significance is only at the first variables
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 5/Apr/2025
% Dependencies: Matlab R2017b, MEDA v1.8
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
N = size(X, 1);
M = size(X, 2);

nFactors = size(F,2);                 % number of factors

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Model','linear'); 
addParameter(p,'Preprocessing',2);
addParameter(p,'Permutations',1000);  
addParameter(p,'Ts',2); 
addParameter(p,'Ordinal',zeros(1,size(F,2))); 
addParameter(p,'Random',zeros(1,size(F,2))); 
addParameter(p,'Mtc',3); 
addParameter(p,'Fmtc',0); 
addParameter(p,'Coding',zeros(1,size(F,2))); 
addParameter(p,'Nested',[]); 
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
model = p.Results.Model;
prep = p.Results.Preprocessing;
nPerm = p.Results.Permutations;
ts = p.Results.Ts;
ordinal = p.Results.Ordinal;
random = p.Results.Random;
mtc = p.Results.Mtc;
fmtc = p.Results.Fmtc;
coding = p.Results.Coding;
nested = p.Results.Nested;

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(nPerm), [1 1]), 'Dimension Error: parameter ''Permutations'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ts), [1 1]), 'Dimension Error: parameter ''Ts'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ordinal), [1 size(F,2)]), 'Dimension Error: parameter ''Ordinal'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(random), [1 size(F,2)]), 'Dimension Error: parameter ''Random'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(mtc), [1 1]), 'Dimension Error: parameter ''Mtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(fmtc), [1 1]), 'Dimension Error: parameter ''Fmtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(coding), [1 size(F,2)]), 'Dimension Error: parameter ''Coding'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);


%% Univariate Inference

if ts
    [~, parglmo, tsFactors, tsInteractions, SSQX, SSQinter, SSQFactors, SSQInteractions, SSQresiduals] = parglmMC(X, F, 'Permutations', ceil(nPerm/M), 'Model', model, 'Preprocessing', prep, 'Ts', ts, 'Ordinal', ordinal, 'Random', random, 'Mtc', mtc, 'Fmtc', fmtc, 'Coding', coding, 'Nested', nested);
else
    [~, parglmo, tsFactors, tsInteractions] = parglmMC(X, F, 'Permutations', ceil(nPerm/M), 'Model', model, 'Preprocessing', prep, 'Ts', ts, 'Ordinal', ordinal, 'Random', random, 'Mtc', mtc, 'Fmtc', fmtc, 'Coding', coding, 'Nested', nested);
end

parglmo.univp = parglmo.p;

if fmtc
    mtcc = parglmo.nFactors + parglmo.nInteractions;        % correction for the number of tests
else
    mtcc = 1;
end


%% Univariate inference sorts in terms of the p-value, but in vASCA we sort according to the statistic


for f = 1 : parglmo.nFactors
    [~,parglmo.ordFactors(f,:)] = sort(tsFactors(1,f,:),'descend');
end

for i = 1 : parglmo.nInteractions
    [~,parglmo.ordInteractions(i,:)] = sort(tsInteractions(1,i,:),'descend');
end


%% Incremental ASCAs

for var = 1 : M 

    for f = 1 : nFactors
        parglmo.factors{f}.refF = []; 
        parglmo.factors{f}.refI = [];

        ord = parglmo.ordFactors(f,:);
        MSq = sum(SSQresiduals(1,ord(1:var)))/parglmo.Rdf;
        if ts==2 % reference MSE
            SSref = 0;
            Dfref = 0;
            for f2 = 1 : parglmo.nFactors % Nested
                if ~isempty(find(f==parglmo.factors{f2}.factors & random(f2) == 1)) && (MSq < sum(SSQFactors(1,f2,ord(1:var)))/parglmo.df(f2)) % when a nested factor is random and significantly larger than the background noise
                    SSref = SSref + sum(SSQFactors(1,f2,ord(1:var)));
                    Dfref = Dfref + parglmo.df(f2);
                    parglmo.factors{f}.refF = [parglmo.factors{f}.refF f2];
                end
            end
            for i = 1 : parglmo.nInteractions
                if ~isempty(find(f==parglmo.interactions{i}.factors))
                    rest = setdiff(parglmo.interactions{i}.factors,f);
                    if (sum(random(rest) == 1) > 0) && (MSq < sum(SSQInteractions(1,i,ord(1:var)))/parglmo.dfint(i)) % when an interaction is random and larger than the background noise
                        SSref = SSref + sum(SSQInteractions(1,i,ord(1:var)));
                        Dfref = Dfref + parglmo.dfint(i);
                        parglmo.factors{f}.refI = [parglmo.factors{f}.refI i];
                    end
                end
            end
            if SSref == 0
                FFactors(1,f,ord(var)) = (sum(SSQFactors(1,f,ord(1:var)))/parglmo.df(f))/MSq;
            else
                FFactors(1,f,ord(var)) = (sum(SSQFactors(1,f,ord(1:var)))/parglmo.df(f))/(SSref/Dfref);
            end
        else
            FFactors(1,f,ord(var)) = (sum(SSQFactors(1,f,ord(1:var)))/parglmo.df(f))/MSq;
        end
    end
 
    for i = 1 : parglmo.nInteractions
        parglmo.interactions{i}.refI = [];

        ord = parglmo.ordInteractions(i,:);
        MSq = sum(SSQresiduals(1,ord(1:var)))/parglmo.Rdf;
        if ts==2 % reference MSE
            SSref = 0;
            Dfref = 0;
            for i2 = 1 : parglmo.nInteractions
                [~,ia,ib] = intersect(parglmo.interactions{i}.factors,parglmo.interactions{i2}.factors);
                if (length(ia) == length(parglmo.interactions{i}.factors) & length(parglmo.interactions{i2}.factors) > length(parglmo.interactions{i}.factors))
                    rest = setdiff(parglmo.interactions{i2}.factors,parglmo.interactions{i}.factors);
                    if (sum(random(parglmo.interactions{i2}.factors(rest)) == 1) > 0) && (MSq < sum(SSQInteractions(1,i2,ord(1:var)))/parglmo.dfint(i2)) % when the higher order interaction is random and larger than the background noise
                        SSref = SSref + sum(SSQInteractions(1,i2,ord(1:var)));
                        Dfref = Dfref + parglmo.dfint(i2);
                        parglmo.interactions{i}.refI = [parglmo.interactions{i}.refI i2];
                    end
                end
            end
            if SSref == 0
                FInteractions(1,i,ord(var)) = (sum(SSQInteractions(1,i,ord(1:var)))/parglmo.dfint(i))/MSq;
            else
                FInteractions(1,i,ord(var)) = (sum(SSQInteractions(1,i,ord(1:var)))/parglmo.dfint(i))/(SSref/Dfref);
            end
        else
            FInteractions(1,i,ord(var)) = (sum(SSQInteractions(1,i,ord(1:var)))/parglmo.dfint(i))/MSq;
        end
    end
    
    % Permutations
    parfor j = 1 : nPerm*mtcc

        Ff = zeros(1,parglmo.nFactors);
        for f = 1 : parglmo.nFactors
            [~,ord] = sort(tsFactors(1+j,f,:),'descend');
            if ts==2 % reference MSE
                SSref = 0;
                Dfref = 0;
                for f2 = 1 : parglmo.nFactors
                    if ~isempty(find(f2==parglmo.factors{f}.refF))
                        SSref = SSref + sum(SSQFactors(1+j,f2,ord(1:var)));
                        Dfref = Dfref + parglmo.df(f2);
                    end
                end

                for i = 1 : parglmo.nInteractions
                    if ~isempty(find(i==parglmo.factors{f}.refI))
                        SSref = SSref + sum(SSQInteractions(1+j,i,ord(1:var)));
                        Dfref = Dfref + parglmo.dfint(i);
                    end
                end
                if SSref == 0
                    Ff(f) =  (sum(SSQFactors(1+j,f,ord(1:var)))/parglmo.df(f))/(sum(SSQresiduals(1+j,ord(1:var)))/parglmo.Rdf);
                else
                    Ff(f) = (sum(SSQFactors(1+j,f,ord(1:var)))/parglmo.df(f))/(SSref/Dfref);
                end
            else
                Ff(f) = (sum(SSQFactors(1+j,f,ord(1:var)))/parglmo.df(f))/(sum(SSQresiduals(1+j,ord(1:var)))/parglmo.Rdf);
            end
        end
        FFactors(1+j,:,var) = Ff;

        Fi = zeros(1,parglmo.nInteractions);
        for i = 1 : parglmo.nInteractions
            [~,ord] = sort(tsInteractions(1+j,i,:),'descend');
            if ts==2 % reference MSE
                SSref = 0;
                Dfref = 0;
                for i2 = 1 : parglmo.nInteractions
                    if ~isempty(find(i2==parglmo.interactions{i}.refI))
                        SSref = SSref + sum(SSQInteractions(1+j,i2,ord(1:var)));
                        Dfref = Dfref + parglmo.dfint(i2);
                    end
                end
                if SSref == 0
                    Fi(i)  = (sum(SSQInteractions(1+j,i,ord(1:var)))/parglmo.dfint(i))/(sum(SSQresiduals(1+j,ord(1:var)))/parglmo.Rdf);
                else
                    Fi(i) = (sum(SSQInteractions(1+j,i,ord(1:var)))/parglmo.dfint(i))/(SSref/Dfref);
                end
            else
                Fi(i) = (sum(SSQInteractions(1+j,i,ord(1:var)))/parglmo.dfint(i))/(sum(SSQresiduals(1+j,ord(1:var)))/parglmo.Rdf);
            end
        end
        FInteractions(1+j,:,var) = Fi;

    end

end

% reorder SSQ permutes to faciliate the computation of p-values
% Permutations
if ~ts
    for j = 1 : nPerm*mtcc
        for f = 1 : parglmo.nFactors
            [~,ord] = sort(tsFactors(1 + j,f,:),'descend');
            SSQFactors(1 + j,f,:) = SSQFactors(1 + j,f,ord);
        end
        for i = 1 : parglmo.nInteractions
            [~,ord] = sort(tsInteractions(1 + j,i,:),'descend');
            SSQInteractions(1 + j,i,:) = SSQInteractions(1 + j,i,ord);
        end
    end
end


% Calculate multivariate p-values
pFactor            = zeros(parglmo.nFactors,M);                % p-values factors
pInteraction       = zeros(parglmo.nInteractions,M);           % p-values interactions

for f = 1 : parglmo.nFactors
    for var = 1 : M
        if ts
            pFactor(f,parglmo.ordFactors(f,var)) = (size(find( FFactors(2:(nPerm*mtcc + 1),f,var) ...
                >= FFactors(1, f,parglmo.ordFactors(f,var))), 1) + 1)/(nPerm*mtcc+1);
        else
            pFactor(f,parglmo.ordFactors(f,var)) = (size(find( sum(SSQFactors(2:(nPerm*mtcc + 1),f,1:var),3) ...
                >= sum(SSQFactors(1, f,parglmo.ordFactors(f,1:var)))), 1) + 1)/(nPerm*mtcc+1);
        end
    end
end
for i = 1 : parglmo.nInteractions
    for var = 1 : M
        if ts
            pInteraction(i,parglmo.ordInteractions(i,var)) = (size(find( FInteractions(2:(nPerm*mtcc + 1),i,var) ...
                >= FInteractions(1, i, parglmo.ordInteractions(i,var))), 1) + 1)/(nPerm*mtcc+1);
        else
            pInteraction(i,parglmo.ordInteractions(i,var)) = (size(find( sum(SSQInteractions(2:(nPerm*mtcc + 1),i,1:var),3) ...
                >= sum(SSQInteractions(1, i, parglmo.ordInteractions(i,1:var)))), 1) + 1)/(nPerm*mtcc+1);
        end
    end
end
        
% Multiple test correction for several factors/interactions
parglmo.p = [pFactor' pInteraction'];
if mtcc > 1
    switch fmtc
        case 1 % Bonferroni 
            parglmo.p = min(1,parglmo.p * mtcc); 

        case 2 % Holm/Hochberg
            [~,indx] = sort(min(parglmo.p),'ascend');
            for ind = 1 : mtcc 
                parglmo.p(:,indx(ind)) = min(1,parglmo.p(:,indx(ind)) * (mtcc-ind+1));
            end

        case 3 % Benjamini & Hochberg
            [mv,indmv] = min(parglmo.p);
            [~,indx] = sort(mv,'ascend');
            parglmo.p(:,indx(mtcc)) = parglmo.p(:,indx(mtcc));
            for ind = mtcc-1 : -1 : 1 
                parglmo.p(indmv(indx(ind)),indx(ind)) = min(1,min(parglmo.p(indmv(indx(ind)),indx(ind)) * mtcc/ind,parglmo.p(indmv(indx(ind+1)))));
                parglmo.p(:,indx(ind)) = parglmo.p(:,indx(ind))*parglmo.p(indmv(indx(ind)),indx(ind))/parglmo.p(indmv(indx(ind)),indx(ind));
            end

        case 4 % Q-value from Benjamini & Hochberg
            [mv,indmv] = min(parglmo.p);
            [~,indx] = sort(mv,'ascend');
            parglmo.p(:,indx(mtcc)) = parglmo.p(:,indx(mtcc));
            for ind = mtcc-1 : -1 : 1 
                parglmo.p(indmv(indx(ind)),indx(ind)) = min(1,min(parglmo.p(indmv(indx(ind)),indx(ind)) * mtcc/ind,parglmo.p(indmv(indx(ind+1)),indx(ind+1))));
                parglmo.p(:,indx(ind)) = parglmo.p(:,indx(ind))*parglmo.p(indmv(indx(ind)),indx(ind))/parglmo.p(indmv(indx(ind)),indx(ind));
            end
    end
end


%% ANOVA-like output table

name={'Mean'};
for f = 1 : parglmo.nFactors
    name{end+1} = sprintf('Factor %d',f);
end
for i = 1 : parglmo.nInteractions
     name{end+1} = sprintf('Interaction %s',strrep(num2str(parglmo.interactions{i}.factors),'  ','-'));
end
name{end+1} = 'Residuals';
name{end+1} = 'Total';
      
if parglmo.nInteractions
    SSQ = sum([SSQinter' permute(SSQFactors(1,:,:),[3 2 1]) permute(SSQInteractions(1,:,:),[3 2 1]) SSQresiduals(1,:)' SSQX'],1);
else
    SSQ = sum([SSQinter' permute(SSQFactors(1,:,:),[3 2 1]) SSQresiduals(1,:)' SSQX'],1);
end
par = [mean(parglmo.effects,1) 100];
DoF = [1 parglmo.df parglmo.dfint parglmo.Rdf parglmo.Tdf];
MSQ = SSQ./DoF;
F = [nan max(tsFactors(1,:,:),[],3) max(tsInteractions(1,:,:),[],3) nan nan];
pValue = [nan min(parglmo.p,[],1) nan nan];

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    T.data = [name'; SSQ'; par'; DoF'; MSQ'; F'; pValue'];
    T.labels = {'Source','SumSq','AvPercSumSq','df','MeanSq','MaxF','minPvalue'};
else
    T = table(name', SSQ', par', DoF', MSQ', F', pValue','VariableNames', {'Source','SumSq','AvPercSumSq','df','MeanSq','MaxF','minPvalue'});
end

end

 

    