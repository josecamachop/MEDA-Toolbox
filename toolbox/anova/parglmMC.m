function [T, parglmo] = parglmMC(X, F, varargin)

% Parallel General Linear Model to obtain factor and interaction matrices 
% in a crossed experimental design and permutation test for univariate 
% statistical significance through multiple-test correction methods that 
% allows variable selection. Missing data is considered.
%
% T = parglmMC(X, F)   % minimum call
%
% See also: parglm, parglmVS, asca, apca, createDesign
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
%       1: mean centering
%       2: autoscaling (default)
%
% 'Permutations': [1x1] number of permutations (1000 by default)
%
% 'Ts': [1x1] Use SSQ (0) or the F-value (otherwise, by default) as test statistic  
%
% 'Ordinal': [1xF] whether factors are nominal or ordinal
%       0: nominal
%       1: ordinal
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
% parglmoMC (structure): structure with the factor and interaction
% matrices, p-values (corrected, depending on mtc and fmtc) and explained 
% variance  
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Random data, one factor and two levels, three variables with information 
% on the factor. This example takes long to compute, you  may reduce the 
% number of variables or permutations.
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
% [TVS, parglmoVS] = parglmVS(X, class); % With multivariate variable selection
% rng(S);
% [TMC, parglmoMC] = parglmMC(X, class, 'Permutations', 500); % With variable selection through multiple testing correction
% 
% h = figure; hold on
% plot([1 nVars],-log10([parglmo.p parglmo.p]),'b-.')
% plot(-log10(parglmoVS.p(parglmoVS.ordFactors)),'g-o')
% plot(-log10(parglmoMC.p(parglmoMC.ordFactors)),'k-')
% plot([0,size(X,2)],-log10([0.01 0.01]),'r--')
% legend('ASCA','VASCA','BH-FDR','alpha=0.01','Location','southeast')
% a=get(h,'CurrentAxes');
% set(a,'FontSize',14)
% ylabel('-log10(p-values)','FontSize',18)
% xlabel('Responses in selected order','FontSize',18)
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 18/Nov/2024
%
% Copyright (C) 2024  Universidad de Granada
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
addParameter(p,'Ts',1); 
addParameter(p,'Ordinal',zeros(1,size(F,2))); 
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
mtc = p.Results.Mtc;
fmtc = p.Results.Fmtc;
coding = p.Results.Coding;
nested = p.Results.Nested;

if isequal(model,'linear')
    interactions = [];
end  
    
f = 1:nFactors;
if ~isempty(nested), f(nested(:,2)) = []; end
if isequal(model,'interaction')
    interactions = allinter(f,2);
end    

if isequal(model,'full')
    interactions = allinter(f,length(f));
end    

if isnumeric(model) && isscalar(model) && model >= 2 && model <= nFactors
        interactions = allinter(f,model);
end    

if isnumeric(model) && ~isscalar(model)
    for i = 1:size(model,1)
        interactions{i} = model(i,:);
    end
end    

if iscell(model), interactions = model; end

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(nPerm), [1 1]), 'Dimension Error: parameter ''Permutations'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ts), [1 1]), 'Dimension Error: parameter ''Ts'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(mtc), [1 1]), 'Dimension Error: parameter ''Mtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(fmtc), [1 1]), 'Dimension Error: parameter ''Fmtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(coding), [1 size(F,2)]), 'Dimension Error: parameter ''Coding'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);

%% Main code

nInteractions      = size(interactions,1);              % number of interactions
mtcc                = nFactors + nInteractions;        % correction for the number of tests
tsFactors         = zeros(nPerm*M+1,nFactors,M);       % sum of squares for factors
tsInteractions    = zeros(nPerm*M+1,nInteractions,M);  % sum of squares for interactions
FFactors           = zeros(nFactors,M);       % F-value 
FInteractions      = zeros(nInteractions,M);  % F-value 
pFactor            = zeros(nFactors,M);                % p-values factors
pInteraction       = zeros(nInteractions,M);           % p-values interactions

% In column space
parglmo.factors                = cell(nFactors,1);
parglmo.interactions           = cell(nInteractions,1);

% preprocess the data
[Xs,m,dt] = preprocess2D(X,'Preprocessing',prep);
X = X./(ones(size(X,1),1)*dt);

% Make structure with unchanging 'variables'
parglmo.data           = X;
parglmo.prep           = prep;
parglmo.scale           = dt;
parglmo.design         = F;
parglmo.nFactors      = nFactors;
parglmo.nInteractions = nInteractions;
parglmo.nPerm          = nPerm;
parglmo.ts              = ts;
parglmo.ordinal         = ordinal;
parglmo.mtc             = mtc;
parglmo.fmtc            = fmtc;
parglmo.coding          = coding;
parglmo.nested          = nested; 

% Create Design Matrix
n = 1;
D = ones(size(X,1),1);

for f = 1 : nFactors
    if ordinal(f)
        D(:,n+1) = preprocess2D(F(:,f),'Preprocessing',1);
        parglmo.factors{f}.Dvars = n+1;
        n = n + 1;
    else
        if isempty(nested) || isempty(find(nested(:,2)==f)) % if not nested
            uF = unique(F(:,f));
            parglmo.nLevels(f) = length(uF);
            for i = 2:length(uF)
                D(find(ismember(F(:,f),uF(i))),n+i-1) = 1;
            end
            parglmo.factors{f}.Dvars = n+(1:length(uF)-1);
            if coding(f) == 1
                D(find(ismember(F(:,f),uF(1))),parglmo.factors{f}.Dvars) = 0;
            else
                D(find(ismember(F(:,f),uF(1))),parglmo.factors{f}.Dvars) = -1;
            end
            n = n + length(uF) - 1;
        else % if nested
            ind = find(nested(:,2)==f);
            ref = nested(ind,1);
            urF = unique(F(:,ref));
            parglmo.nLevels(f) = 0;
            parglmo.factors{f}.Dvars = [];
            for j = 1:length(urF)
                rind = find(ismember(F(:,ref),urF(j)));
                uF = unique(F(rind,f));
                parglmo.nLevels(f) = parglmo.nLevels(f) + length(uF);
                for i = 2:length(uF)
                    D(rind(find(ismember(F(rind,f),uF(i)))),n+i-1) = 1;
                end
                parglmo.factors{f}.Dvars = [parglmo.factors{f}.Dvars n+(1:length(uF)-1)];
                if coding(f) == 1
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = 0;
                else
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = -1;
                end
                n = n + length(uF) - 1;
            end   
        end
    end
end

for i = 1 : nInteractions
    Dout = computaDint(interactions{i},parglmo.factors,D);
    D = [D Dout];
    parglmo.interactions{i}.Dvars = n+1:size(D,2);
    parglmo.interactions{i}.factors = interactions{i};
    n = size(D,2);
end

% Degrees of freedom
Tdf = size(X,1);      
Rdf = Tdf-1;
for f = 1 : nFactors
    if ordinal(f)
        df(f) = 1;
    else
        df(f) = length(parglmo.factors{f}.Dvars);
    end
    Rdf = Rdf-df(f);
end
dfInt = [];
for i = 1 : nInteractions
    dfInt(i) = prod(df(parglmo.interactions{i}.factors));
    Rdf = Rdf-dfInt(i);
end
if Rdf < 0
    disp('Warning: degrees of freedom exhausted');
    return
end

% Handle missing data 
[r,c]=find(isnan(X));
Xnan = X;
ru = unique(r);
for i=1:length(ru)
    ind = find(r==ru(i));
    ind2 = find(sum((D-ones(size(D,1),1)*D(r(ind(1)),:)).^2,2)==0);
    for j=1:length(c(ind))
        ind3 = find(isnan(X(ind2,c(ind(j)))));
        if length(ind2)>length(ind3)
            X(r(ind(j)),c(ind(j))) = nanmean(X(ind2,c(ind(j)))); % use conditional mean replacement
        else
            X(r(ind(j)),c(ind(j))) = nanmean(X(:,c(ind(j)))); % use unconditional mean replacement if CMR not possible
        end
    end
end
parglmo.data = X;
parglmo.Xnan = Xnan;

SSQX = sum(X.^2);

% GLM model calibration with LS, only fixed factors
pD =  pinv(D'*D)*D';
B = pD*X;
Xresiduals = X - D*B;
parglmo.D = D;
parglmo.B = B;
parglmo.mean = parglmo.D*parglmo.B(:,1);

% Create Effect Matrices
parglmo.inter = D(:,1)*B(1,:);
SSQinter = sum(parglmo.inter.^2);
SSQresiduals = sum(Xresiduals.^2);

SSQFactors = [];
FFactors = [];
for f = 1 : nFactors
    parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
    SSQFactors(f,:) = sum(parglmo.factors{f}.matrix.^2);
    FFactors(f,:) = (SSQFactors(f,:)/df(f))./(SSQresiduals/Rdf);
end

% Interactions
SSQInteractions = [];
FInteractions = [];
for i = 1 : nInteractions
    parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
    SSQInteractions(i,:) = sum(parglmo.interactions{i}.matrix.^2);
    FInteractions(i,:) = (SSQInteractions(i,:)/dfInt(i))./(SSQresiduals/Rdf);
end
    
if nInteractions
    parglmo.effects = 100*([SSQinter' SSQFactors' SSQInteractions' SSQresiduals']./(SSQX'*ones(1,2+nFactors+nInteractions)));
else
    parglmo.effects = 100*([SSQinter' SSQFactors' SSQresiduals']./(SSQX'*ones(1,2+nFactors+nInteractions)));
end
parglmo.residuals = Xresiduals;

% Select test statistic to order variables
if ts
    tsFactors(1,:,:) = FFactors;
    tsInteractions(1,:,:) = FInteractions;
else
    tsFactors(1,:,:) = SSQFactors;
    tsInteractions(1,:,:) = SSQInteractions;
end
    
% Permutations
parfor j = 1 : (nPerm * M) % Increase the number of permutations to perform MTC
    
    perms = randperm(size(Xnan,1)); % permuted data (permute whole data matrix)
         
    X = Xnan(perms, :);
    [r,c]=find(isnan(X));
    ru = unique(r);
    for i=1:length(ru)
        ind = find(r==ru(i));
        ind2 = find(sum((D-ones(size(D,1),1)*D(r(ind(1)),:)).^2,2)==0);
        for f=1:length(c(ind))
            ind3 = find(isnan(X(ind2,c(ind(f)))));
            if length(ind2)>length(ind3)
                X(r(ind(f)),c(ind(f))) = nanmean(X(ind2,c(ind(f)))); % use conditional mean replacement
            else
                X(r(ind(f)),c(ind(f))) = nanmean(X(:,c(ind(f)))); % use unconditional mean replacement if CMR not possible
            end
        end
    end
      
    B = pD*X;
    Xresiduals = X - D*B;
    SSQresidualsP = sum(Xresiduals.^2);
    
    % Factors
    SSQFactorsP = [];
    FFactorsP = [];
    for f = 1 : nFactors
        factorMatrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQFactorsP(f,:) = sum(factorMatrix.^2);
        FFactorsP(f,:) = (SSQFactorsP(f,:)/df(f))./(SSQresidualsP/Rdf);
    end
    
    % Interactions
    SSQInteractionsP = [];
    FInteractionsP = [];
    for i = 1 : nInteractions
        interactsMatrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        SSQInteractionsP(i,:) = sum(interactsMatrix.^2);
        FInteractionsP(i,:) = (SSQInteractionsP(i,:)/dfInt(i))./(SSQresidualsP/Rdf);
    end
   
    % Select test statistic to order variables
    if ts
        tsFactors(1 + j,:,:) = FFactorsP;
        tsInteractions(1 + j,:,:) = FInteractionsP;
    else
        tsFactors(1 + j,:,:) = SSQFactorsP;
        tsInteractions(1 + j,:,:) = SSQInteractionsP;
    end
end        
   

% Calculate univariate p-values and order variables by relevance
for f = 1 : nFactors
    for var = 1 : M
        pFactor(f,var) = (size(find(tsFactors(2:end,f,var) ...
            >= tsFactors(1,f,var)) ,1) + 1)/(nPerm*M + 1);
    end
    [~,ordFactors(f,:)] = sort(pFactor(f,:),'ascend');
end
for i = 1 : nInteractions
    for var = 1 : M
        pInteraction(i,var) = (size(find(tsInteractions(2:end,i,var) ...
            >= tsInteractions(1,i,var)) ,1) + 1)/(nPerm*M + 1);
    end
    [~,ordInteractions(i,:)] = sort(pInteraction(i,:),'ascend');
end

parglmo.ordFactors = ordFactors;
if nInteractions>0
    parglmo.ordInteractions = ordInteractions;
end

parglmo.p = [pFactor' pInteraction'];

% Multiple test correction
switch mtc
    
    case 1 % Bonferroni
        if fmtc
            parglmo.p = min(1,parglmo.p * M * mtcc);
        else
            parglmo.p = min(1,parglmo.p * M);
        end
        
    case 2 % Holm/Hochberg
        
        if fmtc
            parglmo.p = parglmo.p;
            [~,indx] = sort(parglmo.p(:),'ascend');
            for ind = 1 : length(indx)
                parglmo.p(indx(ind)) = min(1,parglmo.p(indx(ind)) * ((M*mtcc)-ind+1));
            end
        else
            for f = 1 : nFactors
                for var = 1 : M
                    pFactor(f,ordFactors(f,var)) = min(1,pFactor(f,ordFactors(f,var)) * (M-var+1));
                end
            end
            for i = 1 : nInteractions
                for var = 1 : M
                    pInteraction(i,ordInteractions(i,var)) = min(1,pInteraction(i,ordInteractions(i,var)) * (M-var+1));
                end
            end
            parglmo.p = [pFactor' pInteraction'];
        end
        
    case 3 % Benjamini & Hochberg
        
        if fmtc
            [~,indx] = sort(parglmo.p(:),'ascend');
            for ind = length(indx)-1 : -1 : 1
                parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * (M*mtcc)/ind]);
            end
        else
            for f = 1 : nFactors
                for var = M-1 : -1 : 1
                    pFactor(f,ordFactors(f,var)) = min([1,pFactor(f,ordFactors(f,var)) * M/var]);
                end
            end
            for i = 1 : nInteractions
                for var = M-1 : -1 : 1
                    pInteraction(i,ordInteractions(i,var)) = min([1,pInteraction(i,ordInteractions(i,var)) * M/var]);
                end
            end
            parglmo.p = [pFactor' pInteraction'];
        end
        
end

if mtc <0 % Q-value     
    if fmtc
        parglmo.p = parglmo.p;
        [~,indx] = sort(parglmo.p(:),'ascend');
        parglmo.p(indx(end)) = parglmo.p(indx(end))*(-mtc);
        for ind = length(indx)-1 : -1 : 1
            parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * (-mtc*M*mtcc)/ind,parglmo.p(indx(ind+1))]);
        end
    else
        for f = 1 : nFactors
            pFactor(f,ordFactors(f,M)) = pFactor(f,ordFactors(f,M))*(-mtc);
            for var = M-1 : -1 : 1
                pFactor(f,ordFactors(f,var)) = min([1,pFactor(f,ordFactors(f,var)) * (-mtc) * M/var,pFactor(f,ordFactors(f,var+1))]);
            end
        end
        for i = 1 : nInteractions
            pInteraction(i,ordInteractions(i,M)) = pInteraction(i,ordInteractions(i,M))*(-mtc);
            for var = M-1 : -1 : 1
                pInteraction(i,ordInteractions(i,var)) = min([1,pInteraction(i,ordInteraction(i,var)) * (-mtc) * M/var,pInteraction(i,ordInteractions(i,var+1))]);
            end
        end
        parglmo.p = [pFactor' pInteraction'];
    end
end



%% ANOVA-like output table

name={'Mean'};
for f = 1 : nFactors
    name{end+1} = sprintf('Factor %d',f);
end
for i = 1 : nInteractions
    name{end+1} = sprintf('Interaction %s',strrep(num2str(parglmo.interactions{i}.factors),'  ','-'));
end
name{end+1} = 'Residuals';
name{end+1} = 'Total';
      
if nInteractions
    SSQ = sum([SSQinter' SSQFactors' SSQInteractions' SSQresiduals' SSQX'],1);
else
    SSQ = sum([SSQinter' SSQFactors' SSQresiduals' SSQX'],1);
end
par = [mean(parglmo.effects,1) 100];
DoF = [1 df dfInt Rdf Tdf];
MSQ = SSQ./DoF;
F = [nan max(FFactors,[],2)' max(FInteractions,[],2)' nan nan];
pValue = [nan min(parglmo.p) nan nan];

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    T.data = [name'; SSQ'; par'; DoF'; MSQ'; F'; pValue'];
    T.labels = {'Source','SumSq','AvPercSumSq','df','MeanSq','MaxF','minPvalue'};
else
    T = table(name', SSQ', par', DoF', MSQ', F', pValue','VariableNames', {'Source','SumSq','AvPercSumSq','df','MeanSq','MaxF','minPvalue'});
end

end

%% Auxiliary function for interactions

function interactions = allinter(factors,order)
    
    if order > 2
        interactions = allinter(factors,order-1);
        for i = 1:length(interactions)
            for j = factors(find(factors > max(interactions{i})))
                interactions{end+1} = [interactions{i} j];
            end
        end
    else
        interactions = {};
        for i = factors
            for j = factors(find(factors >i))
                interactions{end+1} = [i j];
            end
        end
    end
    
end
    
        
function Dout = computaDint(interactions,factors,D) % Compute coding matrix

    if length(interactions)>1
        deepD = computaDint(interactions(2:end),factors,D);
        Dout = [];
        for k = factors{interactions(1)}.Dvars
            for l = 1:size(deepD,2)
                Dout(:,end+1) = D(:,k).* deepD(:,l);
            end
        end
    else
        Dout = D(:,factors{interactions}.Dvars);
    end

end

