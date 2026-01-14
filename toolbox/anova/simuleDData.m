function Xout = simuleDData(X, F, varargin)

% Simulation of designed data. 
%
% Xm = simuleDData(X, F)   % minimum call
%
% See also: powercurve, parglm, asca, apca, parglmVS, parglmMC, createDesign
%
%
% INPUTS:
%
% X: (struct) to generate Population PCs
%       X.N: [1x1] number of rows
%       X.M: [1x1] number of columns
%       X.k: [1x(F+I)] standard deviation coefficients
%
% F: [NxF] design matrix, cell or array, where columns correspond to 
% factors and rows to levels
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
% 'RandomGen': (func or cell with F+I+1 functions) random generator (@randn by default) suggested alternatives
%   - @(N,M)simuleMV(N,M,'LevelCorr',8): multivariate correlated with level 8 (other values may be used)
%   - @rand: uniform (moderately non-normal)
%   - @(N,M)exprnd(1,N,M).^3: very non-normal
%   - {@randn,@rand}: normal for the single factor, uniform for the
%   residuals
%
% 'RandomGenC': (func) random generator in effect size coefficients (@()1 by default)
%    - To generate randomness use, e.g., @()0.1*randn+1 
%
% 'Theta': [1xT] For type equal to 1, theta controls the compromise of 
%   true significance vs random (0:0.1:1 by default). For type equal to 2, 
%   theta controls the number of replicates (1:10 by default)
%
% 'Ordinal': [1xF] whether factors are nominal or ordinal
%       0: nominal (default)
%       1: ordinal
%
% 'Random': [1xF] whether factors are fixed or random
%       0: fixed 
%       1: random (default, so far the simulation is always random)
%
% 'Coding': [1xF] type of coding of factors
%       0: sum/deviation coding (default)
%       1: reference coding (reference is the last level)
%
% 'Nested': [nx2] pairs of neted factors, e.g., if factor 2 is nested in 1,
%   and 3 in 2, then nested = [1 2; 2 3]
%
% 'Stable': [bool] maintain a fixed seed for reproducibility (false by default)
%
%
% OUTPUTS:
%
% Xout: [1xT] Each element is a [NxM] billinear data set for model fitting, 
% where each row is a measurement, each column a variable.
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Simulated Power Curve (type 1) from effect coefficients, two factors 
%   with 4 and 3 levels, and 4 replicates, with interaction. The relative 
%   effect (standard deviation) is .1 and .2 for the first and second 
%   factor and .3 for the interaction, respectively. Please, note this
%   example takes several minutes to execute.
%
% reps = 4;
% levels = {[1,2,3,4],[1,2,3]};
% 
% F = createDesign(levels,'Replicates',reps);
% 
% X.N = size(F,1);
% X.M = 400;
% X.k = [.1,.2,.3];
% 
% Xout = simuleDData(X, F, 'Model',{[1 2]});
%
% parglm(Xout{1},F,'Model',{[1 2]}) % Data with null effect
% parglm(Xout{end},F,'Model',{[1 2]}) % Data with maximum effect
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 14/Jan/2026
% Dependencies: Matlab R2017b, MEDA v1.10
%
% Copyright (C) 2026  University of Granada, Granada
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

nFactors = size(F,2);                 % number of factors

N = X.N;
M = X.M;

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
    if isstruct(X), tip = 1; 
    else, tip = 2;
    end
addParameter(p,'Model','linear');
addParameter(p,'RandomGen',@randn);
addParameter(p,'RamdonGenC',@()1);
addParameter(p,'Theta',[]);
addParameter(p,'Ordinal',zeros(1,size(F,2)));
addParameter(p,'Random',ones(1,size(F,2))); 
addParameter(p,'Coding',zeros(1,size(F,2)));
addParameter(p,'Nested',[]);
addParameter(p,'Stable',false);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
model = p.Results.Model;
randg = p.Results.RandomGen;
randgC = p.Results.RamdonGenC;
theta = p.Results.Theta;
ordinal = p.Results.Ordinal;
random = p.Results.Random;
coding = p.Results.Coding;
nested = p.Results.Nested;
stable = p.Results.Stable;

if stable, rng(0); end

if isempty(theta)
    theta = 0:0.1:1; 
end
theta = sort(theta,'ascend');

if isempty(model), model = 'linear'; end

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
assert (isequal(size(ordinal), [1 size(F,2)]), 'Dimension Error: parameter ''Ordinal'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(random), [1 size(F,2)]), 'Dimension Error: parameter ''Random'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(coding), [1 size(F,2)]), 'Dimension Error: parameter ''Coding'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
nInteractions      = length(interactions);      % number of interactions
nFactors           = size(F,2);                 % number of factors

if isscalar(randg)
    randv = repmat({randg},1,nFactors+nInteractions+1);
else
    randv = randg;
end


% Tranform Design Matrix

Fold = F;
F = zeros(size(F));
for f = 1:size(F,2)
    uF = unique(Fold(:,f),'stable');
    for i = 1:length(uF)
        F(find(ismember(Fold(:,f),uF(i))),f) = i;
    end
end

% Create Coding Matrix

n = 1; % This is necessary for indirect orthogonalization, to avoid leakage of variance to the residuals
D = ones(size(X,1),1);

powercurveo.Dvars = [0];

for f = 1 : nFactors
    if ordinal(f)
        D(:,n+1) = preprocess2D(F(:,f),'Preprocessing',1);
        powercurveo.factors{f}.Dvars = n+1;
        n = n + 1;
        powercurveo.factors{f}.order = 1;
        powercurveo.Dvars(n+1) = 1;
    else
        if isempty(nested) || isempty(find(nested(:,2)==f)) % if not nested
            powercurveo.factors{f}.factors = [];
            uF = unique(F(:,f));
            powercurveo.nLevels(f) = length(uF);
            for i = 2:length(uF)
                D(find(ismember(F(:,f),uF(i))),n+i-1) = 1;
            end
            powercurveo.factors{f}.Dvars = n+(1:length(uF)-1);
            if coding(f) == 1
                D(find(ismember(F(:,f),uF(1))),powercurveo.factors{f}.Dvars) = 0;
            else
                D(find(ismember(F(:,f),uF(1))),powercurveo.factors{f}.Dvars) = -1;
            end
            n = n + length(uF) - 1;
            powercurveo.factors{f}.order = 1;
            powercurveo.Dvars(powercurveo.factors{f}.Dvars) = powercurveo.factors{f}.order;
        else % if nested
            ind = find(nested(:,2)==f);
            ref = nested(ind,1);
            powercurveo.factors{f}.factors = [ref powercurveo.factors{ref}.factors];
            urF = unique(F(:,ref));
            powercurveo.nLevels(f) = 0;
            powercurveo.factors{f}.Dvars = [];
            for j = 1:length(urF)
                rind = find(ismember(F(:,ref),urF(j)));
                uF = unique(F(rind,f));
                powercurveo.nLevels(f) = powercurveo.nLevels(f) + length(uF);
                for i = 2:length(uF)
                    D(rind(find(ismember(F(rind,f),uF(i)))),n+i-1) = 1;
                end
                powercurveo.factors{f}.Dvars = [powercurveo.factors{f}.Dvars n+(1:length(uF)-1)];
                if coding(f) == 1
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = 0;
                else
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = -1;
                end
                n = n + length(uF) - 1;
            end   
            powercurveo.factors{f}.order = powercurveo.factors{ref}.order + 1;
            powercurveo.Dvars(powercurveo.factors{f}.Dvars) = powercurveo.factors{f}.order;
        end
    end
end

for i = 1 : nInteractions
    Dout = computaDint(interactions{i},powercurveo.factors,D);
    D = [D Dout];
    powercurveo.interactions{i}.Dvars = n+1:size(D,2);
    powercurveo.interactions{i}.factors = interactions{i};
    n = size(D,2);
    powercurveo.interactions{i}.order = max(powercurveo.factors{interactions{i}(1)}.order,powercurveo.factors{interactions{i}(2)}.order) + 1;
    powercurveo.Dvars(powercurveo.interactions{i}.Dvars) = powercurveo.interactions{i}.order;
end

pD =  pinv(D'*D)*D';
    
% Degrees of freedom
Tdf = N;   

mdf = 0;
Rdf = Tdf-mdf;
for f = 1 : nFactors
    if ordinal(f)
        df(f) = 1;
    else
        df(f) = length(powercurveo.factors{f}.Dvars);
    end
    Rdf = Rdf-df(f);
end
dfint = [];
for i = 1 : nInteractions
    dfint(i) = prod(df(powercurveo.interactions{i}.factors));
    Rdf = Rdf-dfint(i);
end
if Rdf < 0
    disp('Warning: degrees of freedom exhausted');
    return
end

powercurveo.coeffs = X.k;
powercurveo.rescoef = 1;


%% Compute the Power Curve
 
 
for f = 1 : nFactors
    randg = randv{f};
    if ordinal(f)
        powercurveo.factors{f}.matrix = randg(N,M);
        powercurveo.factors{f}.matrix = sqrt(N)*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');
    else
        powercurveo.factors{f}.matrix = randg(powercurveo.nLevels(f),M);
        powercurveo.factors{f}.matrix = sqrt(powercurveo.nLevels(f))*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');

        if isempty(nested) || isempty(find(nested(:,2)==f)) % order 1
            powercurveo.factors{f}.matrix = powercurveo.factors{f}.matrix(F(:,f),:);
        else % nested
            mati = powercurveo.factors{f}.matrix;
            Fi = F(:,[powercurveo.factors{f}.factors f]);
            powercurveo.factors{f}.matrix = [];
            uF = unique(Fi,'rows');
            for n = 1: size(uF,1)
                ind = find(uF(n,1)==Fi(:,1)&uF(n,2)==Fi(:,2));
                powercurveo.factors{f}.matrix(ind,:) = ones(length(ind),1)*mati(n,:);
            end
        end
    end
end

for i = 1 : nInteractions
    randg = randv{nFactors+i};
    mati = randg(prod(powercurveo.nLevels(powercurveo.interactions{i}.factors)),M);
    mati = sqrt(size(mati,1))*mati/norm(mati,'fro');
    Fi = F(:,powercurveo.interactions{i}.factors);
    powercurveo.interactions{i}.matrix = [];
    uF = unique(Fi,'rows');
    for n = 1: size(uF,1)
        ind = find(ismember(Fi, uF(n,:), 'rows'));
        powercurveo.interactions{i}.matrix(ind,:) = ones(length(ind),1)*mati(n,:);
    end
end

Xstruct = zeros(N,M);
for f = 1 : nFactors
    if ~ordinal(f)
        Xstruct = Xstruct + randgC() * powercurveo.coeffs(f) * powercurveo.factors{f}.matrix;
    end
end

for i = 1 : nInteractions
    Xstruct = Xstruct + randgC() * powercurveo.coeffs(nFactors+i) * powercurveo.interactions{i}.matrix;
end

randg = randv{end};
Xnoise = randg(N,M);
Xnoise = randgC() * powercurveo.rescoef * sqrt(N)*Xnoise/norm(Xnoise,'fro');

for a = 1:length(theta)

    Xm = Xnoise; %Xm = (1-theta(a))*Xnoise;
    Xm = Xm + theta(a)*Xstruct;

    for f = 1 : nFactors
        if ordinal(f)
            Xm = Xm + randgC() * powercurveo.coeffs(f) * powercurveo.factors{f}.matrix;
        end
    end

    Xout{a} = Xm;

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









