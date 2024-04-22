function [D, parglmo, anovast] = codglm(F, parglmi, varargin)
%%%Preguntar anovast.model como parseo opciones dentro de anovast
% Compute coding matrix from a design matrix for General Linear Models.
%
% Related routines: parglm, asca, apca, parglmVS, parglmMC, create_design
%
% D = codglm(F, parglmi)   % minimum call
% [D, parglmo, anovast] = codglm(F, parglmi, anovast)   % complete call
%
%
% INPUTS:
%
% F: [NxF] design matrix, cell or array, where columns correspond to 
% factors and rows to levels.
%
% parglmi (structure): structure with the number of factors and the number 
% of levels for each of them.
%
% Optional INPUTS (parameters):
%
% 'Anovast' (structure): structure with the anova choices.
%
%   'Model': This paremeter is similar to 'model' of anovan. It could be:
%       'linear': only main effects are provided (by default)
%       'interaction': two order interactions are provided
%       'full': all potential interactions are provided
%       [1x1]: maximum order of interactions considered
%       [ix2]: array with two order interactions
%       cell: with each element a vector of factors
%
% 	'Ordinal': [1xF] whether factors are nominal or ordinal
%       0: nominal (default)
%       1: ordinal
%
%   'Coding': [1xF] type of coding of factors
%       0: sum/deviation coding (default)
%       1: reference coding (reference is the last level)
%
%   'Nested': [nx2] pairs of neted factors, e.g., if factor 2 is nested in 1,
%   and 3 in 2, then nested = [1 2; 2 3]
%
%
% OUTPUTS:
%
% D [NxC]: Coding matrix
%
% parglmo (structure): structure with the factor and interaction
% matrices, p-values (corrected, depending on fmtc) and explained variance 
%
% anovast (structure): structure with the anova choices.
%
%
% EXAMPLE OF USE (copy and paste the code in the command line) 
%   Two factors, with 4 and 3 levels and 4 replicates
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
% 
% F = create_design(levels,'Replicates',reps);
% 
% D = codglm(F)
%
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 22/Apr/2024
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

n_factors = parglmi.n_factors;                 % number of factors
levels = parglmi.levels;


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Anovast',[]);             
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
anovast = p.Results.Anovast;
    
if ~isfield(anovast,'model') || isempty(anovast.model), anovast.model = 'linear'; end;

if isequal(anovast.model,'linear')
    interactions = [];
end  
    
if isequal(anovast.model,'interaction')
    interactions = allinter(n_factors,2);
end    

if isequal(anovast.model,'full')
    interactions = allinter(n_factors,n_factors);
end    

if isnumeric(anovast.model) && isscalar(anovast.model) && anovast.model >= 2 && anovast.model <= n_factors
        interactions = allinter(n_factors,anovast.model);
end    

if isnumeric(anovast.model) && ~isscalar(anovast.model)
        interactions = {anovast.model};
end    

if iscell(anovast.model), interactions = anovast.model; end
    
if ~isfield(anovast,'ordinal') || isempty(anovast.ordinal), anovast.ordinal = zeros(1,n_factors); end;
if ~isfield(anovast,'coding') || isempty(anovast.coding), anovast.coding = zeros(1,n_factors); end;
if ~isfield(anovast,'nested'), anovast.nested = []; end;

% Validate dimensions of input data
assert (isequal(size(anovast.ordinal), [1 n_factors]), 'Dimension Error: parameter ''Ordinal''  must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(anovast.coding), [1 n_factors]), 'Dimension Error: parameter ''Coding''  must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
n_interactions      = length(interactions);      % number of interactions

% In column space
parglmo.factors                = cell(n_factors,1);
parglmo.interactions           = cell(n_interactions,1);

% Make structure with general 'variables'
parglmo.n_factors      = n_factors;
parglmo.levels         = levels;
parglmo.n_interactions = n_interactions;

% Create Design Matrix
n = 1;
D = ones(size(F,1),1);

for f = 1 : n_factors
    if anovast.ordinal(f)
        D(:,n+1) = preprocess2D(F(:,f),1);
        parglmo.factors{f}.Dvars = n+1;
        n = n + 1;
        parglmo.factors{f}.order = 1;
    else
        if isempty(anovast.nested) || isempty(find(anovast.nested(:,2)==f)) % if not nested
            uF = unique(levels{f});
            parglmo.n_levels(f) = length(uF);
            for i = 2:length(uF)
                D(find(ismember(F(:,f),uF(i))),n+i-1) = 1;
            end
            parglmo.factors{f}.Dvars = n+(1:length(uF)-1);
            if anovast.coding(f) == 1
                D(find(ismember(F(:,f),uF(1))),parglmo.factors{f}.Dvars) = 0;
            else
                D(find(ismember(F(:,f),uF(1))),parglmo.factors{f}.Dvars) = -1;
            end
            n = n + length(uF) - 1;
            parglmo.factors{f}.order = 1;
        else % if nested
            ind = find(anovast.nested(:,2)==f);
            ref = anovast.nested(ind,1);
            urF = unique(levels{ref});
            parglmo.n_levels(f) = 0;
            parglmo.factors{f}.Dvars = [];
            for j = 1:length(urF)
                rind = find(ismember(levels{f}(:,1),urF(j)));
                uF = unique(levels{f}(rind,2));
                parglmo.n_levels(f) = parglmo.n_levels(f) + length(uF);
                for i = 2:length(uF)
                    D(rind(find(ismember(F(rind,f),uF(i)))),n+i-1) = 1;
                end
                parglmo.factors{f}.Dvars = [parglmo.factors{f}.Dvars n+(1:length(uF)-1)];
                if anovast.coding(f) == 1
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = 0;
                else
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = -1;
                end
                n = n + length(uF) - 1;
            end   
            parglmo.factors{f}.order = parglmo.factors{ref}.order + 1;
        end
    end
end

for i = 1 : n_interactions
    Dout = computaDint(interactions{i},parglmo.factors,D);
    D = [D Dout];
    parglmo.interactions{i}.Dvars = n+1:size(D,2);
    parglmo.interactions{i}.factors = interactions{i};
    n = size(D,2);
    parglmo.interactions{i}.order = max(parglmo.factors{interactions{i}(1)}.order,parglmo.factors{interactions{i}(2)}.order) + 1;
end
        
end

%% Auxiliary function for interactions

function interactions = allinter(nF,order)
    
    if order > 2
        interactions = allinter(nF,order-1);
        for i = 1:length(interactions)
            for j = max(interactions{i})+1:nF
                interactions{end+1} = [interactions{i} j];
            end
        end
    else
        interactions = {};
        for i = 1:nF
            for j = i+1:nF
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

