function [bel,states,stree] = gia(map,varargin)

% Group identification algorithm (gia) to generate groups of variables from 
% a correlation map of variables x variables. It should be noted that this
% algorithm is a heuristic, and not all possible groups may be identified.
%
% bel = gia(map)   % minimum call
% [bel,states,stree] = gia(map,'Gamma',gamma,'MinSize',siz,'Stree',stree)   % complete call
%
%
% INPUTS:
%
% map: [MxM] correlation matrix. Values should be between -1 and 1.
%
% Optional INPUTS (parameters):
%
% 'Gamma': [1x1] correlation threshold to identify groups (0.7 by default)
%
% 'MinSize': [1x1] Integer with the minimum size of groups (2 by default)
%
% 'Stree': [struct] tree with GIA division, if previously executed. This
%   structure makes reiterative GIA computations faster (empty by default)
%   - tree: cell with division tree 
%   - indm: number of variables above a given threshold
%   - index: value for each tree division
%
%
% OUTPUTS:
%
% bel: {Mx1} Cell with the list of groups each variable belongs to.
%
% states: {Sx1} Cell with the groups of variables.
%
% stree: [struct] tree with GIA division, if previously executed. This
%   structure makes reiterative GIA computations faster.
%   - tree: cell with division tree 
%   - indm: number of variables above a given threshold
%   - index: value for each tree division
%
%
% EXAMPLE OF USE: Random data. Check the value in states:
%
% X = simuleMV(20,10,'LevelCorr',8);
% pcs = 1:3;
% map = meda_pca(X,'Pcs',pcs);
% [bel,states] = gia(map,'Gamma',0.3);
%
%
% EXAMPLE OF USE: Checking metaparameter:
%
% X = simuleMV(20,100,8);
% pcs = 1:3;
% map = meda_pca(X,'Pcs',pcs);
% C = [0.05:0.05:0.95];
% 
% [belv{1},statesv{1},stree] = gia(map,'Gamma',C(1));
% S = 0;
% for j=1:length(statesv{1}), S = S + length(statesv{1}{j}); end;
% disp(['There are ',num2str(length(statesv{1})),' groups with mean size ' ,num2str(S/length(statesv{1})), ' for C = ', num2str(C(1))])
% 
% for i=2:length(C)
%   [belv{i},statesv{i}] = gia(map,'Gamma',C(i),'Stree',stree);
%   S = 0;
%   for j=1:length(statesv{i}), S = S + length(statesv{i}{j}); end;
%   disp(['There are ',num2str(length(statesv{i})),' groups with mean size ' ,num2str(S/length(statesv{i})), ' for C = ', num2str(C(i))])
% end
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 22/Apr/24.
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
M = size(map, 1);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Gamma',0.7);  
addParameter(p,'MinSize',2);
addParameter(p,'Stree',{});            
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
gamma = p.Results.Gamma;
siz = p.Results.MinSize;
stree = p.Results.Stree;

% Avoid gamma with just 1
%if gamma==1, gamma = 1 -1e-10; end;

% Validate dimensions of input data
assert (isequal(size(map), [M M]), 'Dimension Error: parameter ''map'' must be M-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(gamma), [1 1]), 'Dimension Error: parameter ''Gamma'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(siz), [1 1]), 'Dimension Error: parameter ''MinSize'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(map<-1-1e-10)), 'Value Error: parameter ''map'' must contain values between -1 and 1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(map>1+1e-10)), 'Value Error: parameter ''map'' must contain values between -1 and 1. Type ''help %s'' for more info.', routine(1).name);
assert (gamma>=0 && gamma<=1, 'Value Error: parameter ''Gamma'' must be between 0 and 1. Type ''help %s'' for more info.', routine(1).name);
assert (siz<M, 'Value Error: parameter ''MinSize'' must be below M. Type ''help %s'' for more info.', routine(1).name);
assert (siz>0, 'Value Error: parameter ''MinSize'' must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(siz), siz), 'Value Error: parameter ''MinSize'' must be an integer. Type ''help %s'' for more info.', routine(1).name);


%% Main code

if gamma == 0
    states = {1:M};
    bel = ones(M,1);
    return
end

if gamma == 1
    states = {};
    bel = cell(M,1);
    return
end

map2 = map;

if isempty(stree)
    indm = find(max(map)>gamma); % reduce map
else
    indm = stree.indm;
end
M2 = M;
M = length(indm);

if isempty(stree)
    
    map = abs(map);
    map = map(indm,indm);
    map_v = map - tril(map) + diag(diag(map)); % triangular inferior part is set to zero
    map_v = map_v(:);

    columns = ones(M,1)*(1:M);
    rows = columns';
    columns_v = columns(:);
    rows_v = rows(:);

    states = {};
    bel = cell(M,1);
    
    ind = find(map_v==max(map_v),1);
    tree = {};
    index = [];
    while map_v(ind)>=gamma
        
        val = map_v(ind);
        
        r = rows_v(ind); % r is always minor than c
        c = columns_v(ind);
        
        if c==r
            if isempty(bel{r})
                states{end+1} = r;
                bel{r} = [bel{r} length(states)];
            end
        else
            bel_r = setdiff(bel{r},bel{c});
            
            mod = [];
            
            for i=1:length(bel_r) % Should c be added to an state of r?
                state_rec = states{bel_r(i)};
                
                j = 1;
                while j<=length(state_rec) && map(state_rec(j),c)> val-1e-10
                    j = j+1;
                end
                
                if j > length(state_rec),
                    states{bel_r(i)} = [state_rec c];
                    bel{c} = [bel{c} bel_r(i)];
                    mod(end+1) = bel_r(i);
                end             
            end
                     
            dup = [];
            modbel = bel{c};
            for i=1:length(mod)
                l = find(ismember(modbel,mod(i)),1);
                for j=[1:l-1 l+1:length(modbel)]
                    if isempty(find(map(states{mod(i)},states{modbel(j)}) < val,1)),
                        dup = [dup modbel(j)];
                        modbel(j) = [];
                        break
                    end
                end,
            end
                
            if ~isempty(dup)
                indd = 1:length(states);
                states(dup) = [];
                for k=1:length(dup)
                    indd((dup(k)+1):end) = indd((dup(k)+1):end)-1;
                end
                indd(dup) = 0;
                for k=1:M
                    bel{k} = indd(bel{k});
                    bel{k}(find(bel{k}==0)) = [];
                end
            end
            
            bel_c = setdiff(bel{c},bel{r});
            
            for i=1:length(bel_c) % Should r be added to an state of c?
                state_rec = states{bel_c(i)};
                
                j = 1;
                while j<=length(state_rec) && map(state_rec(j),r)> val-1e-10,
                    j = j+1;
                end
                
                if j > length(state_rec)
                    states{bel_c(i)} = [state_rec r];
                    bel{r} = [bel{r} bel_c(i)];
                    mod(end+1) = bel_c(i);
                    
                end
            end
            
            if isempty(mod) && isempty(intersect(bel{r},bel{c}))   
                states{end+1} = [r,c];
                bel{r} = [bel{r} length(states)];
                bel{c} = [bel{c} length(states)];
                mod = length(states); 
            end

        end
        
        
        
        map_v(ind) = 0;
        ind = find(map_v==max(map_v),1); 
        
        index(end+1) = val;
        tree{end+1} = {states bel};
    end
else
    j = find(stree.index>=gamma);
    if isempty(j)
        states = {};
        bel = cell(M,1);
        return
    else
        states = stree.tree{j(end)}{1};
        bel = stree.tree{j(end)}{2};
    end
end

vs = [];
om = ceil(log10(M));
j=1;
for i=1:length(states)
    is = 0;
    l = 1;
    while ~is && l<j
        if isequal(stateso{l},sort(indm(states{i})))
            is = true;
        end
        l = l+1;
    end
    if ~is && length(states{i})>=siz
        stateso{j} = sort(indm(states{i}));
        vs(j) = sum(stateso{j}.*(10.^(-1*(1:om:om*length(states{i})))));  
        j = j+1;
    end
end    

[kk,ind] = sort(vs);
if isempty(ind)
    states = {};
else
    states = stateso(ind);
end

bel = cell(M2,1);
for j=1:length(states)
    for i=1:length(states{j})
        bel{states{j}(i)} = [bel{states{j}(i)} j];
    end
end

if isempty(stree)
    stree.tree = tree;
    stree.indm = indm;
    stree.index = index;
end
