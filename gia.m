function [bel,states] = gia(map,gamma,siz)

% Group identification algorithm (gia) to generate groups of variables from 
% a correlation map of variables x variables.
%
% bel = gia(map)   % minimum call
% [bel,states] = gia(map,gamma,siz)   % complete call
%
%
% INPUTS:
%
% map: [MxM] correlation matrix. Values should be between -1 and 1.
%
% gamma: [1x1] correlation threshold to identify groups (0.7 by default)
%
% siz: [1x1] Integer with the minimum size of groups (2 by default)
%
%
% OUTPUTS:
%
% bel: {Mx1} Cell with the list of groups each variable belongs to.
%
% states: {Sx1} Cell with the groups of variables.
%
%
% EXAMPLE OF USE: Random data. Check the value in states:
%
% X = simuleMV(20,10,8);
% pcs = 1:3;
% map = meda_pca(X,pcs);
% [bel,states] = gia(map,0.3);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Apr/16.
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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
if nargin < 2 || isempty(gamma), gamma=0.7; end;
if nargin < 3 || isempty(siz), siz=2; end;

% Avoid gamma with just 1
if gamma==1, gamma = 1 -1e-10; end;

% Validate dimensions of input data
assert (isequal(size(map), [M M]), 'Dimension Error: 1st argument must be M-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(gamma), [1 1]), 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(siz), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(map<-1-1e-10)), 'Value Error: 1st argument must contain values between -1 and 1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(map>1+1e-10)), 'Value Error: 1st argument must contain values between -1 and 1. Type ''help %s'' for more info.', routine(1).name);
assert (gamma>=0 && gamma<=1, 'Value Error: 2nd argument must be between 0 and 1. Type ''help %s'' for more info.', routine(1).name);
assert (siz<M, 'Value Error: 3rd argument must be below M. Type ''help %s'' for more info.', routine(1).name);
assert (siz>0, 'Value Error: 3rd argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(siz), siz), 'Value Error: 3rd argument must be an integer. Type ''help %s'' for more info.', routine(1).name);


%% Main code

if gamma == 0,
    states = {1:M};
    bel = ones(M,1);
    return
end

map = abs(map);
indm = find(max(map)>gamma); % reduce map
map = map(indm,indm);
map_v = map - tril(map) + diag(diag(map)); % triangular inferior part is set to zero
map_v = map_v(:);
M2 = M;
M = size(map, 1);

columns = ones(M,1)*(1:M);
rows = columns';
columns_v = columns(:);
rows_v = rows(:);

states = {};
bel = cell(M,1);


ind = find(map_v==max(map_v),1);
while map_v(ind)>=gamma,
    
    r = rows_v(ind); % r is always minor than c
    c = columns_v(ind);
    
    if c==r,
        if isempty(bel{r})
            states{end+1} = r;
            bel{r} = [bel{r} length(states)];
        end
    else
        bel_r = setdiff(bel{r},bel{c});
        bel_c = setdiff(bel{c},bel{r});
        
        mod = 0;%zeros(length(bel_r),length(bel_c));
        mod2 = 0;
        states_n = {};
        
        for i = 1:length(bel_r), % Should two states be combined?
            for j = 1:length(bel_c),
                que = (map(states{bel_r(i)},states{bel_c(j)}) > gamma);
                if isempty(find(~que)),
                    mod(i,j) = 1;
                    states_n{end+1} = [states{bel_r(i)} states{bel_c(j)}]; 
                    map_v = reshape(map_v,M,M);
                    map_v(states{bel_r(i)},states{bel_c(j)}) = 0;
                    map_v(states{bel_c(j)},states{bel_r(i)}) = 0;
                    map_v = map_v(:);
                end
            end
        end
        
        if sum(sum(mod))
            ind_r = find(sum(mod,2));
            ind_c = find(sum(mod,1));
            
            delet = unique([bel_r(ind_r) bel_c(ind_c)]);
            
            states(delet) = [];
            
            states = [states_n states];
            
            bel = cell(M,1);
            for j=1:length(states),
                for i=1:length(states{j}),
                    bel{states{j}(i)} = [bel{states{j}(i)} j];
                end
            end
        end
        
        bel_r = setdiff(bel{r},bel{c});
        bel_c = setdiff(bel{c},bel{r});
        
        for i=1:length(bel_r), % Should c be added to an state of r?
            state_rec = states{bel_r(i)};
            
            j = 1;
            while j<=length(state_rec) && map(state_rec(j),c)>gamma,
                j = j+1;
            end
            
            if j > length(state_rec),
                for k = 1:length(state_rec),
                    if state_rec(k)<c
                        map_v(state_rec(k) + M*(c-1)) = 0;
                    else
                        map_v(c + M*(state_rec(k)-1)) = 0;
                    end
                end
                
                states{bel_r(i)} = [state_rec c];
                bel{c} = [bel{c} bel_r(i)];
                mod2 = 1;
            end
        end
        
        for i=1:length(bel_c), % Should r be added to an state of c?
            state_rec = states{bel_c(i)};
            
            j = 1;
            while j<=length(state_rec) && map(state_rec(j),r)>gamma,
                j = j+1;
            end
            
            if j > length(state_rec),
                for k = 1:length(state_rec),
                    if state_rec(k)<r
                        map_v(state_rec(k) + M*(r-1)) = 0;
                    else
                        map_v(r + M*(state_rec(k)-1)) = 0;
                    end
                end
                
                states{bel_c(i)} = [state_rec r];
                bel{r} = [bel{r} bel_c(i)];
                mod2 = 1;   
            end
        end
        
        if ~sum(sum(mod)) && ~mod2 && isempty(intersect(bel{r},bel{c})),    % non additions: new state
            states{end+1} = [r,c];
            bel{r} = [bel{r} length(states)];
            bel{c} = [bel{c} length(states)];
        end
    end
    
    map_v(ind) = 0;
    ind = find(map_v==max(map_v),1);
end

vs = [];
om = ceil(log10(M));
j=1;
for i=1:length(states),
    if length(states{i})>=siz,
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
for j=1:length(states),
    for i=1:length(states{j}),
        bel{states{j}(i)} = [bel{states{j}(i)} j];
    end
end

