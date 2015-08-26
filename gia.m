function [bel,states] = gia(map,gamma,siz)

% Group identification algorithm (gia) to generate groups of variables from 
% a correlation map of variables x variables.
%
% [bel,states] = gia(map,gamma,siz)   % minimum call
%
% [bel,states] = gia(map,gamma,siz)   % complete call
%
%
% INPUTS:
%
% map: (MxM) correlation matrix. 
%
% gamma: (1x1) correlation threshold to identify groups (0.7 by default)
%
% siz: (1x1) Minimum size of groups (2 by default)
%
%
% OUTPUTS:
%
% bel: {Mx1} Cell with the list of groups each variable belongs to.
%
% states: {Sx1} Cell with the groups of variables.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/Aug/15.
%
% Copyright (C) 2015  University of Granada, Granada
% Copyright (C) 2015  Jose Camacho Paez
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

%% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
if nargin < 2, gamma=0.7; end; 
if nargin < 3, siz=2; end; 

%% Main code

l = length(map);
map = abs(map);

states = {};
bel = cell(l,1);

columns = ones(l,1)*(1:l);
rows = columns';

map_v = map - tril(map) + diag(diag(map)); % triangular inferior part is set to zero
map_v = map_v(:);
columns_v = columns(:);
rows_v = rows(:);

ind = find(map_v==max(map_v),1);
while map_v(ind)>gamma,
    
    r = rows_v(ind); % r is always minor than c
    c = columns_v(ind);
    
    if c==r,
        if isempty(bel{r})
            states{end+1} = r;
            bel{r} = [bel{r} length(states)];
        end
    else
        int = intersect(bel{r},bel{c});
        
        bel_r = setdiff(bel{r},int);
        bel_c = setdiff(bel{c},int);
        
        mod = zeros(length(bel_r),length(bel_c));
        mod2 = 0;
        states_n = {};
        
        for i = 1:length(bel_r), % Should two states be combined?
            for j = 1:length(bel_c),
                que = (map(states{bel_r(i)},states{bel_c(j)}) > gamma);
                if isempty(find(~que)),
                    mod(i,j) = 1;
                    states_n{end+1} = [states{bel_r(i)} states{bel_c(j)}]; % ERROR: REVISAR
                end
            end
        end
        
        if sum(sum(mod))
            ind_r = find(sum(mod,2));
            ind_c = find(sum(mod,1));
            
            delet = unique([bel_r(ind_r) bel_c(ind_c)]);
            
            states(delet) = [];
            
            states = [states_n states];
            
            bel = cell(l,1);
            for i=1:l,
                for j=1:length(states),
                    if ismember(i,states{j}),
                        bel{i} = [bel{i} j];
                    end
                end
            end
        end
        
        int = intersect(bel{r},bel{c});
        
        bel_r = setdiff(bel{r},int);
        bel_c = setdiff(bel{c},int);
        
        for i=1:length(bel_r), % Should c be added to an state of r?
            state_rec = states{bel_r(i)};
            
            j = 1;
            while j<=length(state_rec) && map(state_rec(j),c)>gamma,
                j = j+1;
            end
            
            if j > length(state_rec),
                for k = 1:length(state_rec),
                    if state_rec(k)<c
                        map_v(state_rec(k) + l*(c-1)) = 0;
                    else
                        map_v(c + l*(state_rec(k)-1)) = 0;
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
                        map_v(state_rec(k) + l*r) = 0;
                    else
                        map_v(r + l*state_rec(k)) = 0;
                    end
                end
                
                states{bel_c(i)} = [state_rec r];
                bel{r} = [bel{r} bel_c(i)];
                mod2 = 1;
            end
        end
        
        if ~sum(sum(mod)) && ~mod2 && isempty(int),    % non additions: new state
            states{end+1} = [r,c];
            bel{r} = [bel{r} length(states)];
            bel{c} = [bel{c} length(states)];
        end
    end
    
    map_v(ind) = 0;
    ind = find(map_v==max(map_v),1);
end

vs = [];
om = ceil(log10(l));
j=1;
for i=1:length(states),
    if length(states{i})>=siz,
        stateso{j} = sort(states{i});
        vs(j) = sum(states{i}.*(10.^(-1*(1:om:om*length(states{i})))));  
        j = j+1;
    end
end    

[kk,ind] = sort(vs);
states = stateso(ind);

bel = cell(l,1);
for i=1:l,
    for j=1:length(states),
        if ismember(i,states{j}),
            bel{i} = [bel{i} j];
        end
    end
end

