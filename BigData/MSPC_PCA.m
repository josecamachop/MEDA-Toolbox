
function [T2,Q] = MSPC_PCA(Lmodel,pc,list,opt)

% Compute and plot T2 (D-statistic) and the Q statistic using ADICOV.
%
% MSPC_PCA(Lmodel,pcs,list) % minimum call
% MSPC_PCA(Lmodel,pcs,list,opt) % complete call
%
%
% INPUTS:
%
% Lmodels: (Fx1) list of struct Lmodels for monitoring.
%
% pcs: (1x1) number of components.
%
% list: {Fx1} list of strings with the names of the files with x (and
%       optionally y) matrices for testing.
%
% opt: (1x1) options for data plotting.
%       0: no plots
%       1: Statistics (default)
%
%
% OUTPUTS:
%
% T2: T2 statistics (Fx1) 
%
% Q: Q statistics (Fx1) 
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Aug/15.
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

if nargin < 3, error('Error in the number of arguments.'); end;
if nargin < 4, opt = 1; end;

%% Main code

Lc =length(Lmodel);

for j=1:Lc,
    
    load(list{j},'x');
    
    if Lmodel{j}.N,

        if Lmodel{j}.weight~=0,
            xs = applyprep2D(x,Lmodel{j}.av,Lmodel{j}.sc,Lmodel{j}.weight);
        else
            xs = applyprep2D(x,Lmodel{j}.av,Lmodel{j}.sc);
        end
        
        xs = sum(xs,1);
        
        [p,sdT] = Lpca(Lmodel{j});
        t = xs*p(:,1:pc);
        e = xs - t*p(:,1:pc)';
        
        t = t./(ones(size(t,1),1)*sdT);
    
        T2(j) = sum((t).^2);
        Q(j) = sum((e).^2);
    else
        T2(j) = 0;
        Q(j) = 0;
    end
    
end

if opt, 
    plot_vec(T2,[],'D-statistic');
    plot_vec(Q,[],'Q-statistic');
end