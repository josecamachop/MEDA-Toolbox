
function [T2,Q] = MSPC_ADICOV(Lmodel,pcs,list,opt)

% Compute and plot T2 (D-statistic) and the Q statistic using ADICOV.
%
% MSPC_ADICOV(Lmodel,pcs,list) % minimum call
% MSPC_ADICOV(Lmodel,pcs,list,opt) % complete call
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
%       1: T2 statistics (default)
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

if nargin < 2, error('Error in the number of arguments.'); end;
if nargin < 3, opt = 1; end;

%% Main code

Lc =length(Lmodel);

for j=1:Lc,
    
    load(list{j},'x');
    
    if Lmodel{j}.N,
        Lm = Lmodel{j};
        Lm.lv = rank(Lm.XX);
        if Lm.type==2,
            [beta,W,P,Q,mat,d] = Lpls(Lm);
        else
            [mat,d] = Lpca(Lm);
        end
        
        onesV = ones(size(x,1),1);

        if Lm.weight~=0,
            xs = applyprep2D(x,Lm.av,Lm.sc,Lm.weight);
        else
            xs = applyprep2D(x,Lm.av,Lm.sc);
        end
        
        ti = ADICOV(Lm.XX,xs,pcs,mat(:,1:pcs),mat(:,1:pcs),onesV);
        ri = ADICOV(Lm.XX,xs,size(mat,2)-pcs,mat(:,pcs+1:end),mat(:,pcs+1:end),onesV);
        
        iti = ADindex(xs,ti,mat(:,1:pcs)*diag(1./sqrt(d(1:pcs))));
        iri = ADindex(xs,ri,mat(:,pcs+1:end));
        T2(j) = iti;
        Q(j) = iri;
    else
        T2(j) = 0;
        Q(j) = 0;
    end
    
end

if opt, 
    plot_vec(T2,[],'D-statistic');
    plot_vec(Q,[],'Q-statistic');
end