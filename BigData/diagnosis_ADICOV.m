function [diff,ord,statesp,tp,statesn,tn]= diagnosis_ADICOV(Lmodel,pcs,file,gamma,opt)

% Diagnosis using ADICOV.
%
% Diagnosis_ADICOV(Lmodel,pcs,file) % minimum call
% Diagnosis_ADICOV(Lmodel,pcs,file,opt) % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) Lmodel for monitoring.
%
% pcs: (1xA) Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs)
%
% file: {str} name of the file with x for diagnosis.
%
% gamma: (1x1) correlation threshold to identify groups (0.7 by default)
%
% opt: (1x1) options for data plotting.
%       0: no plots
%       1: Differential map (default)
%
%
% OUTPUTS:
%
% diff: (MxM) Differential map
%
% statesp: (SPx1) groups of positive variables in the map
%
% tp: (NPxSP) scores for each group
%
% statesn: (SNx1) groups of negative variables in the map
%
% tn: (NNxSN) scores for each group
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

if nargin < 3, error('Error in the number of arguments.'); end;
if nargin < 4, gamma = 0.7; end;
if nargin < 5, opt = 1; end;

Lm = Lmodel;
load(file,'x','obs_l');
if Lm.weight~=0,
    xs = applyprep2D(x,Lm.av,Lm.sc,Lm.weight);
else
    xs = applyprep2D(x,Lm.av,Lm.sc);
end

optM.plot = 0;
pcs1 = pcs;
pcs1(find(pcs>rank(Lm.XX))) = [];
meda_map_L = meda_Lpca(Lm,pcs1,0.1,optM);
[meda_map_L,ord] = seriation(meda_map_L);
pcs2 = pcs;
pcs2(find(pcs>rank(xs))) = [];
meda_map_x = meda_pca(xs(:,ord),pcs2,0,0.1,optM);

diff = meda_map_x - meda_map_L;

statesp = [];
statesn = [];
tp = [];
tn = [];

diffp = diff; % positive diagnosis
diffp(diffp<0) = 0;
if find(diffp)
    [bel,statesp] = gia(diffp,gamma,1);  
    if ~isempty(statesp)                              
        [p,t,bel] = gpca(xs(:,ord),statesp,0,1);

        for i=1:length(statesp),
            inds = find(bel==i,1);
            tp(:,i) = t(:,inds);
        end
    end
end

diffn = diff; % negative diagnosis
diffn(diffn>0) = 0;
if find(diffn)
    [bel,statesn] = gia(diffn,gamma,1);    
    if ~isempty(statesn)
        [p,t,bel] = gpca(Lm.centr(:,ord),statesn,0,1);

        for i=1:length(statesn),
            inds = find(bel==i,1);
            tn(:,i) = t(:,inds);
        end
    end
end

if opt, 
    plot_map(diff);
end