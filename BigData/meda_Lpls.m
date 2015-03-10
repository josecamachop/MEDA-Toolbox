function [meda_map,meda_dis] = meda_Lpls(Lmodel,lvs,thres,opt,label,vars)

% Missing data methods for exploratory data analysis in PLS. The original
% paper is Chemometrics and Intelligent Laboratory Systems 103(1), 2010, pp.
% 8-18. This algorithm follows the suggested computation by Arteaga, which
% makes use of the covariance matrices.
%
% [meda_map,meda_dis] = meda_Lpls(Lmodel) % minimum call
% [meda_map,meda_dis] = meda_Lpls(Lmodel,lvs,thres,opt,label,vars) %complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.XY: (MxL) cross-product matrix between the x-block and the
%           y-block.
%       Lmodel.lv: (1x1) number of LVs.
%
% lvs: (1xA) Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two lvs)
%
% thres: (1x1) threshold for the discretized MEDA matrix (0.1 by default)
%
% opt: (struct) options for data plotting
%       plot:
%           0: no plots
%           1: plot MEDA matrix (default)
%           2: plot discretized MEDA matrix
%       seriated:
%           0: no seriated
%           1: seriated (default)
%       discard:
%           0: no discard
%           1: discard 0 variance variables (default)
%
% label: (Mx1) name of the variables (numbers are used by default), eg.
%   num2str((1:M)')'
%
% vars: (1xS) Subset of variables to plot (1:M by default)
%
%
% OUTPUTS:
%
% meda_map: (MxM) MEDA matrix.
%
% meda_dis: (MxM) discretized MEDA matrix.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 11/Mar/15
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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
s = size(Lmodel.XX);
if s(1) ~= s(2) || ndims(Lmodel.XX)~=2, error('Error in the dimension of the arguments.'); end;
if nargin < 2, lvs=1:Lmodel.lv; end;
if nargin < 3, thres=0.1; end; 
if nargin < 4,  
    opt.plot = 1;
    opt.seriated = 1;
    opt.discard = 1;
end;
if isstruct(opt) % backward v1.0 compatibility
    opt2 = opt;
else
    switch opt,
        case 0
            opt2.plot = 0;
            opt2.seriated = 0;
            opt2.discard = 0;
        case 1
            opt2.plot = 1;
            opt2.seriated = 0;
            opt2.discard = 0;
        case 2
            opt2.plot = 2;
            opt2.seriated = 0;
            opt2.discard = 0;       
        case 3
            opt2.plot = 1;
            opt2.seriated = 1;
            opt2.discard = 0; 
        case 4
            opt2.plot = 2;
            opt2.seriated = 1;
            opt2.discard = 0; 
        otherwise
            error('Error in the dimension of the arguments.'); 
    end
end
if nargin < 5 || isempty(label)
    label=num2str((1:s(2))'); 
else
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end
if nargin < 6, vars = 1:s(2); end;

%% Main code

Lmodel.lv = max(lvs);
[beta,W,P,Q,R] = Lpls(Lmodel);

[meda_map,meda_dis] = meda(Lmodel.XX,R,P,thres);

%% Show results

if opt2.plot,
    
    if opt2.plot==1,
        map1 = meda_map;
    else
        map1 = meda_dis;
    end
    
    if opt2.seriated == 1,
        [map1, ord] = seriation(map1);
    else
        ord = 1:s(2);
    end
    
    if opt2.discard == 1,
        Dmap = diag(map1);
        ind = find(Dmap);
        map1 = map1(ind,ind);
        ord = ord(ind); 
    end
    
    if ~exist('label')
        label = num2str(ord');
    end
    
    plot_map(map1,label(ord,:));
    
end
    


        