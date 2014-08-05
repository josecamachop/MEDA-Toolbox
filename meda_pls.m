
function [meda_map,meda_dis] = meda_pls(x,y,lvs,prepx,prepy,thres,opt,label,vars)

% Missing data methods for exploratory data analysis in PLS. The original
% paper is Chemometrics and Intelligent Laboratory Systems 103(1), 2010, pp.
% 8-18. This algorithm follows the suggested computation by Arteaga in his
% technical report "A Note on MEDA", attached to the toolbox, which
% makes use of the covariance matrices.
%
% [meda_map,meda_dis] = meda_pls(x,y,lvs) % minimum call
% [meda_map,meda_dis] = meda_pls(x,y,lvs,prepx,prepy,thres,opt,label,vars) %complete call
%
%
% INPUTS:
%
% x: (NxM) billinear data set of predictor variables
%
% y: (NxO) billinear data set of predicted variables
%
% lvs: (1xA) Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two lvs)
%
% prepx: (1x1) preprocesing of the x-block
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)  
%
% prepy: (1x1) preprocesing of the y-block
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)  
%
% thres: (1x1) threshold for the discretized MEDA matrix (0.1 by default)
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: plot MEDA matrix (default)
%       2: plot discretized MEDA matrix
%       3: plot MEDA matrix seriated 
%       4: plot discretized MEDA matrix seriated
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
% coded by: José Camacho Páez (josecamacho@ugr.es).
% version: 2.0
% last modification: 03/Jul/14.
%
% Copyright (C) 2014  José Camacho Páez
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
s = size(x);
sy = size(y);
if s(1) < 1 || s(1) ~= sy(1) || s(2) < 1 || ndims(x)~=2, error('Error in the dimension of the arguments.'); end;
if nargin < 4, prepx = 2; end;
if nargin < 5, prepy = 2; end;
if nargin < 6, thres=0.1; end; 
if nargin < 7, opt = 1; end;
if nargin < 8 || isempty(label)
    label=num2str((1:s(2))'); 
else
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end
if nargin < 9, vars = 1:s(2); end;

%% Main code

x2 = preprocess2D(x,prepx);
y2 = preprocess2D(y,prepy);

[beta,W,P] = kernel_pls(x2'*x2,x2'*y2,max(lvs));
W2 = W*inv(P'*W);
W2 = W2(:,lvs);
P = P(:,lvs);

[meda_map,meda_dis] = meda(x2'*x2,W2,P,thres);
    
%% Show results

if opt,
    switch opt,
        case 2
            map1 = meda_dis;
            ord = 1:s(2);
        case 3
            [map1, ord] = seriation(meda_map);
        case 4
            [map1, ord] = seriation(meda_dis);
        otherwise
            map1 = meda_map;
            ord = 1:s(2);
    end
    
    varso = [];
    for i=1:length(vars),
        j=find(vars(i)==ord,1);
        varso=[varso j];
    end
    varso = sort(varso);
    
    if ~exist('label')
        label = num2str(ord');
    end
    
    plot_map(map1(varso,varso),label(ord(varso),:));
end

        