
function T2 = leverage_pca(cal,pcs,test,prep,opt,label)

% Compute and plot leverages in PCA.
%
% leverage_pca(cal,pcs) % minimum call
% leverage_pca(cal,pcs,test,prep,opt,label) % complete call
%
% INPUTS:
%
% cal: (LxM) billinear data set for model fitting.
%
% pcs: (1xA) Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs)
%
% test: (NxM) data set with test observations. These data are preprocessed 
%    in the same way than calibration data and are only used in opts
%    0,1,3,4.
%
% prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)
%
% opt: (1x1) options for data plotting.
%       0: no plots
%       1: D-statistic (leverage) in the observations (default)
%       2: Captured sum-of-squares in the variables 
%       3: D-statistic (leverage) in the observations with control limits
%       4: D-statistic (leverage) in the observations with control limits in
%           Phase I
%
% label: name of the observations (opt 1, dimension ((L+N)x1) or 
%   variables (opt 2, dimension (Mx1)) (numbers are used by default), eg.
%   num2str((1:L+N))')'
%
%
% OUTPUTS:
%
% T2: D-statisic (opt 0, 1, 3 or 4, dimension {1x(L+N)}) or Sum of squares
% captured per variable (opt 2, dimension {1x(M)})
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 09/Nov/15.
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

if nargin < 2, error('Error in the number of arguments.'); end;
if nargin < 3, x = cal; else x = [cal;test]; end;
s = size(x);
if s(1) < 1 || s(2) < 1 || ndims(x)~=2, error('Error in the dimension of the arguments.'); end;
sp = length(pcs);
if nargin < 4, prep = 2; end;
if nargin < 5, opt = 1; end;
if nargin < 6, label = []; end

%% Main code

T2 = [];

if ~isempty(pcs) && ~isempty(find(pcs))
    
    [calp,m,dt] = preprocess2D(cal,prep);
    [P,T] = pca_pp(calp,max(pcs));
    P = P(:,pcs);
    T = T(:,pcs);
    
    if exist('test')&~isempty(test),
        testp = (test - ones(size(test,1),1)*m)./(ones(size(test,1),1)*dt);
        TT = testp*P;
    else
        TT = [];
    end


    switch opt,
        case 2
            T2 = sum((T*P').^2,1)';
            if ~isempty(test)
                T2 = [T2 sum((TT*P').^2,1)'];
            end
            if iscell(label),
                label = {label{:} label{:}};
            else
                label = [label label];
            end
        otherwise
            [Ts,notused,dtT2] = preprocess2D(T,2);
            T2c = diag(Ts*Ts');
            if ~isempty(test)
                T2 = diag(TT*diag(1./(dtT2.^2))*TT');
            else
                T2 = T2c;
            end
    end;
    
    if opt,    
        if opt<2,
            plot_vec(T2,label,'D-statistic');
        elseif opt<3
            plot_vec(T2,label,'Captured Sum-of-squares');
        elseif opt<4,
            plot_vec(T2,label,'D-statistic',(ones(size(T2,1),1)*[hot_lim(length(pcs),length(T2c),0.05) hot_lim(length(pcs),length(T2c),0.01)])');
        else
            plot_vec(T2,label,'D-statistic',(ones(size(T2,1),1)*[hot_lim(length(pcs),length(T2c),0.05,1) hot_lim(length(pcs),length(T2c),0.01,1)])');
        end
    end
end

