
function [T,TT] = scores_pls(cal,y,lvs,test,prepx,prepy,opt,label,classes)

% Compute and plot scores in PLS.
%
% scores_pls(cal,y,lvs) % minimum call
% scores_pls(cal,y,lvs,test,prepx,prepy,opt,label,classes) % complete call
%
% INPUTS:
%
% cal: (LxM) billinear data set for model fitting
%
% y: (LxO) billinear data set of predicted variables
%
% lvs: (1xA) Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two lvs)
%
% test: (NxM) billinear data set for test. These data are preprocessed in
%   the same way than calibration data.
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
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: score plot (default)
%       2: score plot with empty marks
%
% label: ((L+N)x1) name of the variables (numbers are used by default), eg.
%   num2str((1:L+N))')', use ' ' to avoid labels.
%
% classes: ((L+N)x1) vector with the assignment of the variables to classes, 
%   numbered from 1 onwards (1 class by default), eg. ones(L+N),1)
%
%
% OUTPUTS:
%
% T: (LxA) calibration scores.
%
% TT: (NxA) test scores.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
%
% Copyright (C) 2014  University of Granada, Granada
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
if nargin < 4, x = cal; else x = [cal;test]; end;
s = size(x);
if s(1) < 1 || s(2) < 1 || ndims(x)~=2, error('Error in the dimension of the arguments.'); end;
sp = length(lvs);
if sp < 2, error('Error in the dimension of the arguments.'); end;
if nargin < 5, prepx = 2; end;
if nargin < 6, prepy = 2; end; 
if nargin < 7, opt = 1; end;
if nargin < 8 || isempty(label)
    label=num2str((1:s(1))'); 
elseif ~isequal(label,' '),
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(1), error('Error in the dimension of the arguments.'); end;
end
if nargin < 9 || isempty(classes)
    classes = [ones(1,size(cal,1)) 2*ones(1,size(test,1))];
else
    if ndims(classes)==2 & find(size(classes)==max(size(classes)))==2, classes = classes'; end
    if size(classes,1)~=s(1), error('Error in the dimension of the arguments.'); end;
end

%% Main code

[calp,m,dt] = preprocess2D(cal,prepx);
yp = preprocess2D(y,prepy);

[beta,W,P] = kernel_pls(calp'*calp,calp'*yp,max(lvs));
W2 = W*inv(P'*W);
T = calp*W2;

if exist('test')&~isempty(test),
    testp = (test - ones(size(test,1),1)*m)./(ones(size(test,1),1)*dt);
    TT = testp*W2;
else
    TT = [];
end

if opt,
    T = [T;TT];
    for i=1:length(lvs)-1,
        for j=i+1:length(lvs),
            plot_scatter([T(:,lvs(i)),T(:,lvs(j))],label,classes,{sprintf('LV %d',lvs(i)),sprintf('LV %d',lvs(j))},opt-1);
        end      
    end
end
        