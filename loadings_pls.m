
function [P,Q,R] = loadings_pls(cal,y,lvs,prepx,prepy,opt,label,classes)

% Compute and plot loadings in PLS.
%
% loadings_pls(cal,y,lvs) % minimum call
% loadings_pls(cal,y,lvs,prepx,prepy,opt,label,classes) % complete call
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
%       1: loading plot of P (default)
%       2: loading plot of P with empty marks
%       3: loading plot of W*inv(P'*W)
%       4: loading plot of W*inv(P'*W) with empty marks
%
% label: (Mx1) name of the variables (numbers are used by default), eg.
%   num2str((1:M)')'
%
% classes: (Mx1) vector with the assignment of the variables to classes, 
%   numbered from 1 onwards (1 class by default), eg. ones(M,1)
%
%
% OUTPUTS:
%
% P: (MxA) x-block scores.
% Q: (OxA) y-block scores.
% R: (MxA) W*inv(P'*W).
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
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

if nargin < 3, error('Error in the number of arguments.'); end;
s = size(cal);
if s(1) < 1 || s(2) < 1 || ndims(cal)~=2, error('Error in the dimension of the arguments.'); end;
sp = length(lvs);
if sp < 2, error('Error in the dimension of the arguments.'); end;
if nargin < 4, prepx = 2; end;
if nargin < 5, prepy = 2; end; 
if nargin < 6, opt = 1; end;
if nargin < 7 || isempty(label)
    label=num2str((1:s(2))'); 
else
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end
if nargin < 8 || isempty(classes)
    classes = ones(s(2),1); 
else
    if ndims(classes)==2 & find(size(classes)==max(size(classes)))==2, classes = classes'; end
    if size(classes,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end

%% Main code

[calp,m,dt] = preprocess2D(cal,prepx);
yp = preprocess2D(y,prepy);

[beta,W,P,Q] = kernel_pls(calp'*calp,calp'*yp,max(lvs));
R = W*inv(P'*W);

if opt,
    for i=1:length(lvs)-1,
        for j=i+1:length(lvs),
            if opt ==1 || opt ==2,
                plot_scatter([P(:,lvs(i)),P(:,lvs(j))],label,classes,{sprintf('LV %d',lvs(i)),sprintf('LV %d',lvs(j))},opt-1);   
            elseif opt == 3 || opt ==4,
                plot_scatter([R(:,lvs(i)),R(:,lvs(j))],label,classes,{sprintf('LV %d',lvs(i)),sprintf('LV %d',lvs(j))},opt-3);                
            end    
        end
    end
end
        