function err = spcBootstrap(theta,PCreps,varargin)

% Simulated Power Curves w/ Bootrapping for uncertainty estimation. Takes
% input PCrep from the powercurves routine, and generates an estimate for
% uncertainty according to the input percentile, where it is a numerical
% value between 0 and 1, or the maximum deviation where the input
% percentile is 'max'
%
% Related routines: parglm, asca, apca, parglmVS, parglmMC, createDesign,
% powercurve
%
% err = spcBootstrap(PCreps)   % minimum call
% 
%
% INPUTS
%
% PCreps: and M x N x O array where M are the incremental deltas from SPCS,
% N are the number of factors/interactions to simulate, and O are the
% permutations used in the analysis for each level. A boolean array of 0 or
% 1.
%
% bstrpReps: the number of bootstrap sampling routines to take place for
% uncertainty estimation.
%
% plotResults: true for plots, false for no plots
%
% prcntl: a value between 0 and 1 to specify the uncertainty level. For
% example: prcntl = 0.05 will bound the results between the bottom 0.05 and
% top 0.95 observable combinations of the data for uncertainty estimation.
% Set to 'max' for the maximum and minimum observable differences.
%
%
% OUTPUTS
%
% err: an M x N x 2 array containing the upper and lower bounds in the 1st
% and second slices respectively.
%
% Example of use:
%
% reps = 1;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
% 
% F = createDesign(levels,'Replicates',reps);
% 
% X = [];
% for i = 1:length(levels{1})
%     for j = 1:length(levels{2})
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = randn(reps, vars) + 0.1*repmat(randn(1,vars),reps,1);
%     end
% end
% 
% PCreps = powercurve(X, F,'Model',{[1 2]},'Type',2);
% spcBootstrap(PCreps,100,true,0.05,false,['A','B']);
%
% coded by: Michael Sorochan Armstorng (mdarmstr@ugr.es)
% last modification: 11/Nov/2024
%
% Copyright (C) 2024  University of Granada, Granada
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

sz = size(PCreps);
noDelta = sz(1);
noFctrs = sz(2);
noPerms = sz(3);

if length(varargin) < 1
    bstrpReps = 1000;
else
    bstrpReps = varargin{1};
end
    
if length(varargin) < 2
    plotResults = true;
else
    plotResults = varargin{2};
end
    
if length(varargin) < 3
    prcntl = 0.05;
else
    prcntl = varargin{3};
end
    
if length(varargin) < 4
    colorVec = colormap(lines(noFctrs));
elseif varargin{4} == false && plotResults == true
    colorVec = colormap(lines(noFctrs));
else
    colorVec = varargin{4};
end

if length(varargin) < 5
elseif length(varargin) <= 5 && plotResults == false
    warning('Legend provided, but plotVec is set to false. Setting plotVec to true')
    plotResults = true;
    if colorVec == false
        warning('colorVec set to false. Setting to lines colormap')
        colorVec = colormap(lines(noFctrs));
    end
end

PCuncert = zeros(noDelta,noFctrs,bstrpReps);

for ii = 1:bstrpReps
    for jj = 1:noFctrs
        bstrp = randi(noPerms,noPerms,1);
        PCuncert(:,jj,ii) = sum(PCreps(:,jj,bstrp),3) ./ noPerms;
    end
end

if strcmp(prcntl,'max')
    mins = min(PCuncert,[],3);
    maxs = max(PCuncert,[],3);
    avgs = mean(PCuncert,3);
else
    mins = quantile(PCuncert,prcntl,3);
    maxs = quantile(PCuncert,1-prcntl,3);
    avgs = mean(PCuncert,3);
end

curve1 = maxs;
curve2 = mins;

err = cat(3,curve1,curve2);

x = theta';

if plotResults == true
   hold on;
   for ii = 1:size(PCreps,2)
       patch([x; flip(x)],[curve2(:,ii); flip(curve1(:,ii))],colorVec(ii,:), 'FaceAlpha',0.5, 'EdgeColor','none');
   end
   if length(varargin) < 5
   else
      ax = gca;
      legend(ax, varargin{5})
   end
   plot(x,avgs,'Color','k','HandleVisibility','off')
   xlabel('Effect size (Theta)','FontSize',14)
   ylabel('Power','FontSize',14)
   title('Power Curve','FontSize',14)
   hold off;
end


