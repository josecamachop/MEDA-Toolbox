function err = SPC_bootstrap(theta,PCreps,varargin)

% Simulated Power Curves w/ Bootrapping for uncertainty estimation. Takes
% input PCrep from the powercurves routine, and generates an estimate for
% uncertainty according to the input percentile, where it is a numerical
% value between 0 and 1, or the maximum deviation where the input
% percentile is 'max'
%
% Related routines: parglm, asca, apca, parglmVS, parglmMC, create_design,
% powercurve
%
% err = SPC_bootstrap(PCreps)   % minimum call
% [T, parglmo] = SPC_bootstrap(PCreps,bstrp_reps=1000,plot_results=true,
% ... prcntl=0.05,colorvec=false, legend = ['A','B','A(C)'] % complete call
% 
% INPUTS
%
% PCreps: and M x N x O array where M are the incremental deltas from SPCS,
% N are the number of factors/interactions to simulate, and O are the
% permutations used in the analysis for each level. A boolean array of 0 or
% 1.
%
% bstrp_reps: the number of bootstrap sampling routines to take place for
% uncertainty estimation.
%
% plot_results: true for plots, false for no plots
%
% prcntl: a value between 0 and 1 to specify the uncertainty level. For
% example: prcntl = 0.05 will bound the results between the bottom 0.05 and
% top 0.95 observable combinations of the data for uncertainty estimation.
% Set to 'max' for the maximum and minimum observable differences.
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
% F = create_design(levels,'Replicates',reps);
% 
% X = [];
% for i = 1:length(levels{1})
%     for j = 1:length(levels{2})
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = randn(reps, vars) + 0.1*repmat(randn(1,vars),reps,1);
%     end
% end
% 
% PCreps = powercurve(X, F,'Model',{[1 2]},'Type',2);
% SPC_bootstrap(PCreps,100,true,0.05,false,['A','B']);
%
% coded by: Michael Sorochan Armstorng (mdarmstr@ugr.es)
% last modification: 23/Apr/2024
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
no_delta = sz(1);
no_fctrs = sz(2);
no_perms = sz(3);

if length(varargin) < 1
    bstrp_reps = 1000;
else
    bstrp_reps = varargin{1};
end
    
if length(varargin) < 2
    plot_results = true;
else
    plot_results = varargin{2};
end
    
if length(varargin) < 3
    prcntl = 0.05;
else
    prcntl = varargin{3};
end
    
if length(varargin) < 4
    color_vec = colormap(lines(no_fctrs));
elseif varargin{4} == false && plot_results == true
    color_vec = colormap(lines(no_fctrs));
else
    color_vec = varargin{4};
end

if length(varargin) < 5
elseif length(varargin) <= 5 && plot_results == false
    warning('Legend provided, but plot_vec is set to false. Setting plot_vec to true')
    plot_results = true;
    if color_vec == false
        warning('color_vec set to false. Setting to lines colormap')
        color_vec = colormap(lines(no_fctrs));
    end
end

PC_uncert = zeros(no_delta,no_fctrs,bstrp_reps);

for ii = 1:bstrp_reps
    for jj = 1:no_fctrs
        bstrp = randi(no_perms,no_perms,1);
        PC_uncert(:,jj,ii) = sum(PCreps(:,jj,bstrp),3) ./ no_perms;
    end
end

if strcmp(prcntl,'max')
    mins = min(PC_uncert,[],3);
    maxs = max(PC_uncert,[],3);
    avgs = mean(PC_uncert,3);
else
    mins = quantile(PC_uncert,prcntl,3);
    maxs = quantile(PC_uncert,1-prcntl,3);
    avgs = mean(PC_uncert,3);
end

curve1 = maxs;
curve2 = mins;

err = cat(3,curve1,curve2);

x = theta';

if plot_results == true
   hold on;
   for ii = 1:size(PCreps,2)
       patch([x; flip(x)],[curve2(:,ii); flip(curve1(:,ii))],color_vec(ii,:), 'FaceAlpha',0.5, 'EdgeColor','none');
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


