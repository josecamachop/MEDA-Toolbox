
function fig_h = plot_vec(vec,olabel,slabel,lcont,opt,pmod,fig_h,leg)

% Bar plot.
%
% plot_vec(vec) % minimum call
% plot_vec(vec,olabel,slabel,lcont,opt,pmod,fig_h,leg) % complete call
%
%
% INPUTS:
%
% vec: (Mx1) vector to plot. 
%
% olabel: (Mx1) name of the x-variables (numbers are used by default), eg.
%   num2str((1:M)')'
%
% slabel: (str) y-variable/statistic plotted (nothing by default)
%
% lcont: (2xM) control limits.
%
% opt: (1x1) options for data plotting.
%       0: bar plot (by default)
%       1: line plot
%
% pmod: (str) character string for line plot.
%
% fig_h: (1x1) handle of figure to plot on (nothing by default)
%
% leg: {Nx1} strings in the legend (nothing by default)
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
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

if nargin < 1, error('Error in the number of arguments.'); end;
s = max(size(vec));
if nargin > 2 && ~isempty(olabel)
    if ndims(olabel)==2 & find(size(olabel)==max(size(olabel)))==2, olabel = olabel'; end
    if size(olabel,1)~=s, error('Error in the dimension of the arguments.'); end;
end
if nargin < 5, opt = 0; end

%% Main code

if ~exist('fig_h')
    fig_h=figure;
else
    hold on
end
if ~opt,
    bar(vec);
else
    if ~exist('pmod'), pmod = 'b'; end;
    plot(vec,pmod);
end
axes_h=get(fig_h,'Children');
if length(axes_h)>1, axes_h=axes_h(1); end;
if exist('olabel')&~isempty(olabel), 
    set(axes_h,'XTick',1:s);
    set(axes_h,'XTickLabel',olabel);
end
if exist('slabel')
    ylabel(slabel,'FontSize',16);
end
set(axes_h,'FontSize',14);
if exist('lcont') && ~isempty(lcont)
    hold on
    b = [0.5:(s+1);0.5:(s+1)];
    a = [lcont(1,:);lcont(1,:)];
    plot(b(2:(end-1))',a(:),'r--','LineWidth',2);
    a = [lcont(2,:);lcont(2,:)];
    plot(b(2:(end-1))',a(:),'r--','LineWidth',2);
end    
axis tight
ax=axis;
axis auto
ax2=axis;
axis([ax(1:2) ax2(3:4)])
if exist('leg'), legend(leg); end
    


        