
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

% Set default values
assert (nargin >= 1, 'Error: Missing arguments.');
if nargin < 2 || isempty(olabel), olabel = []; end;
if nargin < 3 || isempty(slabel), slabel = ''; end;
if nargin < 4 || isempty(lcont),  lcont = []; end;
if nargin < 5 || isempty(opt),    opt = 0; end;
if nargin < 6 || isempty(pmod),   pmod = ''; end;
if nargin < 7 || isempty(fig_h),  fig_h = []; end;
if nargin < 8 || isempty(leg),    leg = []; end;

% Validate parameters
if size(vec,1) == 1, vec = vec'; end;
assert (size(vec,2) == 1, 'Dimension Error: vec must be n-by-1.');
N = size(vec, 1);

if ~isempty(olabel)
    if size(olabel,2) > size(olabel,1), olabel = olabel'; end;
    if ischar(olabel), olabel = cellstr(olabel); end;
    %if size(olabel,1) == 1, olabel = olabel'; end;
    assert (isequal(size(olabel), [N 1]), 'Dimension Error: olabel must be n-by-1.');
end
if ~isempty(lcont)
    assert (isequal(size(lcont), [2 N]), 'Dimension Error: lcont must be 2-by-n');
end
if ~isempty(fig_h)
    assert (isscalar(fig_h));
end


%% Main code

% Create figure window
if isempty(fig_h)
    fig_h = figure;
else
    hold on
end

%color1 = [182,182,92]./255;
%color2 = [92, 157, 182]./255;
%color3 = [233,72,9]./255;

% Plot bar graph
defaultcolor = [0, 154, 179]./255;   % light blue
if ~opt,
    bar(vec, 'FaceColor', defaultcolor, 'EdgeColor', 'w');
else
    if pmod,
        plot(vec, pmod, 'LineWidth', 3);
    else
        plot(vec, 'LineWidth', 3, 'Color', defaultcolor);
    end
end

% Plot control limits
if ~isempty(lcont)
    hold on
    b = [0.5:(N+1);0.5:(N+1)];
    a = [lcont(1,:);lcont(1,:)];
    plot(b(2:(end-1))',a(:),'r--','LineWidth',2);
    a = [lcont(2,:);lcont(2,:)];
    plot(b(2:(end-1))',a(:),'r--','LineWidth',2);
end    

% Get axes handler
axes_h = get(fig_h,'Children');
if length(axes_h)>1, axes_h = axes_h(1); end;

% Set ticks and labels
label_size = max(min(14,round(300/length(olabel))), 9);
set(axes_h, 'FontSize', label_size);
if ~isempty(olabel)
    set(axes_h,'XTick',1:N);
    label_length = max(cellfun('length', olabel));
    if label_length > 2
        % rotate labels
        set(axes_h,'XTickLabel',[]);
        b = get(axes_h, 'XTick');
        c = get(axes_h, 'YTick');
        text(b, repmat(c(1)-.1*(c(2)-c(1)),length(b),1), olabel,...
            'FontSize',label_size,'HorizontalAlignment','right','rotation',90);
    else
        % dont rotate labels
        set(axes_h,'XTickLabel', olabel);
    end
end
if ~isempty(slabel)
    ylabel(slabel, 'FontSize', 16);
end

% Set legend
if ~isempty(leg), legend(leg); end;

% Set axis
axis tight
ax = axis;
axis auto
ax2 = axis;
axis([ax(1:2) ax2(3:4)])

       